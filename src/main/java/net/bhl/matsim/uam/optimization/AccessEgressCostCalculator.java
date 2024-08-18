package net.bhl.matsim.uam.optimization;

import ch.sbb.matsim.routing.pt.raptor.*;
import net.bhl.matsim.uam.analysis.traveltimes.utils.ThreadCounter;
import net.bhl.matsim.uam.optimization.utils.ScenarioSpecific;
import net.bhl.matsim.uam.optimization.utils.TripItemForOptimization;
import org.matsim.api.core.v01.Scenario;
import org.matsim.api.core.v01.TransportMode;
import org.matsim.api.core.v01.network.Network;
import org.matsim.core.config.Config;
import org.matsim.core.controler.AbstractModule;
import org.matsim.core.controler.Injector;
import org.matsim.core.network.NetworkUtils;
import org.matsim.core.network.algorithms.TransportModeNetworkFilter;
import org.matsim.core.router.AStarLandmarksFactory;
import org.matsim.core.router.RoutingModule;
import org.matsim.core.router.TeleportationRoutingModule;
import org.matsim.core.router.util.*;
import org.matsim.core.scenario.ScenarioUtils;
import org.matsim.core.trafficmonitoring.TravelTimeCalculator;

import java.io.IOException;
import java.util.*;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.matsim.pt.router.TransitRouter;
import org.apache.log4j.Logger;
public class AccessEgressCostCalculator {

    public AccessEgressCostCalculator(List<TripItemForOptimization> tripItemForOptimizations, List<Vertiport> vertiportsCandidates, Config config, ScenarioSpecific scenarioSpecific) {
        this.tripItemForOptimizations = tripItemForOptimizations;
        this.vertiportsCandidates = vertiportsCandidates;
        this.config = config;
        this.scenarioSpecific = scenarioSpecific;
    }
    // This class is used to calculate the access and egress time and distance for each vertiport candidate of each trip
    // Input in args: trip file, config file, vertiport candidate file, output file
    // Format of trip file csv: tripID, personID, originX, originY, destinationX, destinationY, departureTime (in seconds), pTravelTime, pTripLength, pInvehicleTime, pWaitingTime, carTravelTime, carTripLength, tripPurpose, carTravelCost, pTravelCost, carUtility, pUtility, UAMUtilityFix (only related to the traveller itself, including income, age,...), carGeneralizedCost, pGeneralizedCost, Income (€/year) # All times are in seconds, all distances are in meters, all costs are in €
    // Format of vertiport candidate file csv: vertiportID, coordX, coordY, constructionCost (optional)
    // Output file should be in format: .dat
    // Warning: if you make any changes to the TripItemForOptimization class or Vertiport class, you need to run this class again to update the serialized file, even if you just add some space or empty lines.
    private List<TripItemForOptimization> tripItemForOptimizations;
    private List<Vertiport> vertiportsCandidates;
    private HashMap<Integer,Vertiport> vertiportCandidatesMap;
    private Config config;
    private ScenarioSpecific scenarioSpecific;
    private static final int processes = Runtime.getRuntime().availableProcessors();
    private static final Logger log = Logger.getLogger(PreCalculateAccessEgressCost.class);
    private static ArrayBlockingQueue<LeastCostPathCalculator> carRouters = new ArrayBlockingQueue<>(processes);
    private static ArrayBlockingQueue<TransitRouter> ptRouters = new ArrayBlockingQueue<>(processes);
    public void calculateAccessEgressCost() throws IOException, InterruptedException {

        //Create scenario
        Scenario scenario = ScenarioUtils.createScenario(config);
        ScenarioUtils.loadScenario(scenario);
        Network network = scenario.getNetwork();
        RaptorStaticConfig raptorStaticConfig = RaptorUtils.createStaticConfig(config);
        SwissRailRaptorData data = SwissRailRaptorData.create(scenario.getTransitSchedule(), raptorStaticConfig,
                network);

        // CREATE CAR NETWORK
        TransportModeNetworkFilter filter = new TransportModeNetworkFilter(network);
        Network networkCar = NetworkUtils.createNetwork();
        Set<String> modesCar = new HashSet<>();
        modesCar.add(TransportMode.car);
        filter.filter(networkCar, modesCar);

        // LEAST COST PATH CALCULATOR
        TravelTimeCalculator.Builder builder = new TravelTimeCalculator.Builder(network);
        builder.configure(config.travelTimeCalculator());
        TravelTimeCalculator ttc = builder.build();
        TravelTime travelTime = ttc.getLinkTravelTimes();
        TravelDisutility travelDisutility = TravelDisutilityUtils
                .createFreespeedTravelTimeAndDisutility(config.planCalcScore());

        com.google.inject.Injector injector = Injector.createInjector(config, new AbstractModule() {
            @Override
            public void install() {
                bind(LeastCostPathCalculatorFactory.class).to(AStarLandmarksFactory.class);
            }
        });
        LeastCostPathCalculatorFactory pathCalculatorFactory = injector
                .getInstance(LeastCostPathCalculatorFactory.class); // AStarLandmarksFactory

        // Provide routers
        for (int i = 0; i < processes; i++) {
            carRouters.add(pathCalculatorFactory.createPathCalculator(networkCar, travelDisutility, travelTime));
            Map<String, RoutingModule> router = new HashMap<>();
            router.put(TransportMode.pt, new TeleportationRoutingModule(TransportMode.pt,
                    scenario, 0, 1.5));
            ptRouters.add(new SwissRailRaptor(data, new DefaultRaptorParametersForPerson(config),
                    new LeastCostRaptorRouteSelector(),
                    new DefaultRaptorStopFinder(null, new DefaultRaptorIntermodalAccessEgress(), router)));
        }
        ThreadCounter threadCounter = new ThreadCounter();
        ExecutorService es = Executors.newFixedThreadPool(processes);
        List<TripItemForOptimization> uamEnabledTrips = new ArrayList<>();
        for (int i = 0; i < tripItemForOptimizations.size(); i++) {
            TripItemForOptimization currentTrip = tripItemForOptimizations.get(i);
            // Find the neighbouring vertiports for the origin and destination
            VertiportCollector vertiportCollector = new VertiportCollector(currentTrip,networkCar,network,vertiportsCandidates,threadCounter,carRouters, ptRouters,scenarioSpecific);
            vertiportCollector.neighbourVertiportCandidateIdentifier();
            // provide information after each 100 trips
            if (i % 10000 == 0) {
                log.info("Neighbour Vertiport Candidate Indentification: Trip " + i + " is processed.");
            }
            if (!currentTrip.originNeighborVertiportCandidates.isEmpty() && !currentTrip.destinationNeighborVertiportCandidates.isEmpty()) {
                if (currentTrip.originNeighborVertiportCandidates.size() > 1 || currentTrip.destinationNeighborVertiportCandidates.size() > 1) {
                    uamEnabledTrips.add(currentTrip);
                } else if (!(currentTrip.originNeighborVertiportCandidates.get(0).equals(currentTrip.destinationNeighborVertiportCandidates.get(0)))) {
                    uamEnabledTrips.add(currentTrip);
                }
            }
        }
        // Show information in logfile
        log.info("The number of trips that can be served by UAM is: " + uamEnabledTrips.size());
        log.info("Calculating travel times...");


        for (int i = 0; i < uamEnabledTrips.size(); i++) {
            TripItemForOptimization currentTrip = uamEnabledTrips.get(i);

            // For each trip, calculate the access and egress time and distance to all the vertiport candidates

            // Find the neighbouring vertiports for the origin and destination

            // provide information after each 1000 trips
            if (uamEnabledTrips.size()<100 || (i+1) % 10000 == 0) {
                log.info("Calculation completion: " + (i+1) + "/" + uamEnabledTrips.size() + " ("
                        + String.format("%.0f", (double) (i+1) / uamEnabledTrips.size() * 100) + "%).");
            }
            while (threadCounter.getProcesses() >= processes - 1)
                Thread.sleep(200);
            VertiportCollector vertiportCollector = new VertiportCollector(currentTrip,networkCar,network,vertiportsCandidates,threadCounter,carRouters, ptRouters,scenarioSpecific);
            es.execute(vertiportCollector);
        }
        es.shutdown();
        // Make sure that the file is not written before all threads are finished
        while (!es.isTerminated())
            Thread.sleep(200);

    }

}

