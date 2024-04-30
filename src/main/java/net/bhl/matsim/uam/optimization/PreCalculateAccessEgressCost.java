package net.bhl.matsim.uam.optimization;
import net.bhl.matsim.uam.analysis.traveltimes.utils.ThreadCounter;
import net.bhl.matsim.uam.analysis.traveltimes.utils.TripItem;
import net.bhl.matsim.uam.analysis.traveltimes.utils.TripItemReader;
import net.bhl.matsim.uam.config.UAMConfigGroup;
import net.bhl.matsim.uam.optimization.Vertiport;
import net.bhl.matsim.uam.optimization.VertiportCollector;
import net.bhl.matsim.uam.optimization.VertiportReader;
import org.matsim.api.core.v01.Scenario;
import org.matsim.api.core.v01.TransportMode;
import org.matsim.api.core.v01.network.Link;
import org.matsim.api.core.v01.network.Network;
import org.matsim.core.config.Config;
import org.matsim.core.config.ConfigUtils;
import org.matsim.core.controler.AbstractModule;
import org.matsim.core.controler.Injector;
import org.matsim.core.network.NetworkUtils;
import org.matsim.core.network.algorithms.TransportModeNetworkFilter;
import org.matsim.core.router.AStarLandmarksFactory;
import org.matsim.core.router.util.*;
import org.matsim.core.router.util.LeastCostPathCalculator.Path;
import org.matsim.core.scenario.ScenarioUtils;
import org.matsim.core.trafficmonitoring.TravelTimeCalculator;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.*;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
public class PreCalculateAccessEgressCost {
    // This class is used to calculate the access and egress time and distance for each vertiport candidate of each trip
    // Input in args: trip file, config file, vertiport candidate file, output file
    // Format of trip file csv: tripID, personID, originX, originY, destinationX, destinationY, departureTime (in seconds), pTravelTime, pTripLength, pInvehicleTime, pWaitingTime, carTravelTime, carTripLength, tripPurpose, carTravelCost, pTravelCost, carUtility, pUtility, UAMUtilityFix (only related to the traveller itself, including income, age,...), carGeneralizedCost, pGeneralizedCost, Income (€/year) # All times are in seconds, all distances are in meters, all costs are in €
    // Format of vertiport candidate file csv: vertiportID, coordX, coordY, constructionCost (optional)
    // Output file should be in format: .dat
    // Warning: if you make any changes to the TripItem class or Vertiport class, you need to run this class again to update the serialized file, even if you just add some space or empty lines.
    private static String configPath;
    private static String tripFile;
    private static String vertiportCandidateFile;
    private static String outputTripFile;
    private static final int processes = Runtime.getRuntime().availableProcessors();

    private static ArrayBlockingQueue<LeastCostPathCalculator> carRouters = new ArrayBlockingQueue<>(processes);

    public static void main(String[] args) throws IOException {
        // Provide the file via program arguments
        if (args.length > 0) {
            tripFile = args[0];
            configPath=args[1];
            vertiportCandidateFile=args[2];
            outputTripFile=args[3];


        }
        // Read all vertiport candidates
        VertiportReader vertiportReader = new VertiportReader();
        List<Vertiport> vertiportsCandidates = vertiportReader.getVertiports(vertiportCandidateFile);
        List<Vertiport> vertiports = new ArrayList<>();

// READ CONFIG
        Config config = ConfigUtils.loadConfig(configPath, new UAMConfigGroup());
        //Create scenario
        Scenario scenario = ScenarioUtils.createScenario(config);
        ScenarioUtils.loadScenario(scenario);
        Network network = scenario.getNetwork();

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
        }
        LeastCostPathCalculator pathCalculator = pathCalculatorFactory.createPathCalculator(networkCar, travelDisutility, travelTime);
        // Read the trip file and store in a list
        TripItemReader tripItemReader = new TripItemReader();
        List<TripItem> tripItems = tripItemReader.getTripItems(tripFile);
        List<TripItem> uamEnabledTrips = new ArrayList<>();
        for (int i=0;i<tripItems.size();i++) {
            TripItem currentTrip = tripItems.get(i);
            // Find the neighbouring vertiports for the origin and destination
            VertiportCollector vertiportCollector = new VertiportCollector(currentTrip, networkCar, vertiportsCandidates);
            vertiportCollector.neighbourVertiportCandidateIdentifier();
            // provide information after each 100 trips
            if (i % 100 == 0) {
                System.out.println("Neighbour Vertiport Candidate Indentification: Trip " + i + " is processed.");
            }
            if(currentTrip.originNeighborVertiportCandidates.size()>0 && currentTrip.destinationNeighborVertiportCandidates.size()>0){
                if(currentTrip.originNeighborVertiportCandidates.size()>1 || currentTrip.destinationNeighborVertiportCandidates.size()>1){
                    uamEnabledTrips.add(currentTrip);
                } else if (!(currentTrip.originNeighborVertiportCandidates.get(0).equals(currentTrip.destinationNeighborVertiportCandidates.get(0)))){
                    uamEnabledTrips.add(currentTrip);
                }
            }

        }
        // Show information in logfile
        System.out.println("The number of trips that can be served by UAM is: "+uamEnabledTrips.size());
// Create the LeastCostPathCalculator



        for (int i=0;i<uamEnabledTrips.size();i++) {
            TripItem currentTrip = uamEnabledTrips.get(i);

            // For each trip, calculate the access and egress time and distance to all the vertiport candidates

            // Find the neighbouring vertiports for the origin and destination
            VertiportCollector vertiportCollector = new VertiportCollector(currentTrip, networkCar, vertiportsCandidates);
            vertiportCollector.neighbourVertiportCandidateTimeAndDistanceCalculator();
            // provide information after each 1000 trips
            if (i % 100 == 0) {
                System.out.println("Trip " + i + " is processed.");
            }
        }

        try (FileOutputStream fileOut = new FileOutputStream(outputTripFile);
             ObjectOutputStream out = new ObjectOutputStream(fileOut)) {
            out.writeObject(tripItems);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
        // Save the access and egress time and distance for each vertiport candidate of each trip and store in a file

    public static Path estimatePath(Link from, Link to, double departureTime, Network carNetwork,
                                     LeastCostPathCalculator pathCalculator) {
        if (carNetwork.getLinks().get(from.getId()) == null)
            from = NetworkUtils.getNearestLinkExactly(carNetwork, from.getCoord());

        if (carNetwork.getLinks().get(to.getId()) == null)
            to = NetworkUtils.getNearestLinkExactly(carNetwork, to.getCoord());

        return pathCalculator.calcLeastCostPath(from.getFromNode(), to.getToNode(), departureTime, null, null);
    }

    static class CarTravelTimeCalculator {

        private TripItem trip;
        private ThreadCounter threadCounter;
        private Network networkCar;
        private LeastCostPathCalculator plcpccar;

        CarTravelTimeCalculator(ThreadCounter threadCounter, Network network, TripItem trip) {
            this.threadCounter = threadCounter;
            this.networkCar = network;
            this.trip = trip;
        }

        public Map<String, Double> calculateTravelInfo() {
            Map<String, Double> travelInfo = new HashMap<>();

            try {
                plcpccar = carRouters.take();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }

            Link from = NetworkUtils.getNearestLink(networkCar, trip.origin);
            Link to = NetworkUtils.getNearestLink(networkCar, trip.destination);

            try {
                Path path = estimatePath(from, to, trip.departureTime, networkCar, plcpccar);
                double distance = 0;
                for (Link link : path.links) {
                    distance += link.getLength();
                }

                travelInfo.put("distance", distance);
                travelInfo.put("travelTime", path.travelTime);
            } catch (NullPointerException e) {
                // Do nothing; failed trip will show as null in results.
            }

            try {
                carRouters.put(plcpccar);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }

            return travelInfo;
        }
    }
}
