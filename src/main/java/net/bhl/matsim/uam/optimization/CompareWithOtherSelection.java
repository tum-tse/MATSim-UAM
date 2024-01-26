package net.bhl.matsim.uam.optimization;

import com.opencsv.CSVWriter;
import net.bhl.matsim.uam.analysis.traveltimes.utils.ThreadCounter;
import net.bhl.matsim.uam.analysis.traveltimes.utils.TripItem;
import net.bhl.matsim.uam.analysis.traveltimes.utils.TripItemReader;
import net.bhl.matsim.uam.config.UAMConfigGroup;
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
import org.matsim.core.scenario.ScenarioUtils;
import org.matsim.core.trafficmonitoring.TravelTimeCalculator;

import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ArrayBlockingQueue;

import static net.bhl.matsim.uam.optimization.VertiportOptimizerGreedyForwards.*;

public class CompareWithOtherSelection {
    private static String configPath;
    private static String tripFile;
    private static String vertiportCandidateFile;
    private static final int processes = Runtime.getRuntime().availableProcessors();
    private static String servedTripsIDFile;

    private static ArrayBlockingQueue<LeastCostPathCalculator> carRouters = new ArrayBlockingQueue<>(processes);

    public static void main(String[] args) throws IOException {
        // Provide the file via program arguments
        if (args.length > 0) {
            tripFile = args[0];
            configPath=args[1];
            vertiportCandidateFile=args[2];
            servedTripsIDFile=args[3];



        }
        // Read all vertiport candidates
        VertiportReader vertiportReader = new VertiportReader();
        List<Vertiport> vertiportsCandidatesOBUAM = vertiportReader.getVertiports(vertiportCandidateFile);


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
            VertiportCollector vertiportCollector = new VertiportCollector(currentTrip, networkCar, vertiportsCandidatesOBUAM);
            vertiportCollector.neighbourVertiportCandidateIdentifier();
            // provide information after each 100 trips
            if (i % 10000 == 0) {
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
        System.out.println("The number of trips that can be served by UAM in current Vertiport Selection is: "+uamEnabledTrips.size());
// Create the LeastCostPathCalculator
// write out the trips ID that can be served by UAM
        List<Integer> uamAvailableTripsID = new ArrayList<>();
        for (TripItem tripItem : uamEnabledTrips) {
            uamAvailableTripsID.add(tripItem.tripID);
        }
        FileWriter fileWriter = new FileWriter(servedTripsIDFile);
        CSVWriter csvWriter = new CSVWriter(fileWriter);
        csvWriter.writeNext(new String[]{"tripID"});
        for (Integer tripID : uamAvailableTripsID) {
            csvWriter.writeNext(new String[]{tripID.toString()});
        }
        csvWriter.close();
        fileWriter.close();



        for (int i=0;i<uamEnabledTrips.size();i++) {
            TripItem currentTrip = uamEnabledTrips.get(i);

            // For each trip, calculate the access and egress time and distance to all the vertiport candidates

            // Find the neighbouring vertiports for the origin and destination
            VertiportCollector vertiportCollector = new VertiportCollector(currentTrip, networkCar, vertiportsCandidatesOBUAM);
            vertiportCollector.neighbourVertiportCandidateTimeAndDistanceCalculatorForCompare();
            // provide information after each 10000 trips
            if (i % 10000 == 0) {
                System.out.println("Trip " + i + " is processed.");
            }
        }

        // Generate a List from 0-73
        List<Integer> vertiportIDList = new ArrayList<>();
        for (int i=0;i<vertiportsCandidatesOBUAM.size();i++){
            vertiportIDList.add(i);
        }
        double Score=calculateSelectionScore(vertiportsCandidatesOBUAM,vertiportIDList,uamEnabledTrips);
        System.out.println("The score of the current Vertiport Selection is: "+Score);

        // Save the access and egress time and distance for each vertiport candidate of each trip and store in a file}
    }
    private static LeastCostPathCalculator.Path estimatePath(Link from, Link to, double departureTime, Network carNetwork,
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
                LeastCostPathCalculator.Path path = estimatePath(from, to, trip.departureTime, networkCar, plcpccar);
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
    public static double calculateSelectionScore( List<Vertiport> vertiportCandidate,List<Integer> chosenVertiportID, List< TripItem > deserializedTripItems) throws IOException {
        // 实现适应度函数的具体逻辑
        double sumGeneralizedCost=0.0;
        double sumVertiportConstructionCost=0.0;
        double savedGeneralizedCost=0.0;
        double score=0.0;
        List<TripItem> uamAvailableTrips = new ArrayList<>();
        // calculate the sum of vertiport construction cost
        for (Integer vertiportID:chosenVertiportID) {
            sumVertiportConstructionCost=sumVertiportConstructionCost+vertiportCandidate.get(vertiportID).constructionCost;
        }

        for (TripItem tripItem : deserializedTripItems)
        {
            tripItem.isUAMAvailable = false;
            tripItem.uamUtility=-9999;
            List<Vertiport> originNeighbourVertiports = findAvailableNeighbourVertiports(chosenVertiportID, tripItem.originNeighborVertiportCandidates);
            List<Vertiport> destinationNeighbourVertiports = findAvailableNeighbourVertiports(chosenVertiportID, tripItem.destinationNeighborVertiportCandidates);
            if (originNeighbourVertiports.size() > 0 && destinationNeighbourVertiports.size() > 0) {
                tripItem.isUAMAvailable = true;
                tripItem.originNeighborVertiports = originNeighbourVertiports;
                tripItem.destinationNeighborVertiports = destinationNeighbourVertiports;
                uamAvailableTrips.add(tripItem);}
        }



        for (TripItem tripItem : uamAvailableTrips) {

            for (Vertiport vertiport : tripItem.originNeighborVertiports) {
                tripItem.originNeighborVertiportsTimeAndDistance.put(vertiport,tripItem.originNeighborVertiportCandidatesTimeAndDistance.get(vertiport));
            }
            for (Vertiport vertiport : tripItem.destinationNeighborVertiports) {
                tripItem.destinationNeighborVertiportsTimeAndDistance.put(vertiport,tripItem.destinationNeighborVertiportCandidatesTimeAndDistance.get(vertiport));
            }

            // Initialize the Vertiport Allocation
            // Find the vertiport pair for a trip with the lowest uam generalized cost, the access and ergress vertiport could not be the same
            double lowestUAMGeneralizedCost = Double.MAX_VALUE;
            Vertiport originVertiport = null;
            Vertiport destinationVertiport = null;
            for (Vertiport origin : tripItem.originNeighborVertiports) {
                for (Vertiport destination :  tripItem.destinationNeighborVertiports) {
                    if (origin.ID != destination.ID) {
                        double accessTime = tripItem.originNeighborVertiportsTimeAndDistance.get(origin).get("travelTime");
                        double egressTime = tripItem.destinationNeighborVertiportsTimeAndDistance.get(destination).get("travelTime");
                        double accessDistance = tripItem.originNeighborVertiportsTimeAndDistance.get(origin).get("distance");
                        double egressDistance = tripItem.destinationNeighborVertiportsTimeAndDistance.get(destination).get("distance");
                        double accessCost =0;
                        double egressCost =0;
                        double flightDistance= calculateEuciDistance(origin.coord,destination.coord);
                        double flightTime=flightDistance/flightSpeed+takeOffLandingTime;
                        double flightCost=6.1+ calculateEuciDistance(origin.coord,destination.coord)/1000*0.6;
                        double uamTravelTime=accessTime+egressTime+flightTime+UAM_PROCESS_TIME;
                        if (tripItem.accessMode.equals("car") ){
                            accessCost=accessDistance/1000*0.42;
                        }
                        if (tripItem.egressMode.equals("car") ){
                            egressCost=egressDistance/1000*0.42;
                        }
                        double UAMCost=accessCost+egressCost+flightCost;
                        double UAMGeneralizedCost=UAMCost+uamTravelTime*tripItem.VOT;
                        if (UAMGeneralizedCost < lowestUAMGeneralizedCost) {
                            tripItem.uamTravelTime=uamTravelTime;
                            tripItem.UAMCost=UAMCost;
                            tripItem.UAMUtilityVar=-2.48*UAMCost/100-4.28*flightTime/6000-6.79*(uamTravelTime-flightTime)/6000;
                            tripItem.uamUtility=tripItem.UAMUtilityFix+tripItem.UAMUtilityVar;
                            tripItem.UAMGeneralizedCost=UAMGeneralizedCost;
                            tripItem.accessVertiport = origin;
                            tripItem.egressVertiport = destination;
                            lowestUAMGeneralizedCost = UAMGeneralizedCost;
                        }
                    }
                }
            }
            // determine the probability of mode choice of each trip
            tripItem.uamProbability=calculateModeProbability(tripItem.uamUtility,tripItem.carUtility,tripItem.ptUtility).get(0);
            tripItem.carProbability=calculateModeProbability(tripItem.uamUtility,tripItem.carUtility,tripItem.ptUtility).get(1);
            tripItem.ptProbability=calculateModeProbability(tripItem.uamUtility,tripItem.carUtility,tripItem.ptUtility).get(2);
            double generalizedCostOneTripBefore=tripItem.carGeneralizedCost*calculateModeProbability(-9999,tripItem.carUtility,tripItem.ptUtility).get(1)+tripItem.ptGeneralizedCost*calculateModeProbability(-9999,tripItem.carUtility,tripItem.ptUtility).get(2);
            double generalizedCostOneTripAfter=tripItem.UAMGeneralizedCost*tripItem.uamProbability+tripItem.carGeneralizedCost*tripItem.carProbability+tripItem.ptGeneralizedCost*tripItem.ptProbability;
            double savedGeneralizedCostOneTrip=generalizedCostOneTripBefore-generalizedCostOneTripAfter;
            if (savedGeneralizedCostOneTrip<0){
                savedGeneralizedCostOneTrip=0;
            }
            if (tripItem.tripPurpose.startsWith("H")){
                savedGeneralizedCostOneTrip=savedGeneralizedCostOneTrip*2;
            }
            savedGeneralizedCost=savedGeneralizedCost+savedGeneralizedCostOneTrip;


        }



            /*
            int carCount=0;
            int ptCount=0;
            int uamCount=0;
            // generate 100 scenarios of each trip
            int j=0;
            for (int modeChoiceIterator=0;modeChoiceIterator<100;modeChoiceIterator++){
                j=modeChoiceIterator;
                ModeDecider modeDecider=new ModeDecider(tripItem.carProbability,tripItem.ptProbability,tripItem.uamProbability);
                String mode=modeDecider.decideMode();
                if (mode.equals("car")){
                    generalizedCostOneTrip=generalizedCostOneTrip+tripItem.carGeneralizedCost;
                    carCount++;
                }
                if (mode.equals("pt")){
                    generalizedCostOneTrip=generalizedCostOneTrip+tripItem.ptGeneralizedCost;
                    ptCount++;
                }
                if (mode.equals("uam")){
                    generalizedCostOneTrip=generalizedCostOneTrip+tripItem.UAMGeneralizedCost;
                    uamCount++;
                }
                if (tripItem.carProbability>0.9999 || tripItem.ptProbability>0.9999 || tripItem.uamProbability>0.9999){
                    carCount=carCount*100;
                    ptCount=ptCount*100;
                    uamCount=uamCount*100;
                    break;
                }
            }
            // calculate the average generalized cost of each trip
            generalizedCostOneTrip=generalizedCostOneTrip/(j+1);

             */




        score=savedGeneralizedCost;
        return score;
    }
}
