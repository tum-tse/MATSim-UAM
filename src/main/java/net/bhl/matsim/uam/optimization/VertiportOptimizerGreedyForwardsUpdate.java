package net.bhl.matsim.uam.optimization;

import net.bhl.matsim.uam.optimization.utils.ScenarioSpecific;
import net.bhl.matsim.uam.optimization.utils.TripItemForOptimization;
import org.apache.log4j.Logger;
import org.matsim.api.core.v01.Coord;
import org.matsim.utils.MemoryObserver;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.Callable;
import java.util.concurrent.atomic.AtomicInteger;

public class VertiportOptimizerGreedyForwardsUpdate {
    public static Logger log = Logger.getLogger(VertiportOptimizerGreedyForwardsUpdate.class);
    private static String vertiportCandidateFile;
    private static String fileName;
    private static  int sampleSize;
    private static long[] RANDOM_SEEDS;
    private static int num_of_run;
    private static final int MEMORY_CHECK_INTERVAL = 600;
    public static double flightSpeed; // m/s
    public static double UAM_PROCESS_TIME; // s
    public static double takeOffLandingTime; // s
    public static int NUM_OF_SELECTED_VERTIPORTS;
    private static double UAM_FIX_COST;
    private static double UAM_KM_COST;
    private static double CAR_COST;
    private static double UAM_UTILITY_COST_PARAMETER;
    private static double UAM_UTILITY_FLIGHT_TIME_PARAMETER;
    private static double UAM_UTILITY_WAIT_TIME_PARAMETER;
    private static boolean considerReturnTrip;

    public static double calculateEuciDistance(Coord coord1, Coord coord2) {
        double euciDistance = Math.sqrt(Math.pow(coord1.getX() - coord2.getX(), 2) + Math.pow(coord1.getY() - coord2.getY(), 2));
        return euciDistance;
    }
    public static void main(String[] args) throws IOException {
        MemoryObserver.start(MEMORY_CHECK_INTERVAL);
        // Provide the file via program arguments
        if (args.length > 0) {
            fileName = args[0];
            vertiportCandidateFile = args[1];
            sampleSize = Integer.parseInt(args[2]);
            num_of_run = Integer.parseInt(args[3]);
            // Provide the random seeds via the next arguments
            RANDOM_SEEDS = new long[num_of_run];
            for (int i = 0; i < num_of_run; i++) {
                RANDOM_SEEDS[i] = Long.parseLong(args[4 + i]);
            }

        }
        ScenarioSpecific scenarioSpecific = new ScenarioSpecific("Munich_A");
        scenarioSpecific.buildScenario();
        NUM_OF_SELECTED_VERTIPORTS = scenarioSpecific.num_of_selected_vertiports;
        UAM_FIX_COST = scenarioSpecific.uam_fix_cost;
        UAM_KM_COST = scenarioSpecific.uam_km_cost;
        CAR_COST = scenarioSpecific.car_km_cost;
        flightSpeed = scenarioSpecific.flight_speed;
        UAM_PROCESS_TIME = scenarioSpecific.uam_process_time;
        takeOffLandingTime = scenarioSpecific.uam_take_off_landing_time;
        UAM_UTILITY_COST_PARAMETER = scenarioSpecific.uam_utility_cost_parameter;
        UAM_UTILITY_FLIGHT_TIME_PARAMETER = scenarioSpecific.uam_utility_flight_time_parameter;
        UAM_UTILITY_WAIT_TIME_PARAMETER = scenarioSpecific.uam_utility_waiting_time_parameter;
        considerReturnTrip = scenarioSpecific.consider_return_trip;


        log.info("Loading the vertiport candidates...");
        VertiportReader vertiportReader = new VertiportReader();
        List<Vertiport> vertiportsCandidates = vertiportReader.getVertiports(vertiportCandidateFile);

        log.info("Finished loading the vertiport candidates.");
        log.info("Loading the trips...");
        List<TripItemForOptimization> deserializedTripItemForOptimizations = deserializeTripItems(fileName);
        log.info("Finished loading the trips.");

        for (int run_index = 0; run_index < num_of_run; run_index++) {
            log.info("The " + run_index + "th run starts.");
            double maxScore = Double.NEGATIVE_INFINITY;
            Random random = new Random(RANDOM_SEEDS[run_index]);
            HashMap<List<Integer>, Double> vertiportPairsScore = new HashMap<>();
            // Calculate the score of each two vertiports selection
            log.info("Calculating the score of each two vertiports selection...");
            for (TripItemForOptimization tripItemForOptimization : deserializedTripItemForOptimizations) {
                tripItemForOptimization.tempSavedGeneralizedCostsMap.clear();
                tripItemForOptimization.savedGeneralizedCost = 0;
                tripItemForOptimization.UAMUtilityVar=Double.NEGATIVE_INFINITY;
                tripItemForOptimization.isUAMAvailable=false;
            }
            HashMap<List<Integer>, Double> vertiportPairsAndScores = new HashMap<>();
            for (int i = 0; i < vertiportsCandidates.size(); i++) {
                for (int j = i + 1; j < vertiportsCandidates.size(); j++) {
                    List<Integer> vertiportPair = new ArrayList<>();
                    vertiportPair.add(i);
                    vertiportPair.add(j);
                    vertiportPairsAndScores.put(vertiportPair, 0.0);
                }
            }
            // Calculate the score of each two vertiports pair, use parallel computing
            // get the available processors
            int processors = Runtime.getRuntime().availableProcessors();
            int batchSize = vertiportPairsAndScores.size() / processors;
            ExecutorService executor = Executors.newFixedThreadPool(processors);
            List<Future<HashMap<List<Integer>, Double>>> futures = new ArrayList<>();
            AtomicInteger count = new AtomicInteger();
            for (int i = 0; i < processors; i++) {
                int start = i * batchSize;
                int end = (i == processors - 1) ? vertiportPairsAndScores.size() : (i + 1) * batchSize;
                HashMap<List<Integer>, Double> subMap = new HashMap<>();
                int index = 0;
                for (Map.Entry<List<Integer>, Double> entry : vertiportPairsAndScores.entrySet()) {
                    if (index >= start && index < end) {
                        subMap.put(entry.getKey(), entry.getValue());
                    }
                    index++;
                }
                Callable<HashMap<List<Integer>, Double>> task = () -> {
                    for (Map.Entry<List<Integer>, Double> entry : subMap.entrySet()) {
                        List<Integer> vertiportPair = entry.getKey();
                        double score = calculateSelectionScore(vertiportPair, deserializedTripItemForOptimizations, random);
                        entry.setValue(score);
                        count.getAndIncrement();
                        if (count.get() % 1000 == 0) {
                            log.info("Finish calculation of "+ count +" vertiport pairs.");
                        }

                    }
                    return subMap;
                };
                futures.add(executor.submit(task));
            }
            // Collect the results in the futures to the vertiportPairsScore
            for (Future<HashMap<List<Integer>, Double>> future : futures) {
                try {
                    HashMap<List<Integer>, Double> subMap = future.get();
                    vertiportPairsScore.putAll(subMap);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
            executor.shutdown();
            List<Integer> currentSelectedVertiportsID = new ArrayList<>();
            List<Integer> remainVetiportsCandidatesID = new ArrayList<>();
            for (int i = 0; i < vertiportsCandidates.size(); i++) {
                remainVetiportsCandidatesID.add(i);
            }

            while (currentSelectedVertiportsID.size() != NUM_OF_SELECTED_VERTIPORTS) {
                // select the first two vertiports
                if (currentSelectedVertiportsID.isEmpty()) {
                    // select the first vertiport pair
                    Integer maxA = null;
                    Integer maxB = null;
                    maxScore = Double.NEGATIVE_INFINITY;
                    for (Map.Entry<List<Integer>, Double> entry : vertiportPairsScore.entrySet()) {
                        if (entry.getValue() > maxScore) {
                            maxScore = entry.getValue();
                            maxA = entry.getKey().get(0);
                            maxB = entry.getKey().get(1);
                        }
                    }
                    currentSelectedVertiportsID.add(maxA);
                    currentSelectedVertiportsID.add(maxB);
                    remainVetiportsCandidatesID.remove(maxA);
                    remainVetiportsCandidatesID.remove(maxB);
                    calculateSelectionScore(currentSelectedVertiportsID,currentSelectedVertiportsID ,deserializedTripItemForOptimizations, random);
                    for (TripItemForOptimization tripItemForOptimization : deserializedTripItemForOptimizations) {
                        if (tripItemForOptimization.tempSavedGeneralizedCostsMap.containsKey(currentSelectedVertiportsID)) {
                            tripItemForOptimization.savedGeneralizedCost = tripItemForOptimization.tempSavedGeneralizedCostsMap.get(currentSelectedVertiportsID);
                        }
                        tripItemForOptimization.tempSavedGeneralizedCostsMap.clear();
                    }
                    log.info("The first two vertiport pair is: " + maxA + " and " + maxB);
                }
                else {
                    int processors2 = Runtime.getRuntime().availableProcessors();
                    int batchSize2 = remainVetiportsCandidatesID.size() / processors2;
                    ExecutorService executor2 = Executors.newFixedThreadPool(processors2);
                    HashMap<List<Integer>,Double> allPossibleVertiportPairsAndScores = new HashMap<>();
                    for (Integer integer : remainVetiportsCandidatesID) {
                        List<Integer> vertiportPair = new ArrayList<>();
                        vertiportPair.add(integer);
                        vertiportPair.addAll(currentSelectedVertiportsID);
                        allPossibleVertiportPairsAndScores.put(vertiportPair, 0.0);
                    }
                    List<Future<HashMap<List<Integer>, Double>>> futures2 = new ArrayList<>();
                    for (int i = 0; i < processors; i++) {
                        int start = i * batchSize2;
                        int end = (i == processors - 1) ? allPossibleVertiportPairsAndScores.size() : (i + 1) * batchSize2;
                        HashMap<List<Integer>, Double> subMap = new HashMap<>();
                        int index = 0;
                        for (Map.Entry<List<Integer>, Double> entry : allPossibleVertiportPairsAndScores.entrySet()) {
                            if (index >= start && index < end) {
                                subMap.put(entry.getKey(), entry.getValue());
                            }
                            index++;
                        }
                        Callable<HashMap<List<Integer>, Double>> task = () -> {
                            for (Map.Entry<List<Integer>, Double> entry : subMap.entrySet()) {
                                List<Integer> vertiportPair = entry.getKey();
                                List<Integer> selectionDifference = new ArrayList<>();
                                selectionDifference.add(vertiportPair.get(0))  ;
                                double score = calculateSelectionScore(vertiportPair, selectionDifference, deserializedTripItemForOptimizations, random);
                                entry.setValue(score);
                            }
                            return subMap;
                        };
                        futures2.add(executor2.submit(task));
                    }
                    // Collect the results in the futures to the allPossibleVertiportPairsAndScores
                    for (Future<HashMap<List<Integer>, Double>> future : futures2) {
                        try {
                            HashMap<List<Integer>, Double> subMap = future.get();
                            allPossibleVertiportPairsAndScores.putAll(subMap);
                        } catch (Exception e) {
                            e.printStackTrace();
                        }
                    }
                    executor2.shutdown();
                    // Find the vertiport pair with the highest score
                    Integer [] max=null;
                    maxScore = Double.NEGATIVE_INFINITY;
                    for (Map.Entry<List<Integer>, Double> entry : allPossibleVertiportPairsAndScores.entrySet()) {
                        if (entry.getValue() > maxScore) {
                            maxScore = entry.getValue();
                            max = new Integer[entry.getKey().size()];
                            entry.getKey().toArray(max);
                        }
                    }
                    int newVertiportID = 0;
                    // find the new selected vertiport, which is in the max array but not in the currentSelectedVertiportsID
                    for (int i = 0; i < Objects.requireNonNull(max).length; i++) {
                        if (!currentSelectedVertiportsID.contains(max[i])) {
                            newVertiportID = max[i];
                            currentSelectedVertiportsID.add(max[i]);
                            remainVetiportsCandidatesID.remove(max[i]);
                            break;
                        }
                    }
                    List<Integer> newSelectedVertiportsID = new ArrayList<>();
                    newSelectedVertiportsID.add(newVertiportID);
                    for (TripItemForOptimization tripItemForOptimization : deserializedTripItemForOptimizations) {
                        if (tripItemForOptimization.tempSavedGeneralizedCostsMap.containsKey(newSelectedVertiportsID)) {
                        tripItemForOptimization.savedGeneralizedCost = tripItemForOptimization.tempSavedGeneralizedCostsMap.get(newSelectedVertiportsID);}
                        tripItemForOptimization.tempSavedGeneralizedCostsMap.clear();
                    }
                    log.info("The " + currentSelectedVertiportsID.size() + "th vertiport is "+ newVertiportID);
                }

            }

            log.info("The selected vertiports are: " + currentSelectedVertiportsID);
            log.info("The score of the selected vertiports is: " + maxScore);
        }
    }





    public static List<TripItemForOptimization> deserializeTripItems(String fileName) {
        List<TripItemForOptimization> tripItemForOptimizations = new ArrayList<>();

        try  {
            FileInputStream fileIn = new FileInputStream(fileName);
            ObjectInputStream objectIn = new ObjectInputStream(fileIn)  ;
            tripItemForOptimizations = (List<TripItemForOptimization>) objectIn.readObject();


        } catch (Exception e) {
            e.printStackTrace();
        }

        return tripItemForOptimizations;
    }
    public static List<Vertiport> findAvailableNeighbourVertiports(List<Integer> intList, List<Vertiport> vertiportList) {
        List<Vertiport> duplicates = new ArrayList<>();

        for (Vertiport vertiport : vertiportList) {
            if (intList.contains(vertiport.ID)) {
                duplicates.add(vertiport);
            } else {
                continue;
            }
        }

        return duplicates;
    }

    public static double calculateSelectionScore(List<Integer> chosenVertiportID, List<Integer> selectionDifference, List<TripItemForOptimization> deserializedTripItemForOptimizations, Random random) {
        List<TripItemForOptimization> uamAvailableTrips = new ArrayList<>();
        List<TripItemForOptimization> differenceTrips = new ArrayList<>();
        double totalSavedGeneralizedCost = 0;
        for (TripItemForOptimization tripItemForOptimization : deserializedTripItemForOptimizations) {
            boolean added = false;
            for (Vertiport vertiport : tripItemForOptimization.originNeighborVertiportCandidates) {
                if (selectionDifference.contains(vertiport.ID)) {
                    differenceTrips.add(tripItemForOptimization);
                    added = true;
                    break;
                }
            }
            if (added) {
                continue;
            }
            for (Vertiport vertiport : tripItemForOptimization.destinationNeighborVertiportCandidates) {
                if (selectionDifference.contains(vertiport.ID)) {
                    differenceTrips.add(tripItemForOptimization);
                    break;
                }
            }
        }

        for (TripItemForOptimization tripItemForOptimization : differenceTrips) {
            List<Vertiport> originNeighbourVertiports = findAvailableNeighbourVertiports(chosenVertiportID, tripItemForOptimization.originNeighborVertiportCandidates);
            List<Vertiport> destinationNeighbourVertiports = findAvailableNeighbourVertiports(chosenVertiportID, tripItemForOptimization.destinationNeighborVertiportCandidates);
            if (!originNeighbourVertiports.isEmpty() && !destinationNeighbourVertiports.isEmpty()) {
                tripItemForOptimization.isUAMAvailable = true;
                tripItemForOptimization.originNeighborVertiports = originNeighbourVertiports;
                tripItemForOptimization.destinationNeighborVertiports = destinationNeighbourVertiports;
                uamAvailableTrips.add(tripItemForOptimization);
            } else {
                tripItemForOptimization.isUAMAvailable = false;
                tripItemForOptimization.UAMUtilityVar = Double.NEGATIVE_INFINITY;
                tripItemForOptimization.tempSavedGeneralizedCostsMap.put(selectionDifference, 0.0);
            }
        }


        for (TripItemForOptimization tripItemForOptimization : uamAvailableTrips) {
            calculateTripSavedCost(tripItemForOptimization, selectionDifference,random);
        }

        for (TripItemForOptimization tripItemForOptimization : deserializedTripItemForOptimizations) {
            if (differenceTrips.contains(tripItemForOptimization)) {
                totalSavedGeneralizedCost += tripItemForOptimization.tempSavedGeneralizedCostsMap.get(selectionDifference);
            }
            else {
                totalSavedGeneralizedCost += tripItemForOptimization.savedGeneralizedCost;
            }
        }
        return totalSavedGeneralizedCost;
    }

    public static synchronized void calculateTripSavedCost(TripItemForOptimization tripItemForOptimization, List<Integer> selectionDifference,Random random) {
        for (Vertiport vertiport : tripItemForOptimization.originNeighborVertiports) {
            tripItemForOptimization.originNeighborVertiportsTimeAndDistance.put(vertiport, tripItemForOptimization.originNeighborVertiportCandidatesTimeAndDistance.get(vertiport));
        }
        for (Vertiport vertiport : tripItemForOptimization.destinationNeighborVertiports) {
            tripItemForOptimization.destinationNeighborVertiportsTimeAndDistance.put(vertiport, tripItemForOptimization.destinationNeighborVertiportCandidatesTimeAndDistance.get(vertiport));
        }
        double objectiveFunctionBefore;
        double objectiveFunctionAfter;
        // Initialize the Vertiport Allocation
        // Find the vertiport pair for a trip with the lowest uam generalized cost, the access and ergress vertiport could not be the same
        double lowestUAMGeneralizedCost = Double.MAX_VALUE;
        for (Vertiport origin : tripItemForOptimization.originNeighborVertiports) {
            for (Vertiport destination :  tripItemForOptimization.destinationNeighborVertiports) {
                if (origin.ID != destination.ID) {
                    double accessTime = tripItemForOptimization.originNeighborVertiportsTimeAndDistance.get(origin).get("travelTime");
                    double egressTime = tripItemForOptimization.destinationNeighborVertiportsTimeAndDistance.get(destination).get("travelTime");
                    double accessDistance = tripItemForOptimization.originNeighborVertiportsTimeAndDistance.get(origin).get("distance");
                    double egressDistance = tripItemForOptimization.destinationNeighborVertiportsTimeAndDistance.get(destination).get("distance");
                    double accessCost =0;
                    double egressCost =0;
                    double flightDistance= calculateEuciDistance(origin.coord,destination.coord);
                    double flightTime=flightDistance/flightSpeed+takeOffLandingTime;
                    double flightCost=UAM_FIX_COST+ calculateEuciDistance(origin.coord,destination.coord)/1000*UAM_KM_COST;
                    double uamTravelTime=accessTime+egressTime+flightTime+UAM_PROCESS_TIME;
                    if (tripItemForOptimization.accessMode.equals("car") ){
                        accessCost=accessDistance/1000*CAR_COST;
                    }
                    if (tripItemForOptimization.egressMode.equals("car") ){
                        egressCost=egressDistance/1000*CAR_COST;
                    }
                    double UAMCost=accessCost+egressCost+flightCost;
                    double UAMGeneralizedCost=UAMCost+uamTravelTime* tripItemForOptimization.VOT;
                    if (UAMGeneralizedCost < lowestUAMGeneralizedCost) {
                        tripItemForOptimization.uamTravelTime=uamTravelTime;
                        tripItemForOptimization.UAMCost=UAMCost;
                        tripItemForOptimization.UAMUtilityVar=UAM_UTILITY_COST_PARAMETER*UAMCost/100+UAM_UTILITY_FLIGHT_TIME_PARAMETER*flightTime/6000+UAM_UTILITY_WAIT_TIME_PARAMETER*(uamTravelTime-flightTime)/6000;
                        tripItemForOptimization.uamUtility= tripItemForOptimization.UAMUtilityFix+ tripItemForOptimization.UAMUtilityVar;
                        tripItemForOptimization.UAMGeneralizedCost=UAMGeneralizedCost;
                        tripItemForOptimization.accessVertiport = origin;
                        tripItemForOptimization.egressVertiport = destination;
                        lowestUAMGeneralizedCost = UAMGeneralizedCost;
                    }
                }
            }
        }
        // determine the probability of mode choice of each trip

        // Create a double list to store the probability of each mode
        ModeDecider modeDecider=new ModeDecider(tripItemForOptimization.uamUtility,tripItemForOptimization.carUtility,tripItemForOptimization.ptUtility,random);
        Double [] modeSamples=modeDecider.sample(sampleSize);

        objectiveFunctionBefore= tripItemForOptimization.currentGeneralizedCost;
        objectiveFunctionAfter= modeSamples[1]* tripItemForOptimization.carGeneralizedCost+modeSamples[2]* tripItemForOptimization.ptGeneralizedCost+modeSamples[0]* tripItemForOptimization.UAMGeneralizedCost;

        double savedGeneralizedCostOneTrip=objectiveFunctionBefore-objectiveFunctionAfter;
        if (savedGeneralizedCostOneTrip<0){
            savedGeneralizedCostOneTrip=0;
        }
        if (tripItemForOptimization.tripPurpose.startsWith("H") && considerReturnTrip){
            savedGeneralizedCostOneTrip=savedGeneralizedCostOneTrip*2;
        }

        // Include your logic here that was previously inside the loop over all trips
        tripItemForOptimization.tempSavedGeneralizedCostsMap.put(selectionDifference, savedGeneralizedCostOneTrip);
    }
    public static double calculateSelectionScore( List<Integer> chosenVertiportID, List<TripItemForOptimization> deserializedTripItemForOptimizations,Random random) {
        List<TripItemForOptimization> uamAvailableTrips = new ArrayList<>();
        double score = 0;

        for (TripItemForOptimization tripItemForOptimization : deserializedTripItemForOptimizations) {
            tripItemForOptimization.isUAMAvailable = false;
            tripItemForOptimization.uamUtility = Double.NEGATIVE_INFINITY;
            List<Vertiport> originNeighbourVertiports = findAvailableNeighbourVertiports(chosenVertiportID, tripItemForOptimization.originNeighborVertiportCandidates);
            List<Vertiport> destinationNeighbourVertiports = findAvailableNeighbourVertiports(chosenVertiportID, tripItemForOptimization.destinationNeighborVertiportCandidates);
            if (!originNeighbourVertiports.isEmpty() && !destinationNeighbourVertiports.isEmpty()) {
                tripItemForOptimization.isUAMAvailable = true;
                tripItemForOptimization.originNeighborVertiports = originNeighbourVertiports;
                tripItemForOptimization.destinationNeighborVertiports = destinationNeighbourVertiports;
                uamAvailableTrips.add(tripItemForOptimization);
            }
        }


        for (TripItemForOptimization tripItemForOptimization : uamAvailableTrips) {
            for (Vertiport vertiport : tripItemForOptimization.originNeighborVertiports) {
                tripItemForOptimization.originNeighborVertiportsTimeAndDistance.put(vertiport, tripItemForOptimization.originNeighborVertiportCandidatesTimeAndDistance.get(vertiport));
            }
            for (Vertiport vertiport : tripItemForOptimization.destinationNeighborVertiports) {
                tripItemForOptimization.destinationNeighborVertiportsTimeAndDistance.put(vertiport, tripItemForOptimization.destinationNeighborVertiportCandidatesTimeAndDistance.get(vertiport));
            }

            // Initialize the Vertiport Allocation
            // Find the vertiport pair for a trip with the lowest uam generalized cost, the access and ergress vertiport could not be the same
            double lowestUAMGeneralizedCost = Double.MAX_VALUE;
            for (Vertiport origin : tripItemForOptimization.originNeighborVertiports) {
                for (Vertiport destination : tripItemForOptimization.destinationNeighborVertiports) {
                    if (origin.ID != destination.ID) {
                        double accessTime = tripItemForOptimization.originNeighborVertiportsTimeAndDistance.get(origin).get("travelTime");
                        double egressTime = tripItemForOptimization.destinationNeighborVertiportsTimeAndDistance.get(destination).get("travelTime");
                        double accessDistance = tripItemForOptimization.originNeighborVertiportsTimeAndDistance.get(origin).get("distance");
                        double egressDistance = tripItemForOptimization.destinationNeighborVertiportsTimeAndDistance.get(destination).get("distance");
                        double accessCost = 0;
                        double egressCost = 0;
                        double flightDistance = calculateEuciDistance(origin.coord, destination.coord);
                        double flightTime = flightDistance / flightSpeed + takeOffLandingTime;
                        double flightCost = UAM_FIX_COST + calculateEuciDistance(origin.coord, destination.coord) / 1000 * UAM_KM_COST;
                        double uamTravelTime = accessTime + egressTime + flightTime + UAM_PROCESS_TIME;
                        if (tripItemForOptimization.accessMode.equals("car")) {
                            accessCost = accessDistance / 1000 * CAR_COST;
                        }
                        if (tripItemForOptimization.egressMode.equals("car")) {
                            egressCost = egressDistance / 1000 * CAR_COST;
                        }
                        double UAMCost = accessCost + egressCost + flightCost;
                        double UAMGeneralizedCost = UAMCost + uamTravelTime * tripItemForOptimization.VOT;
                        if (UAMGeneralizedCost < lowestUAMGeneralizedCost) {
                            tripItemForOptimization.uamTravelTime = uamTravelTime;
                            tripItemForOptimization.UAMCost = UAMCost;
                            tripItemForOptimization.UAMUtilityVar = UAM_UTILITY_COST_PARAMETER*UAMCost/100+UAM_UTILITY_FLIGHT_TIME_PARAMETER*flightTime/6000+UAM_UTILITY_WAIT_TIME_PARAMETER*(uamTravelTime-flightTime)/6000;
                            tripItemForOptimization.uamUtility = tripItemForOptimization.UAMUtilityFix + tripItemForOptimization.UAMUtilityVar;
                            tripItemForOptimization.UAMGeneralizedCost = UAMGeneralizedCost;
                            tripItemForOptimization.accessVertiport = origin;
                            tripItemForOptimization.egressVertiport = destination;
                            lowestUAMGeneralizedCost = UAMGeneralizedCost;
                        }
                    }
                }
            }
            // determine the probability of mode choice of each trip

            // Create a double list to store the probability of each mode
            ModeDecider modeDecider = new ModeDecider(tripItemForOptimization.uamUtility, tripItemForOptimization.carUtility, tripItemForOptimization.ptUtility, random);
            Double[] modeSamples = modeDecider.sample(sampleSize);
            double generalizedCostOneTripBefore = tripItemForOptimization.currentGeneralizedCost;
            double generalizedCostOneTripAfter = tripItemForOptimization.carGeneralizedCost * modeSamples[1] + tripItemForOptimization.ptGeneralizedCost * modeSamples[2] + tripItemForOptimization.UAMGeneralizedCost * modeSamples[0];
            double savedGeneralizedCostOneTrip = generalizedCostOneTripBefore - generalizedCostOneTripAfter;
            if (savedGeneralizedCostOneTrip < 0) {
                savedGeneralizedCostOneTrip = 0;
            }
            if (tripItemForOptimization.tripPurpose.startsWith("H")) {
                savedGeneralizedCostOneTrip = savedGeneralizedCostOneTrip * 2;
            }
            score += savedGeneralizedCostOneTrip;
        }
        return score;
    }



}