package net.bhl.matsim.uam.optimization;

import net.bhl.matsim.uam.optimization.utils.ScenarioSpecific;
import net.bhl.matsim.uam.optimization.utils.TripItemForOptimization;
import net.bhl.matsim.uam.optimization.utils.TripItemReaderForOptimization;
import org.apache.log4j.Logger;
import org.matsim.api.core.v01.Coord;
import org.matsim.core.config.Config;
import org.matsim.core.config.ConfigUtils;
import org.matsim.utils.MemoryObserver;

import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

public class VertiportOptimizerGreedyForwardsNew {
    public static Logger log = Logger.getLogger(VertiportOptimizerGreedyForwardsNew.class);
    private static String vertiportCandidateFile;
    private static String tripItemFile;
    private static  int sampleSize;
    private static long RANDOM_SEED;
    private static String configPath;
    private static final int MEMORY_CHECK_INTERVAL = 600;
    public static double flightSpeed; // m/s
    public static double UAM_PROCESS_TIME; // s
    public static double takeOffLandingTime; // s
    public static int NUM_OF_SELECTED_VERTIPORTS;
    private static double UAM_FIX_COST;
    private static double UAM_KM_COST;
    private static double CAR_COST;
    private static boolean considerReturnTrip;
    public static String scenarioName;
    public static double calculateEuciDistance(Coord coord1, Coord coord2) {
        double euciDistance = Math.sqrt(Math.pow(coord1.getX() - coord2.getX(), 2) + Math.pow(coord1.getY() - coord2.getY(), 2));
        return euciDistance;
    }
    public static void main(String[] args) throws Exception {
        MemoryObserver.start(MEMORY_CHECK_INTERVAL);
        // Provide the file via program arguments
        if (args.length > 0) {
            tripItemFile = args[0];
            configPath=args[1];
            vertiportCandidateFile = args[2];
            sampleSize = Integer.parseInt(args[3]);
            scenarioName = args[4];
            RANDOM_SEED = Long.parseLong(args[5]);
        }
        // Build the scenario of Munich
        ScenarioSpecific scenarioSpecific = new ScenarioSpecific(scenarioName);
        scenarioSpecific.buildScenario();

        // load the config
        Config config = ConfigUtils.loadConfig(configPath);

        // load the scenario specific parameters
        NUM_OF_SELECTED_VERTIPORTS = scenarioSpecific.num_of_selected_vertiports;
        UAM_FIX_COST = scenarioSpecific.uam_fix_cost;
        UAM_KM_COST = scenarioSpecific.uam_km_cost;
        CAR_COST = scenarioSpecific.car_km_cost;
        flightSpeed = scenarioSpecific.flight_speed;
        UAM_PROCESS_TIME = scenarioSpecific.uam_process_time;
        takeOffLandingTime = scenarioSpecific.uam_take_off_landing_time;
        considerReturnTrip = scenarioSpecific.consider_return_trip;

        log.info("Loading the vertiport candidates...");
        VertiportReader vertiportReader = new VertiportReader();
        List<Vertiport> vertiportsCandidates = vertiportReader.getVertiports(vertiportCandidateFile);
        log.info("Finished loading the vertiport candidates.");


        // Load the trip items
        log.info("Loading the trips...");
        TripItemReaderForOptimization tripItemReaderForOptimization = new TripItemReaderForOptimization(scenarioSpecific);
        List<TripItemForOptimization> tripItems = tripItemReaderForOptimization.getTripItemsForOptimization(tripItemFile);
        log.info("Finished loading the trips.");

        // Pre-calculate the access and egress time and distance for each vertiport candidate of each trip
        log.info("Start pre-calculating the access and egress time and distance for each vertiport candidate of each trip...");
        AccessEgressCostCalculator accessEgressCostCalculator = new AccessEgressCostCalculator(tripItems, vertiportsCandidates, config, scenarioSpecific);
        accessEgressCostCalculator.calculateAccessEgressCost();
        log.info("Finished pre-calculating the access and egress time and distance for each vertiport candidate of each trip.");

        // Initialize the UAM enabled trips and UAM utility
        log.info("Start the GRD Algorithm...");



            // start record the time
            long startTime = System.currentTimeMillis();
            double maxScore = Double.NEGATIVE_INFINITY;
            Random random = new Random(RANDOM_SEED);
            ConcurrentHashMap<List<Integer>, Double> vertiportPairsScore = new ConcurrentHashMap<>();
            // Calculate the score of each two vertiports selection
            log.info("Calculating the score of each two vertiports selection...");
            for (TripItemForOptimization tripItemForOptimization : tripItems) {
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
            AtomicInteger count = new AtomicInteger();
            ExecutorService executor = Executors.newFixedThreadPool(processors);
        try {
            List<Future<HashMap<List<Integer>, Double>>> futures = new ArrayList<>();
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
                        double score = calculateSelectionScore(vertiportPair, vertiportsCandidates, tripItems, random, scenarioSpecific);
                        entry.setValue(score);
                        count.getAndIncrement();
                        if (vertiportPairsAndScores.size()<100 || (count.get()) % 1000 == 0) {
                            log.info("Calculation completion: " + count.get() + "/" + vertiportPairsAndScores.size() + " ("
                                    + String.format("%.0f", (double) (count.get()) / vertiportPairsAndScores.size() * 100) + "%).");
                        }

                    }
                    return subMap;
                };
                futures.add(executor.submit(task));
            }
            // 等待所有任务完成
            for (Future<HashMap<List<Integer>, Double>> future : futures) {
                try {
                    HashMap<List<Integer>, Double> subMap = future.get();
                    vertiportPairsScore.putAll(subMap);
                } catch (ExecutionException e) {
                    Throwable cause = e.getCause();  // 获取实际异常
                    cause.printStackTrace();
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                }
            }
        } finally {
            executor.shutdown();
            try {
                if (!executor.awaitTermination(800, TimeUnit.MILLISECONDS)) {
                    executor.shutdownNow();
                }
            } catch (InterruptedException e) {
                executor.shutdownNow();
                Thread.currentThread().interrupt();  // 保留中断状态
            }
        }


        List<Integer> currentSelectedVertiportsID = new CopyOnWriteArrayList<>();
            List<Integer> remainVetiportsCandidatesID = new CopyOnWriteArrayList<>();

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
                    calculateSelectionScore(currentSelectedVertiportsID, vertiportsCandidates ,tripItems, random,scenarioSpecific);
                    log.info("The first two vertiport pair is: " + maxA + " and " + maxB);
                }
                else {
                    // select the next vertiport
                    Integer maxVertiportID = null;
                    maxScore = Double.NEGATIVE_INFINITY;

                    for (Integer vertiportID : remainVetiportsCandidatesID) {
                        double score = 0;
                        for (Integer selectedVertiportID : currentSelectedVertiportsID) {
                            List<Integer> vertiportsPair = new ArrayList<>();
                            vertiportsPair.add(selectedVertiportID);
                            vertiportsPair.add(vertiportID);
                            score = score + getPairScore(vertiportsPair, vertiportPairsScore);
                        }
                        if (score > maxScore) {
                            maxScore = score;
                            maxVertiportID = vertiportID;
                        }
                    }
                    currentSelectedVertiportsID.add(maxVertiportID);
                    remainVetiportsCandidatesID.remove(maxVertiportID);
                    log.info("The " + currentSelectedVertiportsID.size() + "th vertiport is: " + maxVertiportID);

                }

            }

            log.info("The selected vertiports are: " + currentSelectedVertiportsID);
            log.info("The score of the selected vertiports is: " + calculateSelectionScore(currentSelectedVertiportsID, vertiportsCandidates, tripItems, random, scenarioSpecific));
            long endTime = System.currentTimeMillis();
            log.info("The " + "run finishes. The time used is: " + (endTime - startTime) / 1000 + "s.");

    }




    public static double getPairScore(List<Integer> vertiportsPair, ConcurrentHashMap<List<Integer>, Double> vertiportPairsScore) {
        // sort the vertiport pair
        Collections.sort(vertiportsPair);
        return vertiportPairsScore.get(vertiportsPair);

    }
    public static List<Vertiport> findAvailableNeighbourVertiports(List<Integer> intList, List<Vertiport> vertiportList) {
        List<Vertiport> duplicates = new CopyOnWriteArrayList<>();

        for (Vertiport vertiport : vertiportList) {
            if (intList.contains(vertiport.ID)) {
                duplicates.add(vertiport);
            } else {
                continue;
            }
        }

        return duplicates;
    }

    public static double calculateSelectionScore(List<Integer> chosenVertiportID, List<Vertiport> vertiportsCandidates, List<TripItemForOptimization> deserializedTripItemForOptimizations, Random random, ScenarioSpecific scenarioSpecific) throws Exception {
        double totalSavedGeneralizedCost = 0;
        for (TripItemForOptimization tripItemForOptimization : deserializedTripItemForOptimizations) {
                List<Vertiport> originNeighbourVertiports = findAvailableNeighbourVertiports(chosenVertiportID, tripItemForOptimization.originNeighborVertiportCandidates);
                List<Vertiport> destinationNeighbourVertiports = findAvailableNeighbourVertiports(chosenVertiportID, tripItemForOptimization.destinationNeighborVertiportCandidates);

                if (!originNeighbourVertiports.isEmpty() && !destinationNeighbourVertiports.isEmpty()) {
                    if (originNeighbourVertiports.size() > 1 || destinationNeighbourVertiports.size() > 1 || originNeighbourVertiports.get(0).ID != destinationNeighbourVertiports.get(0).ID) {
                        tripItemForOptimization.isUAMAvailable = true;
                        tripItemForOptimization.originNeighborVertiports = originNeighbourVertiports;
                        tripItemForOptimization.destinationNeighborVertiports = destinationNeighbourVertiports;
                        totalSavedGeneralizedCost += calculateTripSavedCost(tripItemForOptimization, vertiportsCandidates, random, scenarioSpecific);
                    }
                    else {
                        tripItemForOptimization.isUAMAvailable = false;
                        tripItemForOptimization.UAMUtilityVar = Double.NEGATIVE_INFINITY;
                    }
                }
                else {

                tripItemForOptimization.isUAMAvailable = false;
                tripItemForOptimization.UAMUtilityVar = Double.NEGATIVE_INFINITY;
            };
        }

        return totalSavedGeneralizedCost;
    }

    public static double calculateTripSavedCost(TripItemForOptimization tripItemForOptimization, List<Vertiport> vertiportsCandidates, Random random,ScenarioSpecific scenarioSpecific) {
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
                        tripItemForOptimization.UAMGeneralizedCost=UAMGeneralizedCost;
                        tripItemForOptimization.UAMUtilityVar= scenarioSpecific.calculateUAMUtilityVAR(flightTime,uamTravelTime-flightTime,UAMCost);
                        tripItemForOptimization.uamUtility= tripItemForOptimization.UAMUtilityFix+ tripItemForOptimization.UAMUtilityVar;

                        tripItemForOptimization.accessVertiport = origin;
                        tripItemForOptimization.egressVertiport = destination;
                        lowestUAMGeneralizedCost = UAMGeneralizedCost;
                    }
                }
            }
        }

        // Create a double list to store the probability of each mode
        ModeDecider modeDecider=new ModeDecider(tripItemForOptimization.uamUtility*3,tripItemForOptimization.carUtility*3,tripItemForOptimization.ptUtility*3,random);
        Double [] modeSamples=modeDecider.sample(sampleSize);
        tripItemForOptimization.uamProbability=modeSamples[0];
        tripItemForOptimization.carProbability=modeSamples[1];
        tripItemForOptimization.ptProbability=modeSamples[2];

        objectiveFunctionBefore= tripItemForOptimization.currentGeneralizedCost;
        objectiveFunctionAfter= modeSamples[1]*tripItemForOptimization.carGeneralizedCost+modeSamples[2]*tripItemForOptimization.ptGeneralizedCost+modeSamples[0]* tripItemForOptimization.UAMGeneralizedCost;

        double savedGeneralizedCostOneTrip=objectiveFunctionBefore-objectiveFunctionAfter;
        if (savedGeneralizedCostOneTrip<0){
            savedGeneralizedCostOneTrip=0;
        }
        if (tripItemForOptimization.tripPurpose.startsWith("H") && considerReturnTrip){
            savedGeneralizedCostOneTrip=savedGeneralizedCostOneTrip*2;
        }

        // Include your logic here that was previously inside the loop over all trips
       return savedGeneralizedCostOneTrip;
    }


}