package net.bhl.matsim.uam.optimization;

import net.bhl.matsim.uam.optimization.utils.*;
import org.apache.log4j.Logger;
import org.matsim.api.core.v01.Coord;
import org.matsim.core.config.Config;
import org.matsim.core.config.ConfigUtils;
import org.matsim.utils.MemoryObserver;
import smile.clustering.DBSCAN;
import smile.clustering.KMeans;
import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.util.*;
import java.util.concurrent.*;
import java.util.stream.Collectors;

public class SimulatedAnnealingForPartD {
    // This class provides the simulated annealing algorithm to solve the vertiport siting problem.

    public static final Logger log = Logger.getLogger(SimulatedAnnealingForPartD.class);

    public static String tripItemFile;
    public static String vertiportCandidateFile;
    public static double flightSpeed; // m/s
    public static double UAM_PROCESS_TIME; // s
    public static double takeOffLandingTime; // s
    public static int NUM_OF_SELECTED_VERTIPORTS;
    private static final int MEMORY_CHECK_INTERVAL=600;
    private static double UAM_FIX_COST;
    private static double UAM_KM_COST;
    private static double CAR_COST;
    private static boolean considerReturnTrip;
    private static  int sampleSize;
    private static long RANDOM_SEED;
    public static double CAR_EMISSION_FACTOR; // kg/km
    public static double PT_EMISSION_FACTOR; // kg/km
    public static double UAM_EMISSION_FACTOR; // kg/km
    public static double CARBON_EQUIVALENCE_FACTOR; // Euro/kgCO2
    private static int SIMULATION_HOURS;
    public static String configPath;
    public static double NEIGHBOR_DISTANCE;

    public static String scenarioConfigurations;
    public static void main(String[] args) throws Exception {
        MemoryObserver.start(MEMORY_CHECK_INTERVAL);
        // Provide the file via program arguments
        if (args.length > 0) {
            tripItemFile = args[0];
            configPath=args[1];
            vertiportCandidateFile = args[2];
            scenarioConfigurations = args[3];
        }
        // Build the scenario of Munich
        log.info("Building the scenario...");
        ScenarioSpecific scenarioSpecific = new ScenarioSpecific(scenarioConfigurations);
        scenarioSpecific.buildScenario();
        ClusteringConfiguration clusteringConfiguration = new ClusteringConfiguration(scenarioConfigurations);
        clusteringConfiguration.buildScenario();
        OptimizationConfiguration optimizationConfiguration = new OptimizationConfiguration(scenarioConfigurations);
        optimizationConfiguration.buildScenario();
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
        CAR_EMISSION_FACTOR = scenarioSpecific.car_emission_factor;
        PT_EMISSION_FACTOR = scenarioSpecific.pt_emission_factor;
        SIMULATION_HOURS = scenarioSpecific.simulation_hours;
        CARBON_EQUIVALENCE_FACTOR = scenarioSpecific.carbon_equivalent_cost;
        RANDOM_SEED =  scenarioSpecific.random_seed;
        sampleSize = scenarioSpecific.sampleSize;
        NEIGHBOR_DISTANCE = scenarioSpecific.neighbouring_distance;

        log.info("Loading the vertiport candidates...");
        VertiportReader vertiportReader = new VertiportReader();
        List<Vertiport> vertiportsCandidates = vertiportReader.getVertiportsWithNeighbors(vertiportCandidateFile);
        log.info("Finished loading the vertiport candidates.");

        // Cluster the vertiport candidates with DBSCAN clustering algorithm
        log.info("Start clustering the vertiport candidates...");
        //int[] clusterAssigment = clusterDataDBSCAN(vertiportsCandidates, EPSILON, MIN_PTS);
        int[] clusterAssigment = clusterDataKMeans(vertiportsCandidates, clusteringConfiguration.num_of_clusters, clusteringConfiguration.num_of_cluster_iterations, clusteringConfiguration.tolerance);
        List<Vertiport> clusteredVertiports = summarizeClusters(vertiportsCandidates, clusterAssigment);
        Map<Integer, List<Vertiport>> clusterResults = saveClusterResults(vertiportsCandidates, clusterAssigment);
        findVertiportsNeighbours(clusteredVertiports, NEIGHBOR_DISTANCE);
        log.info("Finished clustering the vertiport candidates.");

        log.info("Loading the trips...");
        TripItemReaderForOptimization tripItemReaderForOptimization = new TripItemReaderForOptimization(scenarioSpecific);
        List<TripItemForOptimization> tripItems = tripItemReaderForOptimization.getTripItemsForOptimization(tripItemFile);
        log.info("Finished loading the trips.");

        // Pre-calculate the access and egress time and distance for each vertiport candidate of each trip
        log.info("Start pre-calculating the access and egress time and distance for each vertiport candidate of each trip...");
        AccessEgressCostCalculator accessEgressCostCalculator = new AccessEgressCostCalculator(tripItems, clusteredVertiports, config, scenarioSpecific);
        accessEgressCostCalculator.calculateAccessEgressCost();
        log.info("Finished pre-calculating the access and egress time and distance for each vertiport candidate of each trip.");

        // Initialize the UAM enabled trips and UAM utility
        log.info("Start the simulated annealing algorithm...");

            for (TripItemForOptimization tripItemForOptimization : tripItems) {
                tripItemForOptimization.isUAMAvailable = false;
                tripItemForOptimization.uamUtility = Double.NEGATIVE_INFINITY;
                tripItemForOptimization.savedGeneralizedCost = 0;
            }

            // record the start time
            long startTime = System.currentTimeMillis();

            Random random = new Random(RANDOM_SEED);
            List<Integer> currentSolutionID = generateRandomSolution(random, clusteredVertiports, NUM_OF_SELECTED_VERTIPORTS);
            double currentEnergey = calculateFitness(currentSolutionID,  clusteredVertiports, tripItems,random,scenarioSpecific);
            for (TripItemForOptimization tripItemForOptimization : tripItems) {
                if (!tripItemForOptimization.tempSavedGeneralizedCosts.isEmpty()) {
                    tripItemForOptimization.savedGeneralizedCost = tripItemForOptimization.tempSavedGeneralizedCosts.get(0);
                }
            }
            // copy the tempConcurrentSaturationRates to concurrentSaturationRates and tempMaxSaturationRate to maxSaturationRate
            for (Vertiport vertiport : clusteredVertiports) {
                for (int i = 0; i < SIMULATION_HOURS; i++) {
                    vertiport.concurrentSaturationRates.put(i, vertiport.tempConcurrentSaturationRates.get(i));
                }
                vertiport.maxSaturationRate = vertiport.tempMaxSaturationRate;
            }

            List<Integer> bestSolutionID = new ArrayList<>(currentSolutionID);
            double bestEnergy = currentEnergey;
            double currentTemperature = optimizationConfiguration.initial_temperature;
            int notChangeCount = 0;

            int saturatedVertiportCount = 0;
            double maxSaturationRate = 0;


                for (Integer vertiportID : currentSolutionID) {
                    if (clusteredVertiports.get(vertiportID).maxSaturationRate > 1) {
                        saturatedVertiportCount++;
                    }
                    if (clusteredVertiports.get(vertiportID).maxSaturationRate > maxSaturationRate) {
                        maxSaturationRate = clusteredVertiports.get(vertiportID).maxSaturationRate;
                    }
                }
                log.info("Initial Solution: " + currentSolutionID + " Initial Energy: " + currentEnergey + " Saturated Vertiport Count: " + saturatedVertiportCount + " Max Saturation Rate: " + maxSaturationRate);
            for (int iteration = 0; iteration < optimizationConfiguration.max_iteration; iteration++) {

                // set the tempConcurrentSaturationRates and tempMaxSaturationRate for each vertiport as 0
                for (Vertiport vertiport : clusteredVertiports) {
                    for (int i = 0; i < SIMULATION_HOURS; i++) {
                        vertiport.tempConcurrentSaturationRates.put(i, 0.0);
                    }
                    vertiport.tempMaxSaturationRate = 0;
                }

                List<List<Integer>> newSolutionList = generateNewSolution(random, currentSolutionID, clusteredVertiports);
                List<Integer> newSolutionID = newSolutionList.get(0);



                double newEnergy = calculateFitness(newSolutionID, clusteredVertiports,tripItems,random,scenarioSpecific);
                double deltaEnergy = newEnergy - currentEnergey;

                if (deltaEnergy > 0 || Math.exp(deltaEnergy / currentTemperature) > random.nextDouble()) {
                    currentSolutionID = new ArrayList<>(newSolutionID);
                    currentEnergey = newEnergy;
                    saturatedVertiportCount = 0;
                    maxSaturationRate = 0;
                    // update the savedGeneralizedCost for each trip
                    for (TripItemForOptimization tripItemForOptimization : tripItems) {
                        if (!tripItemForOptimization.tempSavedGeneralizedCosts.isEmpty()) {
                            tripItemForOptimization.savedGeneralizedCost = tripItemForOptimization.tempSavedGeneralizedCosts.get(0);
                        }
                    }
                    // copy the tempConcurrentSaturationRates to concurrentSaturationRates and tempMaxSaturationRate to maxSaturationRate
                    for (Vertiport vertiport : clusteredVertiports) {
                        for (int i = 0; i < SIMULATION_HOURS; i++) {
                            vertiport.concurrentSaturationRates.put(i, vertiport.tempConcurrentSaturationRates.get(i));
                        }
                        vertiport.maxSaturationRate = vertiport.tempMaxSaturationRate;
                        if(vertiport.maxSaturationRate>1){
                           saturatedVertiportCount++;}
                        if (vertiport.maxSaturationRate > maxSaturationRate) {
                            maxSaturationRate = vertiport.maxSaturationRate;
                        }
                    }
                }
                // update the best solution
                if (currentEnergey > bestEnergy) {
                    bestSolutionID = new ArrayList<>(currentSolutionID);
                    bestEnergy = currentEnergey;
                    notChangeCount = 0;
                } else {
                    notChangeCount++;
                }

                // update the temperature
                if (currentTemperature > optimizationConfiguration.final_temperature) {
                    currentTemperature = currentTemperature * optimizationConfiguration.annealing_rate;
                }
                // log the information for every 100 iterations
                if (iteration % 1 == 0) {
                    log.info("Iteration: " + iteration + " Current Temperature: " + currentTemperature + " Current Energy: " + currentEnergey + " Best Energy: " + bestEnergy + " Saturated Vertiport Count: " + saturatedVertiportCount + " Max Saturation Rate: " + maxSaturationRate);
                }
                // if the best solution is not updated for more than 1000 iterations, break
                if (notChangeCount > optimizationConfiguration.max_not_change_count) {
                    break;
                }

            }

            log.info("Best Solution: " + bestSolutionID + " Best Energy: " + bestEnergy);
            // record the end time
            long endTime = System.currentTimeMillis();
            log.info( " Time: " + (endTime - startTime) / 1000 + "s");


        // begin of second phase: select the optimal vertiports from each cluster
        log.info("Start the second phase: select the optimal vertiports from each cluster...");
        List<Vertiport> selectedClusteredVertiports = new ArrayList<>();
        // assign clusters in bestSolutionID to selectedClusteredVertiports
        for (int vertiportID : bestSolutionID) {
            selectedClusteredVertiports.add(clusteredVertiports.get(vertiportID));
        }
        List<Vertiport> finalSelectedVertiports = new ArrayList<>();
        double totalConstructionCost = 0;
        for (Vertiport selectedClusteredVertiport : selectedClusteredVertiports) {
             List<Vertiport> subSelectedVertiports = findMinimumCostVertiports(clusterResults.get(selectedClusteredVertiport.ID), selectedClusteredVertiport.capacity);
                finalSelectedVertiports.addAll(subSelectedVertiports);
                for (Vertiport vertiport : subSelectedVertiports) {
                    totalConstructionCost += vertiport.constructionCost;
                }
        }
        List<Integer> finalSelectedVertiportsID = finalSelectedVertiports.stream().map(vertiport -> vertiport.ID).collect(Collectors.toList());
        log.info("Finished the second phase: select the optimal vertiports from each cluster.");
        log.info("Final Selected Vertiports: " + finalSelectedVertiportsID);
        double finalScore = scenarioSpecific.beta_savedCost * bestEnergy + scenarioSpecific.beta_constructionCost * totalConstructionCost;
        log.info("Final Score of the Selected Vertiports: " + finalScore);
    }

        public static List<Vertiport> findMinimumCostVertiports(List<Vertiport> vertiports, int requiredCapacity) {
        // dynamic programming array, dp[i] stores the minimum cost to achieve at least i capacity
        int[] dp = new int[requiredCapacity + 1];
        for (int i = 1; i <= requiredCapacity; i++) {
            dp[i] = Integer.MAX_VALUE;  // initialize the minimum cost to achieve at least i capacity to be infinity
        }
        dp[0] = 0;  // the cost to achieve 0 capacity is 0

        // fill the dp array
        for (Vertiport vertiport : vertiports) {
            for (int capacity = requiredCapacity; capacity >= vertiport.capacity; capacity--) {
                if (dp[capacity - vertiport.capacity] != Integer.MAX_VALUE) {
                    dp[capacity] = (int) Math.min(dp[capacity], dp[capacity - vertiport.capacity] + vertiport.constructionCost);
                }
            }
        }

        // if the minimum cost to achieve the required capacity is infinity, return an empty list
        if (dp[requiredCapacity] == Integer.MAX_VALUE) {
            return new ArrayList<>();
        }

        // backtracking to find the Vertiport set that forms the minimum cost
        List<Vertiport> result = new ArrayList<>();
        int capacity = requiredCapacity;
        while (capacity > 0) {
            for (Vertiport vertiport : vertiports) {
                if (capacity >= vertiport.capacity && dp[capacity] == dp[capacity - vertiport.capacity] + vertiport.constructionCost) {
                    result.add(vertiport);
                    capacity -= vertiport.capacity;
                    break;
                }
            }
        }

        return result;
    }
    public static List<Integer> generateRandomSolution(Random random,List<Vertiport> VertiportCandidates, int numOfSelectedVertiports) {

      List<Integer> selectedVertiportsID = new ArrayList<>();
      while (selectedVertiportsID.size() < numOfSelectedVertiports) {
          int vertiportID = VertiportCandidates.get(random.nextInt(VertiportCandidates.size())).ID;
          if (!selectedVertiportsID.contains(vertiportID)) {
              selectedVertiportsID.add(vertiportID);
          }
      }
      return selectedVertiportsID;
    }

    public static double calculateFitness(List<Integer> chosenVertiportID, List<Vertiport> vertiportsCandidates, List<TripItemForOptimization> deserializedTripItemForOptimizations, Random random, ScenarioSpecific scenarioSpecific) throws Exception {
        ExecutorService executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
        List<Future<Double>> futures = new ArrayList<>();

        for (TripItemForOptimization tripItemForOptimization : deserializedTripItemForOptimizations) {
            Callable<Double> task = () -> {
                tripItemForOptimization.tempSavedGeneralizedCosts.clear();
                List<Vertiport> originNeighbourVertiports = findAvailableNeighbourVertiports(chosenVertiportID, tripItemForOptimization.originNeighborVertiportCandidates);
                List<Vertiport> destinationNeighbourVertiports = findAvailableNeighbourVertiports(chosenVertiportID, tripItemForOptimization.destinationNeighborVertiportCandidates);

                if (!originNeighbourVertiports.isEmpty() && !destinationNeighbourVertiports.isEmpty()) {
                    if (originNeighbourVertiports.size() > 1 || destinationNeighbourVertiports.size() > 1 || originNeighbourVertiports.get(0).ID != destinationNeighbourVertiports.get(0).ID) {
                        tripItemForOptimization.isUAMAvailable = true;
                        tripItemForOptimization.originNeighborVertiports = originNeighbourVertiports;
                        tripItemForOptimization.destinationNeighborVertiports = destinationNeighbourVertiports;
                        calculateTripSavedCost(tripItemForOptimization, vertiportsCandidates, random, scenarioSpecific);
                        return tripItemForOptimization.tempSavedGeneralizedCosts.get(0);
                    }
                }

                tripItemForOptimization.isUAMAvailable = false;
                tripItemForOptimization.UAMUtilityVar = Double.NEGATIVE_INFINITY;
                tripItemForOptimization.tempSavedGeneralizedCosts.add(0.0);
                return 0.0;
            };
            futures.add(executor.submit(task));
        }

        double totalSavedGeneralizedCost = 0;
        for (Future<Double> future : futures) {
            totalSavedGeneralizedCost += future.get();
        }

        executor.shutdown();
        // Update the tempMaxSaturationRate for each vertiport
        for (int vertiportID : chosenVertiportID) {
            Vertiport vertiport = vertiportsCandidates.get(vertiportID);
            vertiport.tempMaxSaturationRate = 0;
            for (int i = 0; i < SIMULATION_HOURS; i++) {
                if (vertiport.concurrentSaturationRates.get(i) > vertiport.tempMaxSaturationRate) {
                    vertiport.tempMaxSaturationRate = vertiport.concurrentSaturationRates.get(i);
                }
            }
        }
        return totalSavedGeneralizedCost;
    }
    public static void calculateTripSavedCost(TripItemForOptimization tripItemForOptimization, List<Vertiport> vertiportsCandidates, Random random,ScenarioSpecific scenarioSpecific) {
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
                    double accessEmisson =0;
                    double egressEmisson =0;
                    double flightDistance= calculateEuciDistance(origin.coord,destination.coord);
                    double flightTime=flightDistance/flightSpeed+takeOffLandingTime;
                    double flightEmission=flightDistance/1000*UAM_EMISSION_FACTOR; // convert m to km
                    double flightCost=UAM_FIX_COST+ calculateEuciDistance(origin.coord,destination.coord)/1000*UAM_KM_COST;
                    double uamTravelTime=accessTime+egressTime+flightTime+UAM_PROCESS_TIME;
                    if (tripItemForOptimization.accessMode.equals("car") ){
                        accessCost=accessDistance/1000*CAR_COST;
                    }
                    if (tripItemForOptimization.egressMode.equals("car") ){
                        egressCost=egressDistance/1000*CAR_COST;
                    }
                    double UAMCost=accessCost+egressCost+flightCost;
                    double UAMGeneralizedCost=UAMCost+uamTravelTime* tripItemForOptimization.VOT+CARBON_EQUIVALENCE_FACTOR*(accessEmisson+egressEmisson+flightEmission);

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
            ModeDecider modeDecider=new ModeDecider(tripItemForOptimization.uamUtility,tripItemForOptimization.carUtility,tripItemForOptimization.ptUtility,random);
            Double [] modeSamples=modeDecider.sample(sampleSize);
            tripItemForOptimization.uamProbability=modeSamples[0];
            tripItemForOptimization.carProbability=modeSamples[1];
            tripItemForOptimization.ptProbability=modeSamples[2];

        // determine the probability of mode choice of each trip
            int arriveVertiportHour = (int) Math.floor((tripItemForOptimization.departureTime + tripItemForOptimization.accessTime) / 3600);
            int leaveVertiportHour = (int) Math.floor((tripItemForOptimization.departureTime + tripItemForOptimization.accessTime + UAM_PROCESS_TIME + tripItemForOptimization.flightTime) / 3600);
            vertiportsCandidates.get(tripItemForOptimization.accessVertiport.ID).tempConcurrentSaturationRates.put(arriveVertiportHour, vertiportsCandidates.get(tripItemForOptimization.accessVertiport.ID).tempConcurrentSaturationRates.get(arriveVertiportHour) + tripItemForOptimization.uamProbability / tripItemForOptimization.accessVertiport.capacity);
            vertiportsCandidates.get(tripItemForOptimization.egressVertiport.ID).tempConcurrentSaturationRates.put(leaveVertiportHour, vertiportsCandidates.get(tripItemForOptimization.egressVertiport.ID).tempConcurrentSaturationRates.get(leaveVertiportHour) + tripItemForOptimization.uamProbability / tripItemForOptimization.egressVertiport.capacity);



            objectiveFunctionBefore= tripItemForOptimization.currentGeneralizedCost+tripItemForOptimization.currentEmission*CARBON_EQUIVALENCE_FACTOR;
            objectiveFunctionAfter= modeSamples[1]*(tripItemForOptimization.carGeneralizedCost+ tripItemForOptimization.carEmission*CARBON_EQUIVALENCE_FACTOR)+modeSamples[2]*(tripItemForOptimization.ptGeneralizedCost+ tripItemForOptimization.ptEmission*CARBON_EQUIVALENCE_FACTOR)+modeSamples[0]* tripItemForOptimization.UAMGeneralizedCost;

        double savedGeneralizedCostOneTrip=objectiveFunctionBefore-objectiveFunctionAfter;
        if (savedGeneralizedCostOneTrip<0){
            savedGeneralizedCostOneTrip=0;
        }
        if (tripItemForOptimization.tripPurpose.startsWith("H") && considerReturnTrip){
            savedGeneralizedCostOneTrip=savedGeneralizedCostOneTrip*2;
        }

        // Include your logic here that was previously inside the loop over all trips
        tripItemForOptimization.tempSavedGeneralizedCosts.add(savedGeneralizedCostOneTrip);
    }

    public static List<List<Integer>> generateNewSolution(Random random, List<Integer> currentSolutionID, List<Vertiport> vertiportsCandidates) {
        List<Integer> newSolutionID = new ArrayList<>(currentSolutionID);
        List<Integer> differenceID = new ArrayList<>();
        List<Integer> notChosenVertiportID = new ArrayList<>();
        List<Integer> saturationVertiportsID = new ArrayList<>();
        List<Integer> notSaturationVertiportsID = new ArrayList<>();
        List<List<Integer>> result = new ArrayList<>();
        List<Vertiport> newSolution = new ArrayList<>();

        for (Integer vertiportID : currentSolutionID) {
            Vertiport vertiport = vertiportsCandidates.get(vertiportID);
            // calculate the max saturation rate of each vertiport
            if (vertiport.maxSaturationRate > 1) {
                saturationVertiportsID.add(vertiport.ID);
            } else {
                notSaturationVertiportsID.add(vertiport.ID); // This should be outside the loop
            }
        }

        for (Vertiport vertiport : vertiportsCandidates) {
            if (!currentSolutionID.contains(vertiport.ID)) {
                notChosenVertiportID.add(vertiport.ID);
            }
        }

        while (!saturationVertiportsID.isEmpty()) {
            Integer vertiportWithHighestSaturationRateID = selectHighestSaturationVertiport(saturationVertiportsID, vertiportsCandidates);
            List<Integer> neighborsID = getNeighborsID(vertiportsCandidates.get(vertiportWithHighestSaturationRateID));

            if (!neighborsID.isEmpty() && !currentSolutionID.containsAll(neighborsID)) { // Check if the neighbors are in the current solution or there is no neighbor
                Integer newVertiportID = randomSelectExcluding(neighborsID, currentSolutionID, random);
                Integer removedVertiportID = currentSolutionID.get(random.nextInt(currentSolutionID.size()));
                newSolutionID.remove(removedVertiportID);
                newSolutionID.add(newVertiportID);
                differenceID.add(removedVertiportID);
                differenceID.add(newVertiportID);
                break; // exit the loop once a successful replacement has been made
            } else {
                saturationVertiportsID.remove(vertiportWithHighestSaturationRateID); // Try next vertiport
            }
        }

        if (differenceID.isEmpty()) {
            // Fallback if no replacement was done
            Integer newVertiportID = notChosenVertiportID.get(random.nextInt(notChosenVertiportID.size()));
            Integer removedVertiportID = currentSolutionID.get(random.nextInt(currentSolutionID.size()));
            differenceID.add(removedVertiportID);
            differenceID.add(newVertiportID);
            newSolutionID.remove(removedVertiportID);
            newSolutionID.add(newVertiportID);
        }

        result.add(newSolutionID);
        result.add(differenceID);
        return result;
    }

    private static Integer selectHighestSaturationVertiport(List<Integer> saturationIDs, List<Vertiport> candidates) {
        Integer highestID = saturationIDs.get(0);
        for (Integer id : saturationIDs) {
            if (candidates.get(id).maxSaturationRate > candidates.get(highestID).maxSaturationRate) {
                highestID = id;
            }
        }
        return highestID;
    }

    private static List<Integer> getNeighborsID(Vertiport vertiport) {
        return vertiport.neighbors.stream().map(v -> v.ID).collect(Collectors.toList());
    }

    private static Integer randomSelectExcluding(List<Integer> options, List<Integer> exclude, Random random) {
        List<Integer> validOptions = options.stream().filter(opt -> !exclude.contains(opt)).collect(Collectors.toList());
        return validOptions.get(random.nextInt(validOptions.size()));
    }

    public static double calculateEuciDistance(Coord coord1, Coord coord2) {
        double euciDistance = Math.sqrt(Math.pow(coord1.getX() - coord2.getX(), 2) + Math.pow(coord1.getY() - coord2.getY(), 2));
        return euciDistance;
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


    public static int[] clusterDataDBSCAN(List<Vertiport> vertiports, double epsilon, int minPts) {
        double[][] locations = vertiports.stream()
                .map(v -> new double[]{v.coord.getX(), v.coord.getY()})
                .toArray(double[][]::new);
        DBSCAN<double[]> dbscan = DBSCAN.fit(locations, minPts, epsilon);
        return dbscan.y;
    }

    public static int[] clusterDataKMeans(List<Vertiport> vertiports, int k, int maxIterations, double tolerance) {
        double[][] locations = vertiports.stream()
                .map(v -> new double[]{v.coord.getX(), v.coord.getY()})
                .toArray(double[][]::new);


        KMeans kmeans = KMeans.fit(locations,k, maxIterations, tolerance);
        return kmeans.y;
    }
    public static List<Vertiport> summarizeClusters(List<Vertiport> vertiports, int[] assignments) {
        Map<Integer, List<Vertiport>> clusters = new HashMap<>();

        // group the vertiports by cluster
        for (int i = 0; i < assignments.length; i++) {
            int clusterId = assignments[i];
            Vertiport vertiport = vertiports.get(i);
            clusters.computeIfAbsent(clusterId, k -> new ArrayList<>()).add(vertiport);
        }

        // compute the capacity and average construction cost for each cluster
        List<Vertiport> summarizedClusters = new ArrayList<>();
        clusters.forEach((clusterId, clusterVerts) -> {
            int totalCapacity = 0;
            double weightedCost = 0;
            double sumX = 0;
            double sumY = 0;
            for (Vertiport v : clusterVerts) {
                totalCapacity += v.capacity;
                weightedCost += v.constructionCost * v.capacity;
                sumX += v.coord.getX();
                sumY += v.coord.getY();
            }
            double averageCost = totalCapacity > 0 ? weightedCost / totalCapacity : 0;
            double centroidX = sumX / clusterVerts.size();
            double centroidY = sumY / clusterVerts.size();

            Vertiport summarizedVertiport = new Vertiport();
            summarizedVertiport.ID = clusterId;
            summarizedVertiport.capacity = totalCapacity;
            summarizedVertiport.constructionCost = averageCost;
            summarizedVertiport.coord = new Coord(centroidX, centroidY);
            summarizedVertiport.concurrentSaturationRates = new ConcurrentHashMap<>();
            summarizedVertiport.tempConcurrentSaturationRates = new ConcurrentHashMap<>();
            for (int i = 0; i < SIMULATION_HOURS; i++) {
                summarizedVertiport.concurrentSaturationRates.put(i, 0.0);
                summarizedVertiport.tempConcurrentSaturationRates.put(i, 0.0);
            }
            summarizedVertiport.maxSaturationRate = 0;
            summarizedVertiport.tempMaxSaturationRate = 0;
            summarizedVertiport.neighbors = new ArrayList<>();
            summarizedClusters.add(summarizedVertiport);
        });

        return summarizedClusters;
    }

    public static Map<Integer, List<Vertiport>> saveClusterResults(List<Vertiport> vertiports, int[] assignments) {
        Map<Integer, List<Vertiport>> clusters = new HashMap<>();
        for (int i = 0; i < assignments.length; i++) {
            int clusterId = assignments[i];
            List<Vertiport> cluster = clusters.computeIfAbsent(clusterId, k -> new ArrayList<>());
            cluster.add(vertiports.get(i));
        }
        return clusters;
    }

    public static void findVertiportsNeighbours (List<Vertiport> vertiports, double NEIGHBOR_DISTANCE) {
        for (int i = 0; i < vertiports.size(); i++) {
            for (int j = 0; j < vertiports.size(); j++) {
                if (i != j) {
                    if (calculateEuciDistance(vertiports.get(i).coord, vertiports.get(j).coord) < NEIGHBOR_DISTANCE) {
                        vertiports.get(i).neighbors.add(vertiports.get(j));
                    }
                }
            }
        }
    }
}
