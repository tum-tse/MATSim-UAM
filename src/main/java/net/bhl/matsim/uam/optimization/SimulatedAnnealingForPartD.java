package net.bhl.matsim.uam.optimization;

import com.opencsv.CSVWriter;
import net.bhl.matsim.uam.optimization.utils.*;
import org.apache.log4j.Logger;
import org.matsim.api.core.v01.Coord;
import org.matsim.core.config.Config;
import org.matsim.core.config.ConfigUtils;
import org.matsim.utils.MemoryObserver;
import smile.clustering.DBSCAN;
import smile.clustering.KMeans;

import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.*;
import java.util.concurrent.*;
import java.util.stream.Collectors;

public class SimulatedAnnealingForPartD {
    // This class provides the simulated annealing algorithm to solve the vertiport siting problem.

    public static final Logger log = Logger.getLogger(SimulatedAnnealingForPartD.class);
    private static final int MEMORY_CHECK_INTERVAL=600;
    public static String tripItemFile;
    public static String vertiportCandidateFile;
    public static double flightSpeed; // m/s
    public static double UAM_PROCESS_TIME; // s
    public static double takeOffLandingTime; // s
    public static int NUM_OF_SELECTED_VERTIPORTS; // In hierarchical clustering, the number of selected vertiports is the number of clusters
    private static double UAM_FIX_COST;
    private static double UAM_KM_COST;
    private static double CAR_COST; // Euro/km
    public static double PT_COST; // Euro/trip
    private static boolean considerReturnTrip;
    private static  int sampleSize;
    private static long RANDOM_SEED;
    public static double car_utility_mean;
    public static double car_utility_sigma;
    public static double pt_utility_mean;
    public static double pt_utility_sigma;
    public static double CAR_EMISSION_FACTOR; // kg/km
    public static double PT_EMISSION_FACTOR; // kg/km
    public static double UAM_EMISSION_FACTOR; // kg/km
    public static double CARBON_EQUIVALENCE_FACTOR; // Euro/kgCO2
    private static int SIMULATION_HOURS;
    public static String configPath;
    public static double NEIGHBOR_DISTANCE;
    public static double beta_savedCost;
    public static double beta_savedEmission;
    public static double beta_constructionCost;
    public static String scenarioConfigurations;

    public static boolean incrementalSiting; // For incremental siting
    public static String existingVertiportFile; // For incremental siting
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


        // load the scenario specific parameters
        NUM_OF_SELECTED_VERTIPORTS = scenarioSpecific.num_of_selected_vertiports;
        UAM_FIX_COST = scenarioSpecific.uam_fix_cost;
        UAM_KM_COST = scenarioSpecific.uam_km_cost;
        CAR_COST = scenarioSpecific.car_km_cost;
        PT_COST = scenarioSpecific.pt_cost;
        flightSpeed = scenarioSpecific.flight_speed;
        UAM_PROCESS_TIME = scenarioSpecific.uam_process_time;
        takeOffLandingTime = scenarioSpecific.uam_take_off_landing_time;
        considerReturnTrip = scenarioSpecific.consider_return_trip;
        CAR_EMISSION_FACTOR = scenarioSpecific.car_emission_factor;
        PT_EMISSION_FACTOR = scenarioSpecific.pt_emission_factor;
        UAM_EMISSION_FACTOR = scenarioSpecific.uam_emission_factor;
        SIMULATION_HOURS = scenarioSpecific.simulation_hours;
        CARBON_EQUIVALENCE_FACTOR = scenarioSpecific.carbon_equivalent_cost;
        RANDOM_SEED =  scenarioSpecific.random_seed;
        sampleSize = scenarioSpecific.sampleSize;
        NEIGHBOR_DISTANCE = scenarioSpecific.neighbouring_distance;
        beta_savedEmission= scenarioSpecific.beta_savedEmission;
        beta_savedCost= scenarioSpecific.beta_savedCost;
        beta_constructionCost= scenarioSpecific.beta_constructionCost;
        car_utility_mean= scenarioSpecific.car_utility_mean;
        car_utility_sigma= scenarioSpecific.car_utility_sigma;
        pt_utility_mean= scenarioSpecific.pt_utility_mean;
        pt_utility_sigma= scenarioSpecific.pt_utility_sigma;
        incrementalSiting = scenarioSpecific.incrementaSiting;
        existingVertiportFile= scenarioSpecific.existingVertiportFile;


        // load the config
        Config config = ConfigUtils.loadConfig(configPath);
        log.info("Loading the vertiport candidates...");
        VertiportReader vertiportReader = new VertiportReader();
        List<Vertiport> vertiportsCandidates = vertiportReader.getVertiportsWithNeighbors(vertiportCandidateFile);
        HashMap<Integer,Vertiport> vertiportsCandidatesMap = new HashMap<>();
        for (Vertiport vertiport : vertiportsCandidates) {
            vertiportsCandidatesMap.put(vertiport.ID, vertiport);
        }
        List<Vertiport> existingVertiports = new ArrayList<>();
        HashMap<Integer,Vertiport> existingVertiportsMap = new HashMap<>();

        log.info("Finished loading the vertiport candidates.");

        // In incremental siting scenarios, the already selected vertiports are loaded
        if (incrementalSiting) {
        log.info("Loading the existing vertiports...");
        existingVertiports = vertiportReader.getVertiportsWithNeighbors(existingVertiportFile);
        for (Vertiport vertiport : existingVertiports) {
                existingVertiportsMap.put(vertiport.ID, vertiport);
            }
        log.info("Finished loading the already selected vertiports."); }

        // Cluster the vertiport candidates with K-means clustering algorithm
        log.info("Start clustering the vertiport candidates & already selected vertiports...");
        //int[] clusterAssigment = clusterDataDBSCAN(vertiportsCandidates, EPSILON, MIN_PTS);
        int[] clusterAssigment = clusterDataKMeans(vertiportsCandidates, clusteringConfiguration.num_of_clusters, clusteringConfiguration.num_of_cluster_iterations, clusteringConfiguration.tolerance);
        List<Vertiport> clusteredVertiports = summarizeClusters(vertiportsCandidates, clusterAssigment);
        HashMap<Integer, Vertiport> clusteredVertiportsMap = new HashMap<>();
        for (Vertiport vertiport : clusteredVertiports) {
            clusteredVertiportsMap.put(vertiport.ID, vertiport);
        }
        Map<Integer, List<Vertiport>> clusterResults = saveClusterResults(vertiportsCandidates, clusterAssigment);
        // Write out the clustering results
        CSVWriter clusteringWriterForCandidates = new CSVWriter(new FileWriter(clusteringConfiguration.clustering_output_file_candidates),CSVWriter.DEFAULT_SEPARATOR, CSVWriter.NO_QUOTE_CHARACTER, CSVWriter.DEFAULT_ESCAPE_CHARACTER, CSVWriter.DEFAULT_LINE_END);
        String[] clusteringHeader = {"VertiportID", "CoordX","CoordY","ClusterID"};
        clusteringWriterForCandidates.writeNext(clusteringHeader);
        for (Vertiport vertiport : vertiportsCandidates) {
           for (Map.Entry<Integer, List<Vertiport>> entry : clusterResults.entrySet()) {
               if (entry.getValue().contains(vertiport)) {
                   String[] data = {String.valueOf(vertiport.ID), String.valueOf(vertiport.coord.getX()), String.valueOf(vertiport.coord.getY()), String.valueOf(entry.getKey())};
                   clusteringWriterForCandidates.writeNext(data);
               }
           }
        }
        clusteringWriterForCandidates.close();
        List<Vertiport> clusteredExistingAndCandidates = new ArrayList<>();
        List<Vertiport> allClusteredVertiports = new ArrayList<>(clusteredVertiports);
        Map<Integer, List<Vertiport>> clusterResultsExisting = new HashMap<>();
        if (incrementalSiting) {
            int[] clusterAssigmentExisting = clusterDataKMeans(existingVertiports, clusteringConfiguration.num_cluster_in_existing_vertiports, clusteringConfiguration.num_of_cluster_iterations, clusteringConfiguration.tolerance);
            clusteredExistingAndCandidates = summarizeClustersForAlreadySelectedVertiports(existingVertiports, clusterAssigmentExisting);
            allClusteredVertiports.addAll(clusteredExistingAndCandidates);
            clusterResultsExisting = saveClusterResultsForExistingVertiports(existingVertiports, clusterAssigmentExisting);
            // Write out the clustering results for existing vertiports
            CSVWriter clusteringWriterForExisting = new CSVWriter(new FileWriter(clusteringConfiguration.clustering_output_file_existing),CSVWriter.DEFAULT_SEPARATOR, CSVWriter.NO_QUOTE_CHARACTER, CSVWriter.DEFAULT_ESCAPE_CHARACTER, CSVWriter.DEFAULT_LINE_END);
            clusteringWriterForExisting.writeNext(clusteringHeader);
            for (Vertiport vertiport : existingVertiports) {
                for (Map.Entry<Integer, List<Vertiport>> entry : clusterResultsExisting.entrySet()) {
                    if (entry.getValue().contains(vertiport)) {
                        String[] data = {String.valueOf(vertiport.ID), String.valueOf(vertiport.coord.getX()), String.valueOf(vertiport.coord.getY()), String.valueOf(entry.getKey())};
                        clusteringWriterForExisting.writeNext(data);
                    }
                }
            }
            clusteringWriterForExisting.close();
        }
        findVertiportsNeighbours(allClusteredVertiports, NEIGHBOR_DISTANCE);
        log.info("Finished clustering the vertiport candidates.");

        log.info("Loading the trips...");
        TripItemReaderForOptimization tripItemReaderForOptimization = new TripItemReaderForOptimization(scenarioSpecific);
        List<TripItemForOptimization> tripItems = tripItemReaderForOptimization.getTripItemsForOptimization(tripItemFile);
        log.info("Finished loading the trips.");

        // Pre-calculate the access and egress time and distance for each vertiport candidate of each trip
        log.info("Start pre-calculating the access and egress time and distance for each vertiport candidate of each trip...");
        AccessEgressCostCalculator accessEgressCostCalculator = new AccessEgressCostCalculator(tripItems, allClusteredVertiports, config, scenarioSpecific);
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
            double currentEnergey = calculateFitness(currentSolutionID,  allClusteredVertiports, tripItems,clusterResults,random,scenarioSpecific);
            for (TripItemForOptimization tripItemForOptimization : tripItems) {
                if (!tripItemForOptimization.tempSavedGeneralizedCosts.isEmpty()) {
                    tripItemForOptimization.savedGeneralizedCost = tripItemForOptimization.tempSavedGeneralizedCosts.get(0);
                    tripItemForOptimization.savedEmission = tripItemForOptimization.tempSavedEmission.get(0);
                }
            }
            // copy the tempConcurrentSaturationRates to concurrentSaturationRates and tempMaxSaturationRate to maxSaturationRate
            for (Vertiport vertiport : clusteredVertiports) {
                for (int i = 0; i < SIMULATION_HOURS*4; i++) { // convert the simulation hours to 15 minutes
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
                HashMap<Integer,List<Double>> iterationRecord=new HashMap<>(); // record the indicator through the iterations


            for (int iteration = 0; iteration < optimizationConfiguration.max_iteration; iteration++) {

                // set the tempConcurrentSaturationRates and tempMaxSaturationRate for each vertiport as 0
                for (Vertiport vertiport : clusteredVertiports) {
                    for (int i = 0; i < SIMULATION_HOURS*4; i++) {
                        vertiport.tempConcurrentSaturationRates.put(i, 0.0);
                    }
                    vertiport.tempMaxSaturationRate = 0;
                }

                List<List<Integer>> newSolutionList;
                if (incrementalSiting) {
                    newSolutionList = generateNewSolution(random, currentSolutionID, clusteredVertiports, clusteredExistingAndCandidates);
                } else {
                    newSolutionList = generateNewSolution(random, currentSolutionID, clusteredVertiportsMap);
                }
                List<Integer> newSolutionID = newSolutionList.get(0);


                double newEnergy = 0;
                if (incrementalSiting){
                    newEnergy=calculateFitness(newSolutionID,allClusteredVertiports,tripItems,clusterResults,random,scenarioSpecific);}
                else{
                    newEnergy = calculateFitness(newSolutionID, clusteredVertiports, tripItems,clusterResults,random,scenarioSpecific);
                }
                double deltaEnergy = newEnergy - currentEnergey;

                if (deltaEnergy > 0 || Math.exp(deltaEnergy / currentTemperature) > random.nextDouble()) {
                    currentSolutionID = new ArrayList<>(newSolutionID);
                    currentEnergey = newEnergy;
                    saturatedVertiportCount = 0;
                    maxSaturationRate = 0;
                    // copy the temp values to the real values
                    for (TripItemForOptimization tripItemForOptimization : tripItems) {
                        if (!tripItemForOptimization.tempSavedGeneralizedCosts.isEmpty()) {
                            tripItemForOptimization.savedGeneralizedCost = tripItemForOptimization.tempSavedGeneralizedCosts.get(0);
                            tripItemForOptimization.savedEmission = tripItemForOptimization.tempSavedEmission.get(0);
                            tripItemForOptimization.savedTravelTime = tripItemForOptimization.tempSavedTravelTime.get(0);
                        }
                        tripItemForOptimization.carProbability = tripItemForOptimization.tempCarProbability;
                        tripItemForOptimization.ptProbability = tripItemForOptimization.tempPTProbability;
                        tripItemForOptimization.uamProbability = tripItemForOptimization.tempUAMProbability;
                        tripItemForOptimization.isUAMAvailable = tripItemForOptimization.tempIsUAMAvailable;
                    }
                    // copy the tempConcurrentSaturationRates to concurrentSaturationRates and tempMaxSaturationRate to maxSaturationRate
                    for (Vertiport vertiport : clusteredVertiports) {
                        for (int i = 0; i < SIMULATION_HOURS*4; i++) {
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
                if (iteration % 100 == 0) {
                    log.info("Iteration: " + iteration + " Current Temperature: " + currentTemperature + " Current Energy: " + currentEnergey + " Best Energy: " + bestEnergy + " Saturated Vertiport Count: " + saturatedVertiportCount + " Max Saturation Rate: " + maxSaturationRate);
                }
                iterationRecord.put(iteration, Arrays.asList(currentTemperature,bestEnergy,(double)saturatedVertiportCount, maxSaturationRate));
                // if the best solution is not updated for more than 1000 iterations, break
                if (notChangeCount > optimizationConfiguration.max_not_change_count) {
                    break;
                }

            }

            log.info("Best Solution: " + bestSolutionID + " Best Energy: " + bestEnergy);
            // record the end time
            long endTime = System.currentTimeMillis();
            log.info( " Time: " + (endTime - startTime) / 1000 + "s");
            // Write the iteration record to a csv file
            CSVWriter iterationWriter = new CSVWriter(new FileWriter(optimizationConfiguration.iteration_record_file),CSVWriter.DEFAULT_SEPARATOR, CSVWriter.NO_QUOTE_CHARACTER, CSVWriter.DEFAULT_ESCAPE_CHARACTER, CSVWriter.DEFAULT_LINE_END);
            String[] iterationHeader = {"Iteration", "Temperature", "Best Energy", "Saturated Vertiport Count", "Max Saturation Rate"};
            iterationWriter.writeNext(iterationHeader);
            for (Map.Entry<Integer, List<Double>> entry : iterationRecord.entrySet()) {
                String[] data = {String.valueOf(entry.getKey()), String.valueOf(entry.getValue().get(0)), String.valueOf(entry.getValue().get(1)), String.valueOf(entry.getValue().get(2)), String.valueOf(entry.getValue().get(3))};
                iterationWriter.writeNext(data);
            }
            iterationWriter.close();


            List<Vertiport> selectedClusteredVertiports = new ArrayList<>();
        // assign clusters in bestSolutionID to selectedClusteredVertiports
        for (int vertiportID : bestSolutionID) {
            selectedClusteredVertiports.add(clusteredVertiportsMap.get(vertiportID));
        }
        List<Vertiport> finalSelectedVertiports = new ArrayList<>();
        double totalConstructionCost = 0;
        HashMap<Integer,List<Integer>> requiredAndAchievedCapacityMapForCandidates = new HashMap<>();
        for (Vertiport selectedClusteredVertiport : selectedClusteredVertiports) {
                HashMap<List<Vertiport>,HashMap<Integer,Double>> subSelectedVertiportsAndCost = findMinimumCostVertiports(clusterResults.get(selectedClusteredVertiport.ID), (int) (selectedClusteredVertiport.capacity*selectedClusteredVertiport.maxSaturationRate)+1);
                finalSelectedVertiports.addAll(subSelectedVertiportsAndCost.entrySet().iterator().next().getKey());
                totalConstructionCost += subSelectedVertiportsAndCost.entrySet().iterator().next().getValue().values().iterator().next();
                List<Integer> requiredAndAchievedCapacity = new ArrayList<>();
                requiredAndAchievedCapacity.add((int) (selectedClusteredVertiport.capacity*selectedClusteredVertiport.maxSaturationRate)+1);
                requiredAndAchievedCapacity.add(subSelectedVertiportsAndCost.entrySet().iterator().next().getValue().entrySet().iterator().next().getKey());
                requiredAndAchievedCapacityMapForCandidates.put(selectedClusteredVertiport.ID, requiredAndAchievedCapacity);
        }
        HashMap<Integer,List<Integer>> requiredAndAchievedCapacityMapForExisting = new HashMap<>();
        if(incrementalSiting){
            for (Vertiport selectedClusteredVertiport : existingVertiports) {
                int requiredCapacity = (int) (selectedClusteredVertiport.capacity * selectedClusteredVertiport.maxSaturationRate)+1;
                int achievedCapacity = selectedClusteredVertiport.capacity;
                requiredAndAchievedCapacityMapForExisting.put(selectedClusteredVertiport.ID, new ArrayList<>(Arrays.asList(requiredCapacity, achievedCapacity)));
            }
        }

        // Write out the required and achieved capacity for each vertiport cluster
        CSVWriter requiredAndAchievedWriter = new CSVWriter(new FileWriter(scenarioSpecific.outputVertiportBasedIndicatorFile),CSVWriter.DEFAULT_SEPARATOR, CSVWriter.NO_QUOTE_CHARACTER, CSVWriter.DEFAULT_ESCAPE_CHARACTER, CSVWriter.DEFAULT_LINE_END);
        String[] requiredAndAchievedHeader = {"ClusterID", "RequiredCapacity", "AchievedCapacity"};
        requiredAndAchievedWriter.writeNext(requiredAndAchievedHeader);
        for (Map.Entry<Integer, List<Integer>> entry : requiredAndAchievedCapacityMapForCandidates.entrySet()) {
            String[] data = {String.valueOf(entry.getKey()), String.valueOf(entry.getValue().get(0)), String.valueOf(entry.getValue().get(1))};
            requiredAndAchievedWriter.writeNext(data);
        }
        if (incrementalSiting) {
            for (Map.Entry<Integer, List<Integer>> entry : requiredAndAchievedCapacityMapForExisting.entrySet()) {
                String[] data = {String.valueOf(entry.getKey()), String.valueOf(entry.getValue().get(0)), String.valueOf(entry.getValue().get(1))};
                requiredAndAchievedWriter.writeNext(data);
            }
        }
        requiredAndAchievedWriter.close();

        List<Integer> finalSelectedVertiportsID = finalSelectedVertiports.stream().map(vertiport -> vertiport.ID).collect(Collectors.toList());
        log.info("Finished the second phase: select the optimal vertiports from each cluster.");
        log.info("Final Selected Vertiports: " + finalSelectedVertiportsID);
        double finalScore = bestEnergy + beta_constructionCost * totalConstructionCost;
        log.info("Final Score of the Selected Vertiports: " + finalScore);
        // save the selected vertiports as a csv file
        CSVWriter vertiportWriter = new CSVWriter(new FileWriter(scenarioSpecific.outputVertiportFile),CSVWriter.DEFAULT_SEPARATOR, CSVWriter.NO_QUOTE_CHARACTER, CSVWriter.DEFAULT_ESCAPE_CHARACTER, CSVWriter.DEFAULT_LINE_END);
        String[] vertiportHeader = {"ID", "coordX", "coordY", "capacity","constructionCost","ClusterID"};
        vertiportWriter.writeNext(vertiportHeader);
        for (Vertiport vertiport : finalSelectedVertiports) {
            // Find the corresponding cluster ID
            for (Map.Entry<Integer, List<Vertiport>> entry : clusterResults.entrySet()) {
                if (entry.getValue().contains(vertiport)) {
                    String[] data = {String.valueOf(vertiport.ID), String.valueOf(vertiport.coord.getX()), String.valueOf(vertiport.coord.getY()), String.valueOf(vertiport.capacity), String.valueOf(vertiport.constructionCost), String.valueOf(entry.getKey())};
                    vertiportWriter.writeNext(data);
                }
            }
        }
        // if incremental siting, save the already selected vertiports in the same file
        if (incrementalSiting) {
            vertiportWriter.writeNext(new String[]{}); // add an empty line
            for (Vertiport vertiport : existingVertiports) {
                for (Map.Entry<Integer, List<Vertiport>> entry : clusterResultsExisting.entrySet()) {
                    if (entry.getValue().contains(vertiport)) {
                        String[] data = {String.valueOf(vertiport.ID), String.valueOf(vertiport.coord.getX()), String.valueOf(vertiport.coord.getY()), String.valueOf(vertiport.capacity), String.valueOf(vertiport.constructionCost), String.valueOf(entry.getKey())};
                        vertiportWriter.writeNext(data);
                    }
                }
            }
        }
        vertiportWriter.close();
        log.info("Finished saving the selected vertiports as a csv file.");
        log.info("Write the indicators for analysis...");
        writeTripBasedAnalysisResult(scenarioSpecific.outputTripBasedIndicatorFile,tripItems);
    }

    public static HashMap<List<Vertiport>, HashMap<Integer,Double>> findMinimumCostVertiports(List<Vertiport> vertiports, int requiredCapacity) {
        int totalCapacity = vertiports.stream().mapToInt(v -> v.capacity).sum();

        // If the total capacity is less than the required capacity, return all vertiports
        if (requiredCapacity > totalCapacity) {
            HashMap<List<Vertiport>, HashMap<Integer,Double>> result = new HashMap<>();
            HashMap<Integer,Double> capacityAndCost = new HashMap<>();
            capacityAndCost.put(totalCapacity, vertiports.stream().mapToDouble(v -> v.constructionCost).sum());
            result.put(vertiports, capacityAndCost);
            return result;
        }

        // Dynamic programming to find the minimum cost to reach at least required capacity
        double[] dp = new double[totalCapacity + 1];
        for (int i = 1; i <= totalCapacity; i++) {
            dp[i] = Double.MAX_VALUE;  // initialize the dp array with the maximum value
        }
        dp[0] = 0;  // the minimum cost to reach 0 capacity is 0

        // Fill in the dp array
        for (Vertiport vertiport : vertiports) {
            for (int capacity = totalCapacity; capacity >= vertiport.capacity; capacity--) {
                if (dp[capacity - vertiport.capacity] != Double.MAX_VALUE) {
                    dp[capacity] = Math.min(dp[capacity], dp[capacity - vertiport.capacity] + vertiport.constructionCost);
                }
            }
        }

        // find the minimum capacity that can be achieved
        int minCapacityAchieved = requiredCapacity;
        while (minCapacityAchieved <= totalCapacity && dp[minCapacityAchieved] == Double.MAX_VALUE) {
            minCapacityAchieved++;
        }

        // backtracking to find the selected vertiports
        List<Vertiport> selectedVertiports = new ArrayList<>();
        boolean[] used = new boolean[vertiports.size()]; // tracking which vertiports have been used
        int capacity = minCapacityAchieved;
        while (capacity > 0) {
            for (int i = 0; i < vertiports.size(); i++) {
                Vertiport vertiport = vertiports.get(i);
                if (!used[i] && capacity >= vertiport.capacity && dp[capacity] == dp[capacity - vertiport.capacity] + vertiport.constructionCost) {
                    selectedVertiports.add(vertiport);
                    capacity -= vertiport.capacity;
                    used[i] = true; // mark this vertiport as used
                    break;
                }
            }
        }

        HashMap<List<Vertiport>, HashMap<Integer,Double>> result = new HashMap<>();
        HashMap<Integer,Double> capacityAndCost = new HashMap<>();
        capacityAndCost.put(minCapacityAchieved, dp[minCapacityAchieved]);
        result.put(selectedVertiports, capacityAndCost);

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

    public static double calculateFitness(List<Integer> chosenVertiportID, List<Vertiport> vertiportsCandidates, List<TripItemForOptimization> deserializedTripItemForOptimizations, Map<Integer, List<Vertiport>> savedClusterResults, Random random, ScenarioSpecific scenarioSpecific) throws Exception {
        ExecutorService executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
        List<Future<Double>> futures = new ArrayList<>();
        HashMap<Integer,Vertiport> vertiportsCandidatesMap = new HashMap<>();
        for (Vertiport vertiport : vertiportsCandidates) {
            vertiportsCandidatesMap.put(vertiport.ID, vertiport);
        }
        for (TripItemForOptimization tripItemForOptimization : deserializedTripItemForOptimizations) {
            Callable<Double> task = () -> {
                tripItemForOptimization.tempSavedGeneralizedCosts.clear();
                tripItemForOptimization.tempSavedEmission.clear();
                tripItemForOptimization.tempSavedTravelTime.clear();
                List<Vertiport> originNeighbourVertiports = findAvailableNeighbourVertiports(chosenVertiportID, tripItemForOptimization.originNeighborVertiportCandidates);
                List<Vertiport> destinationNeighbourVertiports = findAvailableNeighbourVertiports(chosenVertiportID, tripItemForOptimization.destinationNeighborVertiportCandidates);

                if (!originNeighbourVertiports.isEmpty() && !destinationNeighbourVertiports.isEmpty()) {
                    if (originNeighbourVertiports.size() > 1 || destinationNeighbourVertiports.size() > 1 || originNeighbourVertiports.get(0).ID != destinationNeighbourVertiports.get(0).ID) {
                        tripItemForOptimization.tempIsUAMAvailable = true;
                        tripItemForOptimization.originNeighborVertiports = originNeighbourVertiports;
                        tripItemForOptimization.destinationNeighborVertiports = destinationNeighbourVertiports;
                        calculateTripSavedCost(tripItemForOptimization, vertiportsCandidatesMap, random, scenarioSpecific);
                        return beta_savedCost+tripItemForOptimization.tempSavedGeneralizedCosts.get(0)+beta_savedEmission*tripItemForOptimization.tempSavedEmission.get(0)*CARBON_EQUIVALENCE_FACTOR;
                    }
                }

                tripItemForOptimization.tempIsUAMAvailable = false;
                tripItemForOptimization.UAMUtilityVar = Double.NEGATIVE_INFINITY;
                tripItemForOptimization.tempSavedGeneralizedCosts.add(0.0);
                tripItemForOptimization.tempSavedEmission.add(0.0);
                tripItemForOptimization.tempSavedTravelTime.add(0.0);
                tripItemForOptimization.tempCarProbability=tripItemForOptimization.carProbabilityBefore;
                tripItemForOptimization.tempPTProbability=tripItemForOptimization.ptProbabilityBefore;
                tripItemForOptimization.tempUAMProbability=0;
                return 0.0;
            };
            futures.add(executor.submit(task));
        }

        double totalSGAndSE = 0;
        for (Future<Double> future : futures) {
            totalSGAndSE += future.get();
        }

        executor.shutdown();
        double totalConstructionCost=0;
        List<Vertiport> currentSelectedVertiportUnits = new ArrayList<>();
        // Update the tempMaxSaturationRate for each vertiport
        for (int vertiportID : chosenVertiportID) {
            Vertiport vertiport = vertiportsCandidatesMap.get(vertiportID);
            vertiport.tempMaxSaturationRate = 0;
            for (int i = 0; i < SIMULATION_HOURS*4; i++) {
                if (vertiport.tempConcurrentSaturationRates.get(i) > vertiport.tempMaxSaturationRate) {
                    vertiport.tempMaxSaturationRate = vertiport.tempConcurrentSaturationRates.get(i);
                }
            }
            double requiredCapacity=vertiport.tempMaxSaturationRate*vertiport.capacity;
            List<Vertiport> vertiportUnits = savedClusterResults.get(vertiportID);
            HashMap<List<Vertiport>,HashMap<Integer,Double>> result=findMinimumCostVertiports(vertiportUnits, (int) requiredCapacity+1);
            totalConstructionCost+=result.values().iterator().next().values().iterator().next();
            currentSelectedVertiportUnits.addAll(result.keySet().iterator().next());
        }
        double totalSGAndSEAndCC=totalSGAndSE+beta_constructionCost*totalConstructionCost;
        return totalSGAndSEAndCC;
    }
    public static void calculateTripSavedCost(TripItemForOptimization tripItemForOptimization, HashMap<Integer,Vertiport> vertiportsCandidatesMap, Random random,ScenarioSpecific scenarioSpecific) {
        for (Vertiport vertiport : tripItemForOptimization.originNeighborVertiports) {
            tripItemForOptimization.originNeighborVertiportsTimeAndDistance.put(vertiport, tripItemForOptimization.originNeighborVertiportCandidatesTimeAndDistance.get(vertiport));
        }
        for (Vertiport vertiport : tripItemForOptimization.destinationNeighborVertiports) {
            tripItemForOptimization.destinationNeighborVertiportsTimeAndDistance.put(vertiport, tripItemForOptimization.destinationNeighborVertiportCandidatesTimeAndDistance.get(vertiport));
        }

        double generalizedTravelCostAfter;
        double emissionAfter;
        double travelTimeAfter;
        // Initialize the Vertiport Allocation
        // Find the vertiport pair for a trip with the lowest uam generalized cost, the access and ergress vertiport could not be the same
        double lowestUAMGCAndEmission = Double.MAX_VALUE;
        for (Vertiport origin : tripItemForOptimization.originNeighborVertiports) {
            for (Vertiport destination :  tripItemForOptimization.destinationNeighborVertiports) {
                if (origin.ID != destination.ID) {
                    double accessTime = tripItemForOptimization.originNeighborVertiportsTimeAndDistance.get(origin).get("travelTime");
                    double egressTime = tripItemForOptimization.destinationNeighborVertiportsTimeAndDistance.get(destination).get("travelTime");
                    double accessDistance = tripItemForOptimization.originNeighborVertiportsTimeAndDistance.get(origin).get("distance");
                    double egressDistance = tripItemForOptimization.destinationNeighborVertiportsTimeAndDistance.get(destination).get("distance");
                    double accessMode = tripItemForOptimization.originNeighborVertiportsTimeAndDistance.get(origin).get("accessMode"); // 0: walk, 1: car, 2: pt
                    double egressMode = tripItemForOptimization.destinationNeighborVertiportsTimeAndDistance.get(destination).get("egressMode");
                    double accessCost =0;
                    double egressCost =0;
                    double accessEmisson =0;
                    double egressEmisson =0;
                    double flightDistance= calculateEuciDistance(origin.coord,destination.coord);
                    double flightTime=flightDistance/flightSpeed+takeOffLandingTime;
                    double flightEmission=flightDistance/1000*UAM_EMISSION_FACTOR; // convert m to km
                    double flightCost=UAM_FIX_COST+ calculateEuciDistance(origin.coord,destination.coord)/1000*UAM_KM_COST;
                    double uamTravelTime=accessTime+egressTime+flightTime+UAM_PROCESS_TIME;
                    if (accessMode==1){
                        accessCost=accessDistance/1000*CAR_COST;
                        accessEmisson=accessDistance/1000*CAR_EMISSION_FACTOR;
                    }
                    if (accessMode ==2) {
                        accessCost= PT_COST;
                        accessEmisson=accessDistance/1000*PT_EMISSION_FACTOR;
                    }
                    if (egressMode==1 ){
                        egressCost=egressDistance/1000*CAR_COST;
                        egressEmisson=egressDistance/1000*CAR_EMISSION_FACTOR;
                    }
                    if (egressMode==2){
                        egressCost= PT_COST;
                        egressEmisson=egressDistance/1000*PT_EMISSION_FACTOR;
                    }
                    double UAMTotalCost=accessCost+egressCost+flightCost;
                    double UAMTotalGC= UAMTotalCost+uamTravelTime* tripItemForOptimization.VOT;
                    double UAMTotalEmission=accessEmisson+egressEmisson+flightEmission;
                    double UAMTotalGCAndEmission=beta_savedCost*UAMTotalGC+beta_savedEmission*UAMTotalEmission*CARBON_EQUIVALENCE_FACTOR; // Already convert to monetary unit

                    if (UAMTotalGCAndEmission < lowestUAMGCAndEmission) {
                        tripItemForOptimization.uamTravelTime=uamTravelTime;
                        tripItemForOptimization.UAMCost=UAMTotalCost;
                        tripItemForOptimization.uamEmission=UAMTotalEmission;
                        tripItemForOptimization.UAMGeneralizedCost=UAMTotalGC;
                        tripItemForOptimization.UAMUtilityVar= scenarioSpecific.calculateUAMUtilityVAR(flightTime,uamTravelTime-flightTime,UAMTotalCost);
                        tripItemForOptimization.uamUtility= tripItemForOptimization.UAMUtilityFix+ tripItemForOptimization.UAMUtilityVar;
                        tripItemForOptimization.accessVertiport = origin;
                        tripItemForOptimization.egressVertiport = destination;
                        lowestUAMGCAndEmission = UAMTotalGCAndEmission;
                    }
                }
            }
        }
            // Determine the probability of each mode by Monte Carlo Sampling
            // Create a double list to store the probability of each mode
            ModeDecider modeDecider=new ModeDecider(tripItemForOptimization.uamUtility,tripItemForOptimization.carUtility,tripItemForOptimization.ptUtility,car_utility_mean,car_utility_sigma,pt_utility_mean,pt_utility_sigma,random);
            Double [] modeSamples=modeDecider.sampleWithErrorTerm(sampleSize);
            tripItemForOptimization.tempUAMProbability=modeSamples[0];
            tripItemForOptimization.tempCarProbability=modeSamples[1];
            tripItemForOptimization.tempPTProbability=modeSamples[2];


            int arriveVertiportTimeWindow = (int) Math.floor((tripItemForOptimization.departureTime + tripItemForOptimization.accessTime) / 3600*4);
            int leaveVertiportTimeWindow = (int) Math.floor((tripItemForOptimization.departureTime + tripItemForOptimization.accessTime + UAM_PROCESS_TIME + tripItemForOptimization.flightTime) / 3600);
            vertiportsCandidatesMap.get(tripItemForOptimization.accessVertiport.ID).tempConcurrentSaturationRates.put(arriveVertiportTimeWindow, vertiportsCandidatesMap.get(tripItemForOptimization.accessVertiport.ID).tempConcurrentSaturationRates.get(arriveVertiportTimeWindow) + tripItemForOptimization.tempUAMProbability / tripItemForOptimization.accessVertiport.capacity);
            vertiportsCandidatesMap.get(tripItemForOptimization.egressVertiport.ID).tempConcurrentSaturationRates.put(leaveVertiportTimeWindow, vertiportsCandidatesMap.get(tripItemForOptimization.egressVertiport.ID).tempConcurrentSaturationRates.get(leaveVertiportTimeWindow) + tripItemForOptimization.tempUAMProbability / tripItemForOptimization.egressVertiport.capacity);




            generalizedTravelCostAfter= modeSamples[1]*tripItemForOptimization.carGeneralizedCost+modeSamples[2]*tripItemForOptimization.ptGeneralizedCost+modeSamples[0]* tripItemForOptimization.UAMGeneralizedCost;
            emissionAfter=modeSamples[1]*tripItemForOptimization.carEmission+modeSamples[2]*tripItemForOptimization.ptEmission+modeSamples[0]* tripItemForOptimization.uamEmission;
            travelTimeAfter=modeSamples[1]*tripItemForOptimization.carTravelTime+modeSamples[2]*tripItemForOptimization.ptTravelTime+modeSamples[0]* tripItemForOptimization.uamTravelTime;
            double savedGCOneTrip=Double.max(tripItemForOptimization.currentGeneralizedCost-generalizedTravelCostAfter,0);
            double savedEmissionOneTrip=Double.max(tripItemForOptimization.currentEmission-emissionAfter,0);
            double savedTravelTimeOneTrip=Double.max(tripItemForOptimization.currentTravelTime-travelTimeAfter,0);
            if (tripItemForOptimization.tripPurpose.startsWith("H") && considerReturnTrip){
                savedGCOneTrip=savedGCOneTrip*2;
                savedEmissionOneTrip=savedEmissionOneTrip*2;
                savedTravelTimeOneTrip=savedTravelTimeOneTrip*2;
            }

            tripItemForOptimization.tempSavedGeneralizedCosts.add(savedGCOneTrip);
            tripItemForOptimization.tempSavedEmission.add(savedEmissionOneTrip);
            tripItemForOptimization.tempSavedTravelTime.add(savedTravelTimeOneTrip);
    }

    public static void writeTripBasedAnalysisResult(String outputFileName, List<TripItemForOptimization> tripItems) throws IOException, IOException {
// Build up a csv writer with no quote character
        FileWriter fileWriter = new FileWriter(outputFileName);

        // CSVWriter constructor with separator, quote character, and escape character
        CSVWriter writer = new CSVWriter(fileWriter, CSVWriter.DEFAULT_SEPARATOR, CSVWriter.NO_QUOTE_CHARACTER, CSVWriter.DEFAULT_ESCAPE_CHARACTER, CSVWriter.DEFAULT_LINE_END);

        String[] header = {"TripID", "HH_income", "Age", "TripPurpose", "UAM probability", "Car probability","PT probability", "car_distance","pt_distance", "SavedGeneralizedCost","SavedEmission","SavedTravelTime","Access Vertiport ID","Access Travel Time","Access Distance","Access Mode", "Egress Vertiport ID","Egress Travel Time","Egress Distance","Egress Mode","Flight Distance","Flight Time"};
        writer.writeNext(header);
        for (TripItemForOptimization tripItem : tripItems) {
            if (tripItem.isUAMAvailable) {
            double accessModeIndex = tripItem.originNeighborVertiportsTimeAndDistance.get(tripItem.accessVertiport).get("accessMode") ;
            double egressModeIndex = tripItem.destinationNeighborVertiportsTimeAndDistance.get(tripItem.egressVertiport).get("egressMode") ;
            String accessMode = accessModeIndex == 0 ? "Walk" : accessModeIndex == 1 ? "Car" : "PT";
            String egressMode = egressModeIndex == 0 ? "Walk" : egressModeIndex == 1 ? "Car" : "PT";
            double accessDistance = tripItem.originNeighborVertiportsTimeAndDistance.get(tripItem.accessVertiport).get("distance");
            double egressDistance = tripItem.destinationNeighborVertiportsTimeAndDistance.get(tripItem.egressVertiport).get("distance");
            double accessTime = tripItem.originNeighborVertiportsTimeAndDistance.get(tripItem.accessVertiport).get("travelTime");
            double egressTime = tripItem.destinationNeighborVertiportsTimeAndDistance.get(tripItem.egressVertiport).get("travelTime");
            double flightDistance= calculateEuciDistance(tripItem.accessVertiport.coord,tripItem.egressVertiport.coord);
            double flightTime=flightDistance/flightSpeed+takeOffLandingTime;
            int accessVertiportID = tripItem.accessVertiport.ID;
            int egressVertiportID = tripItem.egressVertiport.ID;
            String[] data = {tripItem.tripID, String.valueOf(tripItem.HH_income),String.valueOf(tripItem.age),String.valueOf(tripItem.tripPurpose),String.valueOf(tripItem.uamProbability), String.valueOf(tripItem.carProbability),String.valueOf(tripItem.ptProbability),String.valueOf(tripItem.carTripLength), String.valueOf(tripItem.ptTripLength),String.valueOf(tripItem.savedGeneralizedCost),String.valueOf(tripItem.savedEmission),String.valueOf(tripItem.savedTravelTime),String.valueOf(accessVertiportID),String.valueOf(accessTime),String.valueOf(accessDistance),accessMode,String.valueOf(egressVertiportID),String.valueOf(egressTime),String.valueOf(egressDistance),egressMode,String.valueOf(flightDistance),String.valueOf(flightTime)};
            writer.writeNext(data);
        }
        else{
            String[] data = {tripItem.tripID, String.valueOf(tripItem.HH_income),String.valueOf(tripItem.age),String.valueOf(tripItem.tripPurpose),String.valueOf(0),String.valueOf(tripItem.carProbability),String.valueOf(tripItem.ptProbability),String.valueOf(tripItem.carTripLength), String.valueOf(tripItem.ptTripLength),"NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA"};
            writer.writeNext(data);
            }
        }
        writer.close();
    }
    public static List<List<Integer>> generateNewSolution(Random random, List<Integer> currentSolutionID, HashMap<Integer,Vertiport> vertiportsCandidatesMap) {
        List<Integer> newSolutionID = new ArrayList<>(currentSolutionID);
        List<Integer> differenceID = new ArrayList<>();
        List<Integer> notChosenVertiportID = new ArrayList<>();
        List<Integer> saturationVertiportsID = new ArrayList<>();
        List<Integer> notSaturationVertiportsID = new ArrayList<>();
        List<List<Integer>> result = new ArrayList<>();


        for (Integer vertiportID : currentSolutionID) {
            Vertiport vertiport = vertiportsCandidatesMap.get(vertiportID);
            // calculate the max saturation rate of each vertiport
            if (vertiport.maxSaturationRate > 1) {
                saturationVertiportsID.add(vertiport.ID);
            } else {
                notSaturationVertiportsID.add(vertiport.ID); // This should be outside the loop
            }
        }

        for (Vertiport vertiport : vertiportsCandidatesMap.values()) {
            if (!currentSolutionID.contains(vertiport.ID)) {
                notChosenVertiportID.add(vertiport.ID);
            }
        }

        while (!saturationVertiportsID.isEmpty()) {
            Integer vertiportWithHighestSaturationRateID = selectHighestSaturationVertiport(saturationVertiportsID, (List<Vertiport>) vertiportsCandidatesMap.values().stream());
            List<Integer> neighborsID = getNeighborsID(vertiportsCandidatesMap.get(vertiportWithHighestSaturationRateID));

            if (!neighborsID.isEmpty() && !currentSolutionID.containsAll(neighborsID) ) { // Check if the neighbors are in the current solution or there is no neighbor
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
    public static List<List<Integer>> generateNewSolution(Random random, List<Integer> currentSolutionID, List<Vertiport> vertiportsCandidates,List<Vertiport> alreadySelectedVertiports) {
        HashMap<Integer,Vertiport> vertiportsCandidatesMap = new HashMap<>();
        for (Vertiport vertiport : vertiportsCandidates) {
            vertiportsCandidatesMap.put(vertiport.ID, vertiport);
        }
        List<Integer> newSolutionID = new ArrayList<>(currentSolutionID);
        List<Integer> differenceID = new ArrayList<>();
        List<Integer> alreadySelectedVertiportID = new ArrayList<>();
        List<Integer> notChosenVertiportID = new ArrayList<>();
        List<Integer> saturationVertiportsID = new ArrayList<>();
        List<Integer> notSaturationVertiportsID = new ArrayList<>();
        List<List<Integer>> result = new ArrayList<>();
        List<Integer> currentSelectedVertiportID= new ArrayList<>();
        for (Vertiport vertiport : alreadySelectedVertiports) {
            alreadySelectedVertiportID.add(vertiport.ID);
        }
        // Combine the already selected vertiports with the current solution
        currentSelectedVertiportID.addAll(currentSolutionID);
        currentSelectedVertiportID.addAll(alreadySelectedVertiportID);

        for (Integer vertiportID : currentSolutionID) {
            Vertiport vertiport = vertiportsCandidatesMap.get(vertiportID);
            // calculate the max saturation rate of each vertiport
            if (vertiport.maxSaturationRate > 1) {
                saturationVertiportsID.add(vertiport.ID);
            } else {
                notSaturationVertiportsID.add(vertiport.ID); // This should be outside the loop
            }
        }
        for (Vertiport vertiport: alreadySelectedVertiports) {
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

            if (!neighborsID.isEmpty() && !new HashSet<>(currentSelectedVertiportID).containsAll(neighborsID)  ) { // Check if the neighbors are in the current solution or there is no neighbor
                Integer newVertiportID = randomSelectExcluding(neighborsID, currentSelectedVertiportID, random);
                Integer removedVertiportID = notSaturationVertiportsID.get(random.nextInt(notSaturationVertiportsID.size()));
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
        return Math.sqrt(Math.pow(coord1.getX() - coord2.getX(), 2) + Math.pow(coord1.getY() - coord2.getY(), 2));
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
                if (vertiport.ID<0) {
                    duplicates.add(vertiport);
                }

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
        HashMap<Integer, List<Vertiport>> clusters = new HashMap<>();

        // Group the vertiports by cluster
        for (int i = 0; i < assignments.length; i++) {
            int clusterId = assignments[i]+1;
            Vertiport vertiport = vertiports.get(i);
            clusters.computeIfAbsent(clusterId, k -> new ArrayList<>()).add(vertiport);
        }

        // Compute the capacity and average construction cost for each cluster
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
            for (int i = 0; i < SIMULATION_HOURS * 4; i++) {
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

    public static List<Vertiport> summarizeClustersForAlreadySelectedVertiports(List<Vertiport> vertiports, int[] assignments) {
        HashMap<Integer, List<Vertiport>> clusters = new HashMap<>();

        // Group the vertiports by cluster
        for (int i = 0; i < assignments.length; i++) {
            int clusterId = (assignments[i]+1)*(-1);
            Vertiport vertiport = vertiports.get(i);
            clusters.computeIfAbsent(clusterId, k -> new ArrayList<>()).add(vertiport);
        }

        // Compute the capacity and average construction cost for each cluster
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
            summarizedVertiport.ID = (clusterId);
            summarizedVertiport.capacity = totalCapacity;
            summarizedVertiport.constructionCost = averageCost;
            summarizedVertiport.coord = new Coord(centroidX, centroidY);
            summarizedVertiport.concurrentSaturationRates = new ConcurrentHashMap<>();
            summarizedVertiport.tempConcurrentSaturationRates = new ConcurrentHashMap<>();
            for (int i = 0; i < SIMULATION_HOURS * 4; i++) {
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
            int clusterId = assignments[i]+1;
            List<Vertiport> cluster = clusters.computeIfAbsent(clusterId, k -> new ArrayList<>());
            cluster.add(vertiports.get(i));
        }
        return clusters;
    }
    public static Map<Integer, List<Vertiport>> saveClusterResultsForExistingVertiports(List<Vertiport> vertiports, int[] assignments) {
        Map<Integer, List<Vertiport>> clusters = new HashMap<>();
        for (int i = 0; i < assignments.length; i++) {
            int clusterId = (assignments[i]+1)*(-1);
            List<Vertiport> cluster = clusters.computeIfAbsent(clusterId, k -> new ArrayList<>());
            cluster.add(vertiports.get(i));
        }
        return clusters;
    }
    public static void findVertiportsNeighbours (List<Vertiport> vertiports, double NEIGHBOR_DISTANCE) {
        // clear the neighbours for all vertiports
        for (Vertiport vertiport : vertiports) {
            vertiport.neighbors.clear();
        }
        for (Vertiport vertiport : vertiports) {
            for (Vertiport neighbour : vertiports) {
                if (vertiport.ID != neighbour.ID) {
                    double distance = calculateEuciDistance(vertiport.coord, neighbour.coord);
                    if (distance <= NEIGHBOR_DISTANCE) {
                        vertiport.neighbors.add(neighbour);
                    }
                }
            }
        }
    }
}
