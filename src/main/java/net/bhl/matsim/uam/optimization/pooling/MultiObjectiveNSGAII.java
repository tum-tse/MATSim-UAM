package net.bhl.matsim.uam.optimization.pooling;

import net.bhl.matsim.uam.optimization.SimulatedAnnealingForPartD;
import net.bhl.matsim.uam.optimization.Vertiport;
import net.bhl.matsim.uam.optimization.pooling.sn.OptimizationResult;
import net.bhl.matsim.uam.optimization.pooling.sn.UAMOptimizationController;
import net.bhl.matsim.uam.optimization.utils.TripItemForOptimization;
import org.apache.log4j.Level;
import org.matsim.api.core.v01.Coord;
import org.matsim.api.core.v01.Id;
import org.matsim.api.core.v01.network.Link;
import org.matsim.api.core.v01.network.Network;
import org.matsim.contrib.dvrp.fleet.DvrpVehicle;
import org.matsim.contrib.dvrp.fleet.ImmutableDvrpVehicleSpecification;
import net.bhl.matsim.uam.infrastructure.UAMStation;
import net.bhl.matsim.uam.infrastructure.UAMVehicle;
import net.bhl.matsim.uam.infrastructure.UAMVehicleType;

import org.apache.log4j.Logger;
import org.matsim.core.network.NetworkUtils;
import org.matsim.core.network.io.MatsimNetworkReader;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.io.FileWriter;

public class MultiObjectiveNSGAII {
    private static final Logger log = Logger.getLogger(MultiObjectiveNSGAII.class);

    // Genetic Algorithm parameters ====================================================================================
    private static final int MAX_GENERATIONS = 10000; // Max number of generations
    //private static final int CROSSOVER_DISABLE_AFTER = 100; // New field to control when to stop crossover
    private static final int POP_SIZE = 200; // Population size
    private static final double MUTATION_RATE = 0.05; // Mutation rate
    private static final double CROSSOVER_RATE = 0.7; // Crossover rate
    private static final int TOURNAMENT_SIZE = 5; // Tournament size for selection
    private boolean ENABLE_LOCAL_SEARCH = false; // Enable local search after each generation
    private boolean ENABLE_PRINT_RESULTS = true; // Enable printing results to the CSVs

    private static final double ALPHA = -2.02 * 0.9101 / 1000; // Weight for changed flight distances
    private static final double BETA = -64.0 / 3600; // Weight for change in travel time
    //private static final double BETA_CRUCIAL_TIME_ChANGE = - 0.1; //TODO: need to reconsider the value
    private static final double PENALTY_FOR_VEHICLE_CAPACITY_VIOLATION = -10000;
    private static final double REVERT_SIGN = -1.0;

    private final long SEED = 4711; // MATSim default Random Seed
    private final Random rand = new Random(SEED);

    // Parameters and constant for the UAM problem =====================================================================
    private int FIRST_UAM_VEHICLE_ID = 1;
    private static final int VALUE_FOR_NO_VEHICLE_AVAILABLE = -1; // For example, using -1 as an indicator of no vehicle available for a trip
    private static final double END_SERVICE_TIME_OF_THE_DAY = 3600*36; // End service time of the day
    public static final double VEHICLE_CRUISE_SPEED = 350000.0 / 3600.0; // Vehicle cruise speed in m/s
    public static final int VEHICLE_CAPACITY = 4; // Vehicle capacity

    // Variables for the UAM problem ===================================================================================
    private double BUFFER_START_TIME = 3600*7; // Buffer start time for the first trip
    private double BUFFER_END_TIME = 3600*7+600; // Buffer end time for the last trip
    private double SEARCH_RADIUS_ORIGIN = 10000; // search radius for origin station
    private double SEARCH_RADIUS_DESTINATION = 10000; // search radius for destination station

    // Data container for the UAM problem ==============================================================================
    //private static List<UAMTrip> trips;
    private List<TripItemForOptimization> subTrips;
    //private static Map<Id<UAMStation>, UAMStation> stations;
    private HashMap<Integer, Vertiport> vertiportsMap;
    private final Map<Id<Vertiport>, List<UAMVehicle>> originStationVehicleMap = new HashMap<>();
    private final Map<Id<DvrpVehicle>, Vertiport> vehicleOriginStationMap = new HashMap<>();
    private final Map<Id<DvrpVehicle>, Vertiport> vehicleDestinationStationMap = new HashMap<>();
    private Map<String, List<UAMVehicle>> tripVehicleMap = new HashMap<>(); // Update tripVehicleMap to use ConcurrentHashMap
    //private static final Map<UAMVehicle, Integer> vehicleOccupancyMap = new HashMap<>();

    // Data container for outputs
    private final PriorityQueue<SolutionFitnessPair> solutionsHeap = new PriorityQueue<>(Comparator.comparingDouble(p -> p.getFitness()[0])); // Modify comparator to use first fitness objective
    private final PriorityQueue<SolutionFitnessPair> repairedSolutionsHeap = new PriorityQueue<>(Comparator.comparingDouble(p -> p.getFitness()[0])); // Modify comparator to use first fitness objective
    private final PriorityQueue<SolutionFitnessPair> bestSolutionsAcrossGenerations = new PriorityQueue<>(
            (a, b) -> Double.compare(a.getFitness()[0], b.getFitness()[0])  // Min heap
    );
    private final int MAX_BEST_SOLUTIONS = 100;  // Adjust as needed

    private static double nonPooledDeadheadingDistance;
    private static int nonPooledFleetSize;

    // Parallel computing
    private static final int numProcessors = Runtime.getRuntime().availableProcessors();
    private static final int bufferDivider = 1;

    // io paths
    private static String outputFile = "src/main/java/net/bhl/matsim/uam/optimization/pooling/output/";
    private double POOLING_TIME_WINDOW = BUFFER_END_TIME - BUFFER_START_TIME;
    private String subFolder = POOLING_TIME_WINDOW + "_" + SEARCH_RADIUS_ORIGIN + "_" + SEARCH_RADIUS_DESTINATION + "/";
    private String outputSubFolder = !outputFile.endsWith("/") ? outputFile + "/" + subFolder : outputFile + subFolder;
    // TODO: Create an initial population of solutions using domain-specific knowledge (in our case is the vehicles which were used to create the initial fleet of the vehicles).

    private final double SHARED_RIDE_TRAVEL_TIME_CHANGE_THRESHOLD = POOLING_TIME_WINDOW; // Threshold for shared ride travel time change

    // For travel time calculator
    private static String tripItemFile;
    private static String configFile;
    private static String vertiportUnitsCandidateFile;
    private static String scenarioConfigurations;
    private static SimulatedAnnealingForPartD.DataInitializer dataInitializer;
    private static Network network;

    public static void setFilePaths(String tripItemFilePath, String configFilePath, String vertiportUnitsCandidateFilePath, String scenarioConfigurationsPath, String outputFilePath) {
        tripItemFile = tripItemFilePath;
        configFile = configFilePath;
        vertiportUnitsCandidateFile = vertiportUnitsCandidateFilePath;
        scenarioConfigurations = scenarioConfigurationsPath;
        outputFile = outputFilePath;
    }

    // Constants for the GA solver =====================================================================================
    private int STABILITY_THRESHOLD = 100;
    private static final int BASE_STABILITY_THRESHOLD = 30;
    private static final int MAX_STABILITY_THRESHOLD = 200;
    private static final double PARETO_CHANGE_THRESHOLD = 1.0 / POP_SIZE;

    private List<SolutionFitnessPair> previousParetoFront = new ArrayList<>();
    private int stableGenerations = 0;
    // Constants for the localSearch solver ============================================================================
    private static final double MAX_ITERATIONS_SHARE_WITHOUT_IMPROVEMENT = 0.5;

    private static final double INITIAL_LOCAL_SEARCH_PROBABILITY = 0.1;
    private static final double FINAL_LOCAL_SEARCH_PROBABILITY = 0.5;

    // Setting for the shareability network ============================================================================
    private static final double MAX_DETOUR_RATIO = 0.3;
    private static final int MAX_CONNECTION_TIME_MINUTES = 30;

    // Settings for demand uncertainty and mode choice model ===========================================================
    // Number of Monte Carlo simulations per choice evaluation
    private static final int NUM_SIMULATIONS = 1000;
    private static final double NON_SHARED_UAM_FARE = 3.0; // euros per km

/*    // Static initializer block
    static {
        try {
            initializeData();
        } catch (IOException | InterruptedException e) {
            throw new RuntimeException(e);
        }
    }*/

    // Static method to initialize data
    private static void initializeData() throws IOException, InterruptedException {
        // Load data
/*        {
            DataLoader dataLoader = new DataLoader();
            dataLoader.loadAllData();
            //vehicles = dataLoader.getVehicles();
            stations = dataLoader.getStations();
        }*/
/*        Network network = NetworkUtils.createNetwork();
        new MatsimNetworkReader(network).readFile(uamScenarioInputPath + "/uam_network.xml.gz");
        UAMXMLReader uamReader = new UAMXMLReader(network);
        uamReader.readFile(uamScenarioInputPath + "/uam_vehicles.xml.gz");
        stations = uamReader.getStations();

        //subTrips = extractSubTrips(dataLoader.getUamTrips());
        String filePath = "scenarios/1-percent/sao_paulo_population2trips.csv";
        trips = readTripsFromCsv(filePath);*/

        dataInitializer = SimulatedAnnealingForPartD.getDataInitializer(tripItemFile, configFile, vertiportUnitsCandidateFile, scenarioConfigurations);

        network = NetworkUtils.createNetwork();
        new MatsimNetworkReader(network).readFile("./uam_routed_network.xml.gz");
    }

    // Constructor
    public MultiObjectiveNSGAII() {
        // The constructor is now empty as initialization is done in the static block
    }

    // Main method for testing
    public static void main(String[] args) throws IOException, InterruptedException {
        initialization(args);
        callAlgorithm(args);
        log.setLevel(Level.ALL);
    }

    public static double[] callAlgorithm(String[] args) {
        MultiObjectiveNSGAII instance = new MultiObjectiveNSGAII();
        return instance.runAlgorithm(args);
    }

    // Main method to run the the specified algorithm ==================================================================
    public static void initialization(String[] args) throws IOException, InterruptedException {
        if (args.length < 5) {
            System.out.println("Necessary: <Trip_Item> <Config> <Vertiport_Unit_Candidate> <Scenario_Configuration> <Result_Output> Optional: <BUFFER_END_TIME> <SEARCH_RADIUS_ORIGIN> <SEARCH_RADIUS_DESTINATION> <ENABLE_LOCAL_SEARCH> <ENABLE_PRINT_RESULTS>");
            System.exit(1);
        }

        // Provide the file via arguments
        setFilePaths(args[0], args[1], args[2], args[3], args[4]);

        // Reinitialize data with new file paths
        initializeData();

        log.setLevel(Level.ERROR);
    }
    public double[] runAlgorithm(String[] args) {
        if (args.length > 5) {
            BUFFER_END_TIME = BUFFER_START_TIME + Double.parseDouble(args[5]) * 60;
            SEARCH_RADIUS_ORIGIN = Double.parseDouble(args[6]);
            SEARCH_RADIUS_DESTINATION = Double.parseDouble(args[7]);
            ENABLE_LOCAL_SEARCH = Boolean.parseBoolean(args[8]);
            ENABLE_PRINT_RESULTS = Boolean.parseBoolean(args[9]);
            outputSubFolder = outputFile.endsWith("/") ? outputFile + args[10] : outputFile + "/" + args[10];
        }
        // Check if the output folder exists, if not create it
        createFolder(outputSubFolder);

        log.info("Scheduler started...");
        List<TripItemForOptimization> trips = dataInitializer.tripItems;
        log.info("The number of UAM trips before filtering: " + trips.size());

        vertiportsMap = dataInitializer.clusteredVertiportCandidatesMap;
        log.info("The number of Vertiports: " + dataInitializer.clusteredVertiportCandidatesMap.size());
        // Randomly select a share of trips from the list of subTrips
        subTrips = trips.stream()
                .filter(trip -> trip.departureTime >= BUFFER_START_TIME && trip.departureTime < BUFFER_END_TIME) // Add the filter
                .filter(trip -> calculateEuclideanDistance(findEuclideanNearestStation(trip, vertiportsMap, true).coord, trip.origin) <= SEARCH_RADIUS_ORIGIN)
                .filter(trip -> calculateEuclideanDistance(findEuclideanNearestStation(trip, vertiportsMap, false).coord, trip.destination) <= SEARCH_RADIUS_DESTINATION)
                .collect(Collectors.toCollection(ArrayList::new));
        log.info("The number of UAM trips after filtering: " + subTrips.size());

        // Initialize the origin station and destination station for each trip
        for (TripItemForOptimization uamTrip : subTrips) {
            uamTrip.accessVertiport = (findNearestStation(uamTrip, vertiportsMap, true));
            uamTrip.egressVertiport = (findNearestStation(uamTrip, vertiportsMap, false));
        }
        saveStationVehicleNumber(subTrips);
        tripVehicleMap = findNearbyVehiclesToTrips(subTrips);

        // Calculate the total deadheading distance for the non-pooled scenario
        calculateNonPooledDeadheadingDistance();

        List<SolutionFitnessPair> population = initializePopulation();
        for (int gen = 0; gen < MAX_GENERATIONS; gen++) {
            STABILITY_THRESHOLD = Math.min(BASE_STABILITY_THRESHOLD + gen / 10, MAX_STABILITY_THRESHOLD);

            population = evolvePopulation(population, gen);
            // Find the best solution in the current population
            SolutionFitnessPair bestSolution = Collections.max(population, Comparator.comparingDouble(p -> p.getFitness()[0]));
            log.info("Generation " + gen + ": Best fitness = " + Arrays.toString(bestSolution.getFitness()));

            if (stableGenerations >= STABILITY_THRESHOLD) {
                log.info("Terminating: Pareto front stable for " + STABILITY_THRESHOLD + " generations");
                break;
            }
        }
        // Fix or improve the final population as best solutions
        population = localSearch(population, 0);
        // Find the best feasible solution at the end of GA execution without altering the original solutions heap //TODO: Find the best feasible solution from all generations
        SolutionFitnessPair bestSolutionFitnessPair = calculatePopulationIndicators(population);
        int[] bestFeasibleSolution = bestSolutionFitnessPair.getSolution();
        double [] bestFeasibleSolutionFitness = bestSolutionFitnessPair.getFitness();
        SolutionIndicatorData bestSolutionIndicators = bestSolutionFitnessPair.getIndicatorData();

        log.info("Best feasible solution: " + Arrays.toString(bestFeasibleSolution));
        log.info("The fitness of the best feasible solution: " + Arrays.toString(bestFeasibleSolutionFitness));

        // Calculate and print the performance indicators to the console
        //printPerformanceIndicators(bestFeasibleSolution, bestSolutionIndicators);
        if (ENABLE_PRINT_RESULTS) {
            // Print statistics to CSV
            printStatisticsToCsv(bestFeasibleSolution, bestSolutionIndicators, outputSubFolder + "trip_statistics.csv");
        }

        // Print the NUMBER_OF_TRIPS_LONGER_THAN
        //log.info("Threshold for trips longer than " + THRESHOLD_FOR_TRIPS_LONGER_THAN_STRING + ": " + NUMBER_OF_TRIPS_LONGER_TAHN);

        return new double[]{BUFFER_END_TIME, SEARCH_RADIUS_ORIGIN, SEARCH_RADIUS_DESTINATION, bestFeasibleSolutionFitness[0]};
    }

    public static void createFolder(String folderPath) {
        File outputFolder = new File(folderPath);
        if (!outputFolder.exists()) {
            if (outputFolder.mkdirs()) {
                log.info("Output directory created: " + outputFolder.getAbsolutePath());
            } else {
                System.out.println("Failed to create output directory: " + outputFolder.getAbsolutePath());
            }
        } else {
            log.info("Output directory already exists: " + outputFolder.getAbsolutePath());
        }
    }

    // Method to calculate the non-pooled deadheading distance
    private void calculateNonPooledDeadheadingDistance() {
        // Create a map of vehicle assignments where each trip is served individually
        Map<Integer, List<TripItemForOptimization>> individualAssignments = new HashMap<>();

        for (int i = 0; i < subTrips.size(); i++) {
            TripItemForOptimization trip = subTrips.get(i);
            // Assign each trip to nearest vehicle/station
            //Vertiport originStation = findNearestStation(trip, vertiportsMap, true);
            //Vertiport destinationStation = findNearestStation(trip, vertiportsMap, false);

            // Create a single-trip list
            List<TripItemForOptimization> tripList = new ArrayList<>();
            tripList.add(trip);

            // Create a new vehicle ID for each trip since we're not pooling
            int vehicleId = i + 1;  // Simple sequential vehicle IDs
            individualAssignments.put(vehicleId, tripList);
        }

        // Use the shareability network to calculate deadheading distance
        UAMOptimizationController optimizer = new UAMOptimizationController(
                individualAssignments,
                MAX_DETOUR_RATIO,                // maxDetourRatio
                1,                  // maxPassengersPerVehicle (1 for non-pooled)
                MAX_CONNECTION_TIME_MINUTES,                 // maxConnectionTimeMinutes
                VEHICLE_CRUISE_SPEED,
                null,
                null
        );

        OptimizationResult result = optimizer.optimize();
        nonPooledFleetSize = result.getFleetSize();
        nonPooledDeadheadingDistance = result.getTotalDeadheadingFlightDistance();
    }
    public double getNonPooledDeadheadingDistance() {
        return nonPooledDeadheadingDistance;
    }
    public int getNonPooledFleetSize() { return nonPooledFleetSize; }

    // GA solver with NSGA-II modifications==============================================================================
    private List<SolutionFitnessPair> evolvePopulation(List<SolutionFitnessPair> population, int currentGeneration) {
        // Apply local search to improve the population after NSGA-II operations, and before offspring generation
        if (ENABLE_LOCAL_SEARCH && shouldApplyLocalSearch(currentGeneration)) {
            population = localSearch(population, currentGeneration);
        }

        List<SolutionFitnessPair> offspring = new ArrayList<>();

        // Generate new offspring
        while (offspring.size() < POP_SIZE) {
            int[] parent1 = selectParent(population);
            int[] parent2 = selectParent(population);
            int[] child;

            if (rand.nextDouble() < CROSSOVER_RATE) {
                child = crossover(parent1, parent2);
            } else {
                child = rand.nextBoolean() ? parent1.clone() : parent2.clone(); // Skip crossover, use parent directly
            }

            child = mutate(child);
            SolutionFitnessPair solution = calculateFitness(child, null, false);
            offspring.add(solution);
        }

        // Combine the old population and the new population
        List<SolutionFitnessPair> combinedPop = new ArrayList<>(population);
        combinedPop.addAll(offspring);

        // Apply local search to improve the combined population
        //combinedPop = localSearch(combinedPop);

        // Perform non-dominated sorting
        List<List<SolutionFitnessPair>> fronts = nonDominatedSort(combinedPop);

        // Select the next generation
        List<SolutionFitnessPair> nextGeneration = new ArrayList<>();
        for (List<SolutionFitnessPair> front : fronts) {
            calculateCrowdingDistance(front);
            if (nextGeneration.size() + front.size() <= POP_SIZE) {
                nextGeneration.addAll(front);
            } else {
                front.sort(Comparator.comparingInt(SolutionFitnessPair::getRank)
                        .thenComparingDouble(p -> p.getCrowdingDistance()).reversed());
                nextGeneration.addAll(front.subList(0, POP_SIZE - nextGeneration.size()));
                break;
            }
        }

        // Check for Pareto front stability
        List<SolutionFitnessPair> currentParetoFront = getNonDominatedSolutions(nextGeneration);
        double paretoChange = calculateParetoFrontChange(previousParetoFront, currentParetoFront);

        if (paretoChange < PARETO_CHANGE_THRESHOLD) {
            stableGenerations++;
        } else {
            stableGenerations = 0;
        }

        previousParetoFront = currentParetoFront;

        return nextGeneration;
    }
    private List<SolutionFitnessPair> getNonDominatedSolutions(List<SolutionFitnessPair> population) {
        return population.stream()
                .filter(solution -> solution.getRank() == 0)
                .collect(Collectors.toList());
    }

    private double calculateParetoFrontChange(List<SolutionFitnessPair> previous, List<SolutionFitnessPair> current) {
        if (previous.isEmpty()) {
            return 1.0;  // Maximum change if there was no previous front
        }

        int dominatingCount = 0;

        for (SolutionFitnessPair currentSolution : current) {
            boolean dominatesAny = false;
            for (SolutionFitnessPair previousSolution : previous) {
                if (dominated(previousSolution, currentSolution)) {
                    dominatesAny = true;
                    break;
                }
            }
            if (dominatesAny) {
                dominatingCount++;
            }
        }

        // Calculate the improvement as the proportion of current solutions that dominate any previous solution
        return (double) dominatingCount / current.size();
    }

    // Initialize population with random assignments
    private List<SolutionFitnessPair> initializePopulation() {
        List<SolutionFitnessPair> population = new ArrayList<>();
        for (int i = 0; i < MultiObjectiveNSGAII.POP_SIZE; i++) {
            int[] individual = generateIndividual();
            SolutionFitnessPair solution = calculateFitness(individual, null, false);
            population.add(solution);
        }
        return population;
    }

    // Generate a random individual
    private int[] generateIndividual() {
        int[] individual = new int[subTrips.size()];
        if (individual.length == 0) {
            log.info("Run: Pooling Time Window " + POOLING_TIME_WINDOW + "Origin Search Radius " + SEARCH_RADIUS_ORIGIN + "Destination Search Radius " + SEARCH_RADIUS_DESTINATION + " generated an empty individual. subTrips size: " + subTrips.size());
        }
        //resetVehicleCapacities(tripVehicleMap); // Reset the vehicle capacity since capacity of vehicles will be updated during each individual generation
        //resetVehicleOccupancy(vehicleOccupancyMap);
        for (int i = 0; i < individual.length; i++) {
            assignAvailableVehicle(i, individual);
        }
        return individual;
    }

    // Selection - Tournament selection with rank and crowding distance
    private int[] selectParent(List<SolutionFitnessPair> population) {
        List<SolutionFitnessPair> tournament = new ArrayList<>();
        for (int i = 0; i < TOURNAMENT_SIZE; i++) {
            tournament.add(population.get(rand.nextInt(population.size())));
        }
        tournament.sort(Comparator.comparingInt(SolutionFitnessPair::getRank)
                .thenComparingDouble(p -> p.getCrowdingDistance()).reversed()); // Sort by rank and then by crowding distance
        return tournament.get(0).getSolution();
    }

    //TODO: consider to repair the solution if it is not feasible due to the violation of vehicle capacity constraint
    // Crossover - Single point crossover //TODO: Implement other types of crossover instead of single point
    private int[] crossover(int[] parent1, int[] parent2) {
        int[] child = new int[parent1.length];

        int crossoverPoint = rand.nextInt(parent1.length);
        for (int i = 0; i < crossoverPoint; i++) {
            child[i] = parent1[i];
        }
        for (int i = crossoverPoint; i < parent2.length; i++) {
            child[i] = parent2[i];
        }
        return child;
    }

    // Mutation - Randomly change vehicle assignment
    private int[] mutate(int[] individual) {
        //resetVehicleCapacities(tripVehicleMap); // Reset the vehicle capacity since capacity of vehicles will be updated during each individual generation
        //resetVehicleOccupancy(vehicleOccupancyMap);
        for (int i = 0; i < individual.length; i++) {
            if (rand.nextDouble() < MUTATION_RATE) {
                assignAvailableVehicle(i, individual);
            }
        }
        return individual;
    }

    // Method to assign an available vehicle to a trip
    private void assignAvailableVehicle(int i, int[] individual) {
        TripItemForOptimization trip = subTrips.get(i);
        List<UAMVehicle> vehicleList = tripVehicleMap.get(trip.tripID);

        /*        //add occupancy constraint
        if (!vehicleList.isEmpty()) {
            Iterator<UAMVehicle> iterator0 = vehicleList.iterator();
            while (iterator0.hasNext()) {
                UAMVehicle vehicle = iterator0.next();
                int occupancy = vehicleOccupancyMap.get(vehicle);
                if (occupancy <= 0) {
                    iterator0.remove();
                }
            }
        }*/

        // synchronize the block that checks for empty vehicle list and adds a new vehicle (i.e., vehicle has capacity) after checking the egress constraint
            // Add a new vehicle for the trip when there is no available vehicle
            if (vehicleList == null || vehicleList.isEmpty()) {
                if (vehicleList == null) {
                    vehicleList = new ArrayList<>();
                }
                UAMVehicle newVehicle = feedDataForVehicleCreation(trip, false);
                vehicleList.add(newVehicle);
                tripVehicleMap.put(trip.tripID, vehicleList);
                //vehicleOccupancyMap.put(vehicle, VEHICLE_CAPACITY);

                // Update tripVehicleMap for all relevant trips
                updateTripVehicleMapForNewVehicle(newVehicle);
        }

        if (!vehicleList.isEmpty()) {

            //add access constraint
            int vehicleIndex = rand.nextInt(vehicleList.size());
            UAMVehicle selectedVehicle = vehicleList.get(vehicleIndex);
            //Integer currentCapacity = vehicleOccupancyMap.get(selectedVehicle);

            individual[i] = Integer.parseInt(selectedVehicle.getId().toString());

            // Decrement capacity and explicitly update tripVehicleMap
            //vehicleOccupancyMap.put(selectedVehicle, vehicleOccupancyMap.get(selectedVehicle) - 1);

        } else {
            // Handle the case when there is no available vehicle
            throw new IllegalArgumentException("Need to handle the case when there is no available vehicle for the trip.");
        }
    }
    private void updateTripVehicleMapForNewVehicle(UAMVehicle newVehicle) {
        Vertiport originStation = vehicleOriginStationMap.get(newVehicle.getId());
        Vertiport destinationStation = vehicleDestinationStationMap.get(newVehicle.getId());

        for (TripItemForOptimization trip : subTrips) {
            if (trip.originNeighborVertiportCandidatesTimeAndDistance.containsKey(originStation)&&trip.destinationNeighborVertiportCandidatesTimeAndDistance.containsKey(destinationStation)) {
                if (calculateEuclideanDistance(trip.origin, originStation.coord) <= SEARCH_RADIUS_ORIGIN &&
                        calculateEuclideanDistance(trip.destination, destinationStation.coord) <= SEARCH_RADIUS_DESTINATION) {

                    List<UAMVehicle> tripVehicles = tripVehicleMap.getOrDefault(trip.tripID, new ArrayList<>());
                    tripVehicles.add(newVehicle);
                    tripVehicleMap.put(trip.tripID, tripVehicles);
                }
            }
        }
    }

    // Objective function ==============================================================================================
    // Calculate fitness for an individual
    private SolutionFitnessPair calculateFitness(int[] individual, SolutionIndicatorData indicatorData, boolean isFinalSolutions) {
        Map<Integer, List<TripItemForOptimization>> vehicleAssignments = new HashMap<>();
        Map<Integer, Integer> vehicleLoadCount = new HashMap<>();
        Map<String, Double> travelTimeChangeMap = new HashMap<>();

        // Create UAM Mode Choice Model instance with our random seed
        UAMModeChoiceModel choiceModel = new UAMModeChoiceModel(rand, NUM_SIMULATIONS);

        // Organize trips by assigned vehicle
        for (int i = 0; i < individual.length; i++) {
            int vehicleId = individual[i];
            /*            if (vehicleId == VALUE_FOR_NO_VEHICLE_AVAILABLE) {
                continue;
            }*/
            vehicleAssignments.computeIfAbsent(vehicleId, k -> new ArrayList<>()).add(subTrips.get(i));

            // Update vehicle load count
            vehicleLoadCount.put(vehicleId, vehicleLoadCount.getOrDefault(vehicleId, 0) + 1);
        }

        double[] fitness = getFitnessForAllAssignment(isFinalSolutions, vehicleAssignments, travelTimeChangeMap, indicatorData, choiceModel);

        // Calculate pooling rate and vehicle capacity rates
        if (isFinalSolutions) {
            // Store fitness in indicatorData
            indicatorData.setFitness(fitness);
            indicatorData.setVehicleAssignments(vehicleAssignments);

            // Get choice model results for all trips
            Map<String, Map<Integer, Boolean>> allChoices = choiceModel.getAllSimulationChoices();
            Map<String, Double> allProbabilities = choiceModel.getAllAcceptanceProbabilities();

            // Calculate averages across all Monte Carlo scenarios
            double totalPooledTrips = 0;
            Map<Integer, Double> averageCapacityCount = new HashMap<>();
            int sharedRidesExceedingThreshold = 0;
            double totalAcceptedTrips = 0;
            Map<Integer, Double> vehicleUsageCount = new HashMap<>();

            // For each Monte Carlo scenario
            int numScenarios = NUM_SIMULATIONS;
            for (int scenario = 0; scenario < numScenarios; scenario++) {
                int scenarioPooledTrips = 0;
                Map<Integer, Integer> scenarioCapacityCount = new HashMap<>();
                int scenarioVehicles = 0;

                // Process each vehicle's assignments
                for (Map.Entry<Integer, List<TripItemForOptimization>> entry : vehicleAssignments.entrySet()) {
                    int vehicleId = entry.getKey();
                    List<TripItemForOptimization> trips = entry.getValue();

                    // Count accepted trips for this vehicle in this scenario
                    int acceptedTripsForVehicle = 0;
                    boolean anyTripAccepted = false;

                    for (TripItemForOptimization trip : trips) {
                        Map<Integer, Boolean> tripChoices = allChoices.get(trip.tripID);
                        if (tripChoices != null && tripChoices.get(scenario) != null && tripChoices.get(scenario)) {
                            acceptedTripsForVehicle++;
                            anyTripAccepted = true;

                            // Check if this is a pooled trip
                            if (trips.size() > 1) {
                                scenarioPooledTrips++;

                                // Check travel time threshold
                                double travelTimeChange = indicatorData.getTravelTimeChanges().get(trip.tripID);
                                if (travelTimeChange > SHARED_RIDE_TRAVEL_TIME_CHANGE_THRESHOLD) {
                                    sharedRidesExceedingThreshold++;
                                }
                            }
                        }
                    }

                    // If no trips were accepted for pooling, count each trip individually
                    if (!anyTripAccepted) {
                        acceptedTripsForVehicle = trips.size();  // Each trip requires its own vehicle
                        scenarioVehicles+=trips.size();
                    }

                    // Update capacity counts if vehicle is used
                    if (acceptedTripsForVehicle > 0) {
                        scenarioCapacityCount.put(acceptedTripsForVehicle,
                                scenarioCapacityCount.getOrDefault(acceptedTripsForVehicle, 0) + 1);
                        scenarioVehicles++;
                    }
                }

                // Accumulate results for this scenario
                totalPooledTrips += scenarioPooledTrips;

                // Update average capacity counts
                for (Map.Entry<Integer, Integer> entry : scenarioCapacityCount.entrySet()) {
                    averageCapacityCount.merge(entry.getKey(),
                            (double)entry.getValue() / numScenarios, Double::sum);
                }

                // Count vehicles used in this scenario
                vehicleUsageCount.merge(scenarioVehicles, 1.0, Double::sum);
            }

            // Calculate final metrics
            double averagePoolingRate = totalPooledTrips / (numScenarios * subTrips.size());
            indicatorData.setPoolingRate(averagePoolingRate);

            // Set vehicle capacity rates
            int totalVehicles = vehicleAssignments.size();
            for (int capacity = 0; capacity <= VEHICLE_CAPACITY; capacity++) {
                double rate = averageCapacityCount.getOrDefault(capacity, 0.0) / totalVehicles;
                indicatorData.getVehicleCapacityRates().put(capacity, rate);
            }

            // Calculate threshold rates
            double totalSharedTrips = totalPooledTrips / numScenarios;
            indicatorData.setSharedRidesExceedingThresholdRate(
                    totalSharedTrips == 0 ? 0 : (double) sharedRidesExceedingThreshold / (numScenarios * totalSharedTrips));
            indicatorData.setTotalSharedRidesExceedingThresholdRate(
                    (double) sharedRidesExceedingThreshold / (numScenarios * subTrips.size()));

            // Calculate average number of vehicles used
            double averageVehiclesUsed = vehicleUsageCount.entrySet().stream()
                    .mapToDouble(e -> e.getKey() * (e.getValue() / numScenarios))
                    .sum();
            indicatorData.setNumberOfUAMVehiclesUsed((int)Math.round(averageVehiclesUsed));
        }

        // Store vehicleLoadCount and travelTimeChangeMap in the solution pair
        return new SolutionFitnessPair(individual, fitness, vehicleLoadCount, travelTimeChangeMap);
    }
    private double[] getFitnessForAllAssignment(boolean isFinalSolutions, Map<Integer, List<TripItemForOptimization>> vehicleAssignments, Map<String, Double> travelTimeChangeMap,
                                                SolutionIndicatorData indicatorData, UAMModeChoiceModel choiceModel) {
        double totalFitness = 0.0;
        double totalFlightDistanceChange = 0.0;
        double totalTimeChange = 0.0;
        double totalViolationPenalty = 0.0;

        // Calculate fitness per vehicle
        for (Map.Entry<Integer, List<TripItemForOptimization>> entry : vehicleAssignments.entrySet()) {
            List<TripItemForOptimization> trips = entry.getValue();
            Vertiport originStationOfVehicle = vehicleOriginStationMap.get(Id.create(entry.getKey().toString(), DvrpVehicle.class));
            Vertiport destinationStationOfVehicle = vehicleDestinationStationMap.get(Id.create(entry.getKey().toString(), DvrpVehicle.class));

            // safety check
            if (trips.isEmpty()) continue;
            if (trips.size() == 1){
                TripItemForOptimization trip = trips.get(0);
                double[] fitnessValue = getFitnessForNonPooledOrBaseTrip(trip, originStationOfVehicle, destinationStationOfVehicle,
                        new double[]{totalFitness, totalFlightDistanceChange, totalTimeChange}, isFinalSolutions, travelTimeChangeMap,
                        indicatorData, null,0);
                continue;
            }

            // Find the base trip (the trip with the latest arrival time at the departure UAM station)
            TripItemForOptimization baseTrip = trips.get(0);
            for (TripItemForOptimization trip : trips) {
                double accessTimeOfBaseTrip = baseTrip.originNeighborVertiportCandidatesTimeAndDistance.get(originStationOfVehicle).get("travelTime");
                double accessTimeOfPooledTrip = trip.originNeighborVertiportCandidatesTimeAndDistance.get(originStationOfVehicle).get("travelTime");
                if ((trip.departureTime + accessTimeOfPooledTrip) > (baseTrip.departureTime + accessTimeOfBaseTrip)) {
                    baseTrip = trip;
                }
            }

            double boardingTimeForAllTrips = baseTrip.departureTime + baseTrip.originNeighborVertiportCandidatesTimeAndDistance.get(originStationOfVehicle).get("travelTime");
            // Calculate fitness based on the proposed pooling option
            for (TripItemForOptimization trip : trips) {
                if(trip.tripID.equals(baseTrip.tripID)){
                    double[] fitnessValue = getFitnessForNonPooledOrBaseTrip(trip, originStationOfVehicle, destinationStationOfVehicle,
                            new double[]{totalFitness, totalFlightDistanceChange, totalTimeChange}, isFinalSolutions, travelTimeChangeMap,
                            indicatorData, choiceModel, trips.size()-1);
                    continue;
                }

                double tripTotalFitness = 0.0;
                double tripTimeChange = 0.0;
                double tripFlightDistanceChange = 0.0;
                double tripTotalTravelTime;
                double originalFlightDistance;

                // calculate change in flight distance
                double flightDistanceChange = getFlightDistanceChange(trip, originStationOfVehicle, destinationStationOfVehicle);
                //totalFitness += ALPHA * flightDistanceChange;
                //tripFlightDistanceChange += flightDistanceChange;
                // calculate saved flight distance
                originalFlightDistance = calculateFlightDistance(trip.accessVertiport, trip.egressVertiport);
                double savedFlightDistance = originalFlightDistance;
                tripTotalFitness += ALPHA * (-1) * savedFlightDistance;
                tripFlightDistanceChange -= savedFlightDistance;
                // calculate change in flight time due to the change in flight distance
                double flightTimeChange = flightDistanceChange / VEHICLE_CRUISE_SPEED;
                tripTotalFitness += BETA * flightTimeChange;
                tripTimeChange += flightTimeChange;
                // calculate additional travel time
                double originalAccessTimeForThePooledTrip = trip.originNeighborVertiportCandidatesTimeAndDistance.get(trip.accessVertiport).get("travelTime");
                double originalArrivalTimeForThePooledTrip = trip.departureTime + originalAccessTimeForThePooledTrip;
                double travelTimeChangeDueToAccessMatching = boardingTimeForAllTrips - originalArrivalTimeForThePooledTrip;
                tripTimeChange += travelTimeChangeDueToAccessMatching;
                if(travelTimeChangeDueToAccessMatching > 0) {
                    tripTotalFitness += BETA * travelTimeChangeDueToAccessMatching;
                } else {
                    tripTotalFitness += BETA * (- travelTimeChangeDueToAccessMatching); //TODO: reconsider for "negative additional travel time" cases
                }
                double additionalTravelTimeDueToEgressMatching = getTravelTimeChangeDueToEgressMatching(trip, destinationStationOfVehicle);
                tripTotalFitness += BETA * additionalTravelTimeDueToEgressMatching;
                tripTimeChange += additionalTravelTimeDueToEgressMatching;

                travelTimeChangeMap.put(trip.tripID, tripTimeChange);

                //demand uncertainty
                // Get acceptance probability using choice model
                double acceptanceProbability;
                {
                    double accessTimeOfPooledTrip = trip.originNeighborVertiportCandidatesTimeAndDistance.get(originStationOfVehicle).get("travelTime");
                    //non-shared UAM related variables
                    double nonSharedInVehicleTime = originalFlightDistance / VEHICLE_CRUISE_SPEED;
                    double nonSharedCost = NON_SHARED_UAM_FARE * originalFlightDistance / 1000.0; // Convert meters to kilometers
                    //shared UAM related variables
                    double sharedInVehicleTime = originalFlightDistance / VEHICLE_CRUISE_SPEED + flightTimeChange;
                    double sharedWaitingTime = travelTimeChangeDueToAccessMatching - (accessTimeOfPooledTrip - originalAccessTimeForThePooledTrip);
                    int coPassengers = trips.size() - 1;
                    double sharedCost = SharedUAMCostCalculator.calculateSharedUAMCost(nonSharedCost, coPassengers, accessTimeOfPooledTrip);
                    acceptanceProbability = choiceModel.calculateAcceptanceProbability(trip.tripID, sharedCost, sharedInVehicleTime, accessTimeOfPooledTrip, sharedWaitingTime, coPassengers,
                            nonSharedCost, nonSharedInVehicleTime, originalAccessTimeForThePooledTrip);
                }
                totalFitness += acceptanceProbability * tripTotalFitness;
                totalFlightDistanceChange += acceptanceProbability * tripFlightDistanceChange;
                totalTimeChange += acceptanceProbability * tripTimeChange;

                // Store the statistics of the final solution
                if(isFinalSolutions){
                    indicatorData.setFlightDistanceChanges(trip.tripID, acceptanceProbability * tripFlightDistanceChange);
                    indicatorData.setTravelTimeChanges(trip.tripID, acceptanceProbability * tripTimeChange);

                    // Calculate departureRedirectionRate
                    double departureRedirectionRate = 0.0;
                    try {
                        double originStationDistance = trip.originNeighborVertiportCandidatesTimeAndDistance.get(originStationOfVehicle).get("distance");
                        double accessVertiportDistance = trip.originNeighborVertiportCandidatesTimeAndDistance.get(trip.accessVertiport).get("distance");

                        if (accessVertiportDistance != 0) {
                            departureRedirectionRate = (originStationDistance - accessVertiportDistance) / accessVertiportDistance;
                        } else {
                            // Handle the case where accessVertiportDistance is 0
                            departureRedirectionRate = originStationDistance > 0 ? Double.POSITIVE_INFINITY : 0.0;
                        }
                    } catch (NullPointerException e) {
                        // Handle the case where one of the get() operations returns null
                        log.warn("Null value encountered when calculating departureRedirectionRate for trip " + trip.tripID);
                        //departureRedirectionRate = 0.0;
                    }
                    // Check if the result is finite
                    if (!Double.isFinite(departureRedirectionRate)) {
                        log.warn("Non-finite departureRedirectionRate calculated for trip " + trip.tripID + ": " + departureRedirectionRate);
                        departureRedirectionRate = 0.0;  // or some other default value
                    }
                    indicatorData.setDepartureRedirectionRate(trip.tripID, acceptanceProbability * departureRedirectionRate);

                    // Calculate arrivalRedirectionRate
                    double arrivalRedirectionRate = 0.0;
                    try {
                        double destinationStationDistance = trip.destinationNeighborVertiportCandidatesTimeAndDistance.get(destinationStationOfVehicle).get("distance");
                        double egressVertiportDistance = trip.destinationNeighborVertiportCandidatesTimeAndDistance.get(trip.egressVertiport).get("distance");

                        if (egressVertiportDistance != 0) {
                            arrivalRedirectionRate = (destinationStationDistance - egressVertiportDistance) / egressVertiportDistance;
                        } else {
                            // Handle the case where egressVertiportDistance is 0
                            arrivalRedirectionRate = destinationStationDistance > 0 ? Double.POSITIVE_INFINITY : 0.0;
                        }
                    } catch (NullPointerException e) {
                        // Handle the case where one of the get() operations returns null
                        log.warn("Null value encountered when calculating arrivalRedirectionRate for trip " + trip.tripID);
                        //arrivalRedirectionRate = 0.0;
                    }
                    // Check if the result is finite
                    if (!Double.isFinite(arrivalRedirectionRate)) {
                        log.warn("Non-finite arrivalRedirectionRate calculated for trip " + trip.tripID + ": " + arrivalRedirectionRate);
                        arrivalRedirectionRate = 0.0;  // or some other default value
                    }
                    indicatorData.setArrivalRedirectionRate(trip.tripID, acceptanceProbability * arrivalRedirectionRate);

                    // total travel time for the trip
                    //TODO: Should the accessTime = boardingTimeForAllTrips - trip.getDepartureTime()?
                    tripTotalTravelTime = trip.originNeighborVertiportCandidatesTimeAndDistance.get(originStationOfVehicle).get("travelTime")
                            + calculateFlightDistance(originStationOfVehicle, destinationStationOfVehicle) / VEHICLE_CRUISE_SPEED
                            + trip.destinationNeighborVertiportCandidatesTimeAndDistance.get(destinationStationOfVehicle).get("travelTime");
                    indicatorData.setTotalTravelTime(trip.tripID, acceptanceProbability * tripTotalTravelTime);

                    // assigned origin station
                    indicatorData.setAssignedAccessStation(trip.tripID, String.valueOf(originStationOfVehicle.ID));
                    // assigned destination station
                    indicatorData.setAssignedEgressStation(trip.tripID, String.valueOf(destinationStationOfVehicle.ID));
                }

            }

            // Calculate capacity violation penalties across Monte Carlo scenarios =====================================
            double assignmentViolationPenalty = 0.0;
            Map<Integer, List<Boolean>> scenarioChoices = new HashMap<>();
            // First, collect all choices for each trip in this vehicle's assignments
            for (TripItemForOptimization trip : trips) {
                Map<Integer, Boolean> tripChoices = choiceModel.getChoicesForTrip(trip.tripID);
                if (tripChoices != null) {
                    // For each scenario
                    for (Map.Entry<Integer, Boolean> tripChoicesEntry : tripChoices.entrySet()) {
                        int scenarioNum = tripChoicesEntry.getKey();
                        boolean accepted = tripChoicesEntry.getValue();
                        scenarioChoices.computeIfAbsent(scenarioNum, k -> new ArrayList<>())
                                .add(accepted);
                    }
                }
            }
            // Calculate violations for each scenario
            for (List<Boolean> scenarioAcceptances : scenarioChoices.values()) {
                int acceptedTrips = 0;
                for (Boolean accepted : scenarioAcceptances) {
                    if (accepted) {
                        acceptedTrips++;
                    }
                }

                // If the number of accepted trips exceeds vehicle capacity
                if (acceptedTrips > VEHICLE_CAPACITY) {
                    assignmentViolationPenalty += PENALTY_FOR_VEHICLE_CAPACITY_VIOLATION *
                            (acceptedTrips - VEHICLE_CAPACITY);
                }
            }
            // Average the penalty across all scenarios
            if (!scenarioChoices.isEmpty()) {
                assignmentViolationPenalty /= scenarioChoices.size();
            }
            totalViolationPenalty += assignmentViolationPenalty;

        }
        return new double[]{totalFitness, REVERT_SIGN*totalFlightDistanceChange, REVERT_SIGN*totalTimeChange, totalViolationPenalty};
    }
    private double[] getFitnessForNonPooledOrBaseTrip(TripItemForOptimization trip, Vertiport originStationOfVehicle, Vertiport destinationStationOfVehicle, double[] fitnessValues,
                                                    boolean isFinalSolutions, Map<String, Double> travelTimeChangeMap, SolutionIndicatorData indicatorData, UAMModeChoiceModel choiceModel,
                                                    int coPassengers) {
        double tripTotalFitness = 0.0;
        double tripTimeChange = 0.0;
        double tripFlightDistanceChange = 0.0;

        // calculate change in flight distance
        double flightDistanceChange = getFlightDistanceChange(trip, originStationOfVehicle, destinationStationOfVehicle);
        tripTotalFitness += ALPHA * flightDistanceChange;
        tripFlightDistanceChange += flightDistanceChange;
        // calculate change in flight time due to the change in flight distance
        double flightTimeChange = flightDistanceChange / VEHICLE_CRUISE_SPEED;
        tripTotalFitness += BETA * flightTimeChange;
        tripTimeChange += flightTimeChange;
        // calculate change in travel time due to access matching
        double travelTimeChangeDueToAccessMatching = trip.originNeighborVertiportCandidatesTimeAndDistance.get(originStationOfVehicle).get("travelTime")
                - trip.originNeighborVertiportCandidatesTimeAndDistance.get(trip.accessVertiport).get("travelTime");
        tripTimeChange += travelTimeChangeDueToAccessMatching;
        if(travelTimeChangeDueToAccessMatching > 0) {
            tripTotalFitness += BETA * travelTimeChangeDueToAccessMatching;
        } else {
            tripTotalFitness += BETA * (- travelTimeChangeDueToAccessMatching);
        }
        // calculate change in travel time due to egress matching
        double travelTimeChangeDueToEgressMatching = getTravelTimeChangeDueToEgressMatching(trip, destinationStationOfVehicle);
        tripTotalFitness += BETA * travelTimeChangeDueToEgressMatching;
        tripTimeChange += travelTimeChangeDueToEgressMatching;

        travelTimeChangeMap.put(trip.tripID, tripTimeChange);

        double acceptanceProbability = 1;
        if(choiceModel != null) {
            //demand uncertainty
            // Get acceptance probability using choice model

            {
                double accessTimeOfPooledTrip = trip.originNeighborVertiportCandidatesTimeAndDistance.get(originStationOfVehicle).get("travelTime");
                double originalFlightDistance = calculateFlightDistance(trip.accessVertiport, trip.egressVertiport);
                double originalAccessTimeForThePooledTrip = trip.originNeighborVertiportCandidatesTimeAndDistance.get(trip.accessVertiport).get("travelTime");
                //non-shared UAM related variables
                double nonSharedInVehicleTime = originalFlightDistance / VEHICLE_CRUISE_SPEED;
                double nonSharedCost = NON_SHARED_UAM_FARE * originalFlightDistance / 1000.0; // Convert meters to kilometers
                //shared UAM related variables
                double sharedInVehicleTime = originalFlightDistance / VEHICLE_CRUISE_SPEED + flightTimeChange;
                double sharedWaitingTime = travelTimeChangeDueToAccessMatching - (accessTimeOfPooledTrip - originalAccessTimeForThePooledTrip);
                double sharedCost = SharedUAMCostCalculator.calculateSharedUAMCost(nonSharedCost, coPassengers, accessTimeOfPooledTrip);
                acceptanceProbability = choiceModel.calculateAcceptanceProbability(trip.tripID, sharedCost, sharedInVehicleTime, accessTimeOfPooledTrip, sharedWaitingTime, coPassengers,
                        nonSharedCost, nonSharedInVehicleTime, originalAccessTimeForThePooledTrip);
            }
        }
        fitnessValues[0] += acceptanceProbability * tripTotalFitness;
        fitnessValues[1] += acceptanceProbability * tripFlightDistanceChange;
        fitnessValues[2] += acceptanceProbability * tripTimeChange;


        // Store the statistics of the final solution
        if(isFinalSolutions){
            // Calculate UAM vehicle kilometer (passenger kilometer for non-pooled or base trip)
            double uamVehicleMeter = acceptanceProbability * calculateFlightDistance(originStationOfVehicle, destinationStationOfVehicle) /*/ 1000.0*/; // Convert meters to kilometers
            indicatorData.setUamVehicleMeter(indicatorData.getUamVehicleMeter() + uamVehicleMeter);

            indicatorData.setFlightDistanceChanges(trip.tripID, acceptanceProbability * tripFlightDistanceChange);

            indicatorData.setTravelTimeChanges(trip.tripID, acceptanceProbability * tripTimeChange);

            // Calculate departureRedirectionRate
            double departureRedirectionRate = 0.0;
            try {
                double originStationDistance = trip.originNeighborVertiportCandidatesTimeAndDistance.get(originStationOfVehicle).get("distance");
                double accessVertiportDistance = trip.originNeighborVertiportCandidatesTimeAndDistance.get(trip.accessVertiport).get("distance");

                if (accessVertiportDistance != 0) {
                    departureRedirectionRate = (originStationDistance - accessVertiportDistance) / accessVertiportDistance;
                } else {
                    // Handle the case where accessVertiportDistance is 0
                    departureRedirectionRate = originStationDistance > 0 ? Double.POSITIVE_INFINITY : 0.0;
                }
            } catch (NullPointerException e) {
                // Handle the case where one of the get() operations returns null
                log.warn("Null value encountered when calculating departureRedirectionRate for trip " + trip.tripID);
                //departureRedirectionRate = 0.0;
            }
            // Check if the result is finite
            if (!Double.isFinite(departureRedirectionRate)) {
                log.warn("Non-finite departureRedirectionRate calculated for trip " + trip.tripID + ": " + departureRedirectionRate);
                departureRedirectionRate = 0.0;  // or some other default value
            }
            indicatorData.setDepartureRedirectionRate(trip.tripID, acceptanceProbability * departureRedirectionRate);

            // Calculate arrivalRedirectionRate
            double arrivalRedirectionRate = 0.0;
            try {
                double destinationStationDistance = trip.destinationNeighborVertiportCandidatesTimeAndDistance.get(destinationStationOfVehicle).get("distance");
                double egressVertiportDistance = trip.destinationNeighborVertiportCandidatesTimeAndDistance.get(trip.egressVertiport).get("distance");

                if (egressVertiportDistance != 0) {
                    arrivalRedirectionRate = (destinationStationDistance - egressVertiportDistance) / egressVertiportDistance;
                } else {
                    // Handle the case where egressVertiportDistance is 0
                    arrivalRedirectionRate = destinationStationDistance > 0 ? Double.POSITIVE_INFINITY : 0.0;
                }
            } catch (NullPointerException e) {
                // Handle the case where one of the get() operations returns null
                log.warn("Null value encountered when calculating arrivalRedirectionRate for trip " + trip.tripID);
                //arrivalRedirectionRate = 0.0;
            }
            // Check if the result is finite
            if (!Double.isFinite(arrivalRedirectionRate)) {
                log.warn("Non-finite arrivalRedirectionRate calculated for trip " + trip.tripID + ": " + arrivalRedirectionRate);
                arrivalRedirectionRate = 0.0;  // or some other default value
            }
            indicatorData.setArrivalRedirectionRate(trip.tripID, acceptanceProbability * arrivalRedirectionRate);

            // total travel time for the trip
            double totalTravelTime = trip.originNeighborVertiportCandidatesTimeAndDistance.get(originStationOfVehicle).get("travelTime")
                    + calculateFlightDistance(originStationOfVehicle, destinationStationOfVehicle) / VEHICLE_CRUISE_SPEED
                    + trip.destinationNeighborVertiportCandidatesTimeAndDistance.get(destinationStationOfVehicle).get("travelTime");
            indicatorData.setTotalTravelTime(trip.tripID, acceptanceProbability * totalTravelTime);

            // assigned origin station
            indicatorData.setAssignedAccessStation(trip.tripID, String.valueOf(originStationOfVehicle.ID));
            // assigned destination station
            indicatorData.setAssignedEgressStation(trip.tripID, String.valueOf(destinationStationOfVehicle.ID));
        }
        return fitnessValues;
    }
    private static double getFlightDistanceChange(TripItemForOptimization trip, Vertiport originStationOfVehicle, Vertiport destinationStationOfVehicle) {
        return calculateFlightDistance(originStationOfVehicle, destinationStationOfVehicle) - calculateFlightDistance(trip.accessVertiport, trip.egressVertiport);
    }
    public static double calculateFlightDistance(Vertiport originStation, Vertiport destStation) {
        return calculateEuclideanDistance(originStation.coord, destStation.coord);
    }
    public static double calculateEuclideanDistance(Coord coord1, Coord coord2) {
        return Math.sqrt(Math.pow(coord1.getX() - coord2.getX(), 2) + Math.pow(coord1.getY() - coord2.getY(), 2));
    }
    private static double getTravelTimeChangeDueToEgressMatching(TripItemForOptimization trip, Vertiport destinationStationOfVehicle) {
        return trip.destinationNeighborVertiportCandidatesTimeAndDistance.get(destinationStationOfVehicle).get("travelTime") - trip.destinationNeighborVertiportCandidatesTimeAndDistance.get(trip.egressVertiport).get("travelTime");
    }

    // Helper methods for GA and NSGA-II ======================================================================
    private static List<List<SolutionFitnessPair>> nonDominatedSort(List<SolutionFitnessPair> population) {
        List<List<SolutionFitnessPair>> fronts = new ArrayList<>();
        Map<SolutionFitnessPair, List<SolutionFitnessPair>> dominationMap = new HashMap<>();
        Map<SolutionFitnessPair, Integer> dominatedCount = new HashMap<>();

        for (SolutionFitnessPair p : population) {
            dominationMap.put(p, new ArrayList<>());
            dominatedCount.put(p, 0);
            for (SolutionFitnessPair q : population) {
                if (dominated(p, q)) {
                    dominationMap.get(p).add(q);
                } else if (dominated(q, p)) {
                    dominatedCount.put(p, dominatedCount.get(p) + 1);
                }
            }
            if (dominatedCount.get(p) == 0) {
                p.rank = 0;
                if (fronts.isEmpty()) {
                    fronts.add(new ArrayList<>());
                }
                fronts.get(0).add(p);
            }
        }

        int i = 0;
        while (i < fronts.size() && !fronts.get(i).isEmpty()) {
            List<SolutionFitnessPair> nextFront = new ArrayList<>();
            for (SolutionFitnessPair p : fronts.get(i)) {
                for (SolutionFitnessPair q : dominationMap.get(p)) {
                    dominatedCount.put(q, dominatedCount.get(q) - 1);
                    if (dominatedCount.get(q) == 0) {
                        q.rank = i + 1;
                        nextFront.add(q);
                    }
                }
            }
            if (!nextFront.isEmpty()) {
                fronts.add(nextFront);
            }
            i++;
        }
        return fronts;
    }

    private Map<String, List<UAMVehicle>> findNearbyVehiclesToTrips(List<TripItemForOptimization> subTrips) {
        Map<String, List<UAMVehicle>> tripVehicleMap = new HashMap<>();
        for (TripItemForOptimization trip : subTrips) {
            for (Vertiport station : vertiportsMap.values()) {
                if (!trip.originNeighborVertiportCandidatesTimeAndDistance.containsKey(station)){
                    continue;
                }
                if (calculateEuclideanDistance(trip.origin, station.coord) <= SEARCH_RADIUS_ORIGIN) {
                    List<UAMVehicle> vehicles = originStationVehicleMap.get(Id.create(station.ID, Vertiport.class));
                    if (vehicles == null){
                        continue;
                    }
                    List<UAMVehicle> existingVehicles = tripVehicleMap.getOrDefault(trip.tripID, new ArrayList<>());

                    //add egress constraint
                    vehicles = vehicles.stream()
                            .filter(vehicle -> {
                                Vertiport vertiport = vehicleDestinationStationMap.get(vehicle.getId());
                                return trip.destinationNeighborVertiportCandidatesTimeAndDistance.containsKey(vertiport) &&
                                        calculateEuclideanDistance(trip.destination, vertiport.coord) <= SEARCH_RADIUS_DESTINATION;
                            })
                            .collect(Collectors.toCollection(ArrayList::new));

                    existingVehicles.addAll(vehicles);

                    tripVehicleMap.put(trip.tripID, existingVehicles);
                }
            }
        }
        return tripVehicleMap;
    }

    private static boolean dominated(SolutionFitnessPair p, SolutionFitnessPair q) {
        boolean betterInAnyObjective = false;
        for (int i = 0; i < p.getFitness().length; i++) {
            if (p.getFitness()[i] < q.getFitness()[i]) {
                betterInAnyObjective = true;
            } else if (p.getFitness()[i] > q.getFitness()[i]) {
                return false;
            }
        }
        return betterInAnyObjective;
    }

    private static void calculateCrowdingDistance(List<SolutionFitnessPair> front) {
        int n = front.size();
        if (n == 0) return;

        for (SolutionFitnessPair p : front) {
            p.crowdingDistance = 0;
        }

        int m = front.get(0).getFitness().length;
        for (int i = 0; i < m; i++) {
            final int objIndex = i;
            front.sort(Comparator.comparingDouble(p -> p.getFitness()[objIndex]));
            front.get(0).crowdingDistance = Double.POSITIVE_INFINITY;
            front.get(n - 1).crowdingDistance = Double.POSITIVE_INFINITY;
            double minValue = front.get(0).getFitness()[objIndex];
            double maxValue = front.get(n - 1).getFitness()[objIndex];
            for (int j = 1; j < n - 1; j++) {
                front.get(j).crowdingDistance += (front.get(j + 1).getFitness()[objIndex] - front.get(j - 1).getFitness()[objIndex]) / (maxValue - minValue);
            }
        }
    }

    // Method to count the vehicles by capacity
    @Deprecated
    private static Map<Integer, Integer> countVehicleCapacities(int[] solution) {
        Map<Integer, Integer> vehicleLoadCount = new HashMap<>();
        for (int vehicleId : solution) {
            vehicleLoadCount.put(vehicleId, vehicleLoadCount.getOrDefault(vehicleId, 0) + 1);
        }

        Map<Integer, Integer> capacityCount = new HashMap<>();
        for (int load : vehicleLoadCount.values()) {
            capacityCount.put(load, capacityCount.getOrDefault(load, 0) + 1);
        }

        return capacityCount;
    }

    // Local search methods ============================================================================================
    //TODO: Have not considered demand uncertainty!
    // Ruin and recreate solution

    private boolean shouldApplyLocalSearch(int currentGeneration) {
        double probability = INITIAL_LOCAL_SEARCH_PROBABILITY +
                (FINAL_LOCAL_SEARCH_PROBABILITY - INITIAL_LOCAL_SEARCH_PROBABILITY) *
                        (currentGeneration / (double)MAX_GENERATIONS);
        return rand.nextDouble() < probability;
    }
    // Adjusted Ruin and Recreate Methods
    // Method to determine the number of trips to ruin based on the current generation
    private int determineRuinDegree(int currentGeneration, int maxGenerations) {
        // Start with higher ruin degree and gradually decrease
        double ruinFactor = (1.0 - ((double) currentGeneration / maxGenerations));
        return (int) Math.ceil(ruinFactor * subTrips.size() / 2);
    }
    /*private static int[] ruinSolution(int[] solution, int currentGeneration, int maxGenerations) {
        int[] ruinedSolution = Arrays.copyOf(solution, solution.length);

        int numTripsToRuin = determineRuinDegree(currentGeneration, maxGenerations);

        Set<Integer> selectedIndices = new HashSet<>();
        while (selectedIndices.size() < numTripsToRuin) {
            int tripIndex = rand.nextInt(ruinedSolution.length);
            if (!selectedIndices.contains(tripIndex)) {
                selectedIndices.add(tripIndex);
                ruinedSolution[tripIndex] = VALUE_FOR_NO_VEHICLE_AVAILABLE; // Mark the trip as unassigned
            }
        }

        return ruinedSolution;
    }
    private static List<SolutionFitnessPair> localSearch(List<SolutionFitnessPair> population, int currentGeneration) {
        int maxGenerations = 100;
        List<SolutionFitnessPair> improvedPopulation = new ArrayList<>();

        for (SolutionFitnessPair solutionPair : population) {
            int[] currentSolution = solutionPair.getSolution();
            double[] currentFitness = solutionPair.getFitness();

            int[] bestSolution = Arrays.copyOf(currentSolution, currentSolution.length);
            double[] bestFitness = Arrays.copyOf(currentFitness, currentFitness.length);

            for (int i = 0; i < maxGenerations; i++) { // Number of iterations for local search
                int[] ruinedSolution = ruinSolution(bestSolution, currentGeneration, maxGenerations);
                int[] recreatedSolution = recreateSolution(ruinedSolution);
                double[] recreatedFitness = calculateFitness(recreatedSolution, false);

                if (dominates(new SolutionFitnessPair(recreatedSolution, recreatedFitness), new SolutionFitnessPair(bestSolution, bestFitness))) {
                    bestSolution = recreatedSolution;
                    bestFitness = recreatedFitness;
                }
            }

            improvedPopulation.add(new SolutionFitnessPair(bestSolution, bestFitness));
        }

        return improvedPopulation;
    }*/
    // Targeted Ruin - Focus on trips that are likely to improve the solution
    private int[] targetedRuin(SolutionFitnessPair solutionPair, int currentGeneration, int maxGenerations) {
        int[] ruinedSolution = Arrays.copyOf(solutionPair.getSolution(), solutionPair.getSolution().length);

        // Reset the vehicle load count and travel time change map
        Map<Integer, Integer> vehicleLoadCount = solutionPair.getVehicleLoadCount();
        Map<String, Double> travelTimeChangeMap = solutionPair.getTravelTimeChangeMap();

        List<Integer> targetTrips = new ArrayList<>();

        for (int i = 0; i < ruinedSolution.length; i++) {
            int vehicleId = ruinedSolution[i];

            boolean isOverCapacity = /*vehicleLoadCount.containsKey(vehicleId) &&*/ vehicleLoadCount.get(vehicleId) > VEHICLE_CAPACITY;
            boolean isSignificantTimeChange = /*vehicleLoadCount.containsKey(vehicleId) &&*/ vehicleLoadCount.get(vehicleId) > 1 &&
                    /*travelTimeChangeMap.containsKey(subTrips.get(i).tripID) &&*/
                    travelTimeChangeMap.get(subTrips.get(i).tripID) > SHARED_RIDE_TRAVEL_TIME_CHANGE_THRESHOLD;

            if (isOverCapacity || isSignificantTimeChange) {
                targetTrips.add(i);
            }
        }

        // Randomly select trips to ruin from the targeted trips
        int numTripsToRuin = Math.min(targetTrips.size(), determineRuinDegree(currentGeneration, maxGenerations));  // Adjust as needed
        Collections.shuffle(targetTrips);

        for (int i = 0; i < numTripsToRuin; i++) {
            int tripIndex = targetTrips.get(i);
            ruinedSolution[tripIndex] = VALUE_FOR_NO_VEHICLE_AVAILABLE; // Mark the trip as unassigned
        }

        return ruinedSolution;
    }
    private int[] recreateSolution(int[] ruinedSolution) {
        int[] recreatedSolution = Arrays.copyOf(ruinedSolution, ruinedSolution.length);

        for (int i = 0; i < recreatedSolution.length; i++) {
            if (recreatedSolution[i] == VALUE_FOR_NO_VEHICLE_AVAILABLE) {
                assignAvailableVehicle(i, recreatedSolution);
            }
        }

        return recreatedSolution;
    }
    private List<SolutionFitnessPair> localSearch(List<SolutionFitnessPair> population, int currentGeneration) {
        List<SolutionFitnessPair> improvedPopulation = new ArrayList<>();
        long startTime = System.currentTimeMillis();
        int maxIterations = 100; // Adjust as needed
        long maxRuntime = 1000; // 1 second, adjust as needed

        for (SolutionFitnessPair solutionPair : population) {
            SolutionFitnessPair bestSolution = solutionPair;
            int iterationsWithoutImprovement = 0;
            int iteration = 0;

            while (iterationsWithoutImprovement < MAX_ITERATIONS_SHARE_WITHOUT_IMPROVEMENT * maxIterations
                    && iteration < maxIterations
                    && (System.currentTimeMillis() - startTime) < maxRuntime) {

                int[] ruinedSolution = targetedRuin(bestSolution, currentGeneration, MAX_GENERATIONS);
                int[] recreatedSolution = recreateSolution(ruinedSolution);
                SolutionFitnessPair newSolution = calculateFitness(recreatedSolution, null, false);

                if (!dominated(newSolution, bestSolution) && !dominated(bestSolution, newSolution)) {
                    // If solutions are non-dominated, consider it an improvement
                    bestSolution = newSolution;
                    iterationsWithoutImprovement = 0;
                } else if (dominated(bestSolution, newSolution)) {
                    bestSolution = newSolution;
                    iterationsWithoutImprovement = 0;
                } else {
                    iterationsWithoutImprovement++;
                }

                iteration++;
            }

            improvedPopulation.add(bestSolution);
        }

        return improvedPopulation;
    }

    // SolutionFitnessPair related methods =============================================================================
    // SolutionFitnessPair class to hold individual solutions and their fitness, along with vehicleLoadCount and travelTimeChangeMap
    private static class SolutionFitnessPair {
        private final int[] solution;
        private final double[] fitness;
        private final Map<Integer, Integer> vehicleLoadCount;
        private final Map<String, Double> travelTimeChangeMap;
        private int rank;
        private double crowdingDistance;
        private SolutionIndicatorData indicatorData; // For saving the best solution's data

        public SolutionFitnessPair(int[] solution, double[] fitness) {
            this.solution = solution;
            this.fitness = fitness;
            this.vehicleLoadCount = null;
            this.travelTimeChangeMap = null;
            this.rank = Integer.MAX_VALUE;
            this.crowdingDistance = 0;
            this.indicatorData = null;
        }

        public SolutionFitnessPair(int[] solution, double[] fitness, Map<Integer, Integer> vehicleLoadCount,
                                   Map<String, Double> travelTimeChangeMap) {
            this.solution = solution;
            this.fitness = fitness;
            this.vehicleLoadCount = vehicleLoadCount;
            this.travelTimeChangeMap = travelTimeChangeMap;
            this.rank = Integer.MAX_VALUE;
            this.crowdingDistance = 0;
            this.indicatorData = null;
        }

        public void setIndicatorData(SolutionIndicatorData indicatorData) {
            this.indicatorData = indicatorData;
        }
        public SolutionIndicatorData getIndicatorData() {
            return indicatorData;
        }

        public int[] getSolution() {
            return solution;
        }

        public double[] getFitness() {
            return fitness;
        }

        public int getRank() {
            return rank;
        }

        public double getCrowdingDistance() {
            return crowdingDistance;
        }

        public Map<Integer, Integer> getVehicleLoadCount() {
            return vehicleLoadCount;
        }

        public Map<String, Double> getTravelTimeChangeMap() {
            return travelTimeChangeMap;
        }
    }

    // Method to find the first feasible solution from the priority queue without altering the original heap
    @Deprecated
    private SolutionFitnessPair findBestFeasibleSolution(List<SolutionFitnessPair> population) {
        // Create a new priority queue that is a copy of the original but sorted in descending order by fitness
        PriorityQueue<SolutionFitnessPair> solutionsHeapCopy = new PriorityQueue<>(
                Comparator.comparingDouble((SolutionFitnessPair p) -> p.getFitness()[0]).reversed()
        );
        solutionsHeapCopy.addAll(population);
        SolutionFitnessPair bestSolutionButMaybeInfeasible = solutionsHeapCopy.peek();

        // Iterate through the copied solutions heap to find a feasible solution
        while (!solutionsHeapCopy.isEmpty()) {
            SolutionFitnessPair solutionPair = solutionsHeapCopy.poll(); // Remove and retrieve the solution with the highest fitness
            int[] candidateSolution = solutionPair.getSolution();
            if (isFeasible(candidateSolution, false)) {
                return solutionPair;
            }
        }
        int[] quickFixedSolution = guaranteeFeasibleSolution(bestSolutionButMaybeInfeasible.getSolution());
        return new SolutionFitnessPair(quickFixedSolution, calculateFitness(quickFixedSolution, null, false).getFitness());
    }
    // Helper method to check if a solution violates vehicle capacity constraints
    private static boolean isFeasible(int[] solution, boolean isPrintCapacityViolation) {
        boolean isFeasible = true;
        int vehicleCapacityViolated = 0;
        Map<Integer, Integer> vehicleLoadCount = new HashMap<>();
        for (int vehicleId : solution) {
            vehicleLoadCount.put(vehicleId, vehicleLoadCount.getOrDefault(vehicleId, 0) + 1);
        }
        for (int load : vehicleLoadCount.values()) {
            if (load > VEHICLE_CAPACITY) {
                isFeasible = false;
                vehicleCapacityViolated++;
            }
        }
        if(isPrintCapacityViolation){
            log.info("Number of vehicles with capacity violation: " + vehicleCapacityViolated);
        }
        return isFeasible;
    }

    public int[] guaranteeFeasibleSolution(int[] solution) {
        boolean isSolutionFeasible = false;
        int iterationCount = 0;
        //final int MAX_ITERATIONS = 1000000; // A very high number, but not infinite

        while (!isSolutionFeasible /*&& iterationCount < MAX_ITERATIONS*/) {
            Map<Integer, List<Integer>> vehicleAssignments = new HashMap<>();

            // Group trips by assigned vehicle
            for (int i = 0; i < solution.length; i++) {
                int vehicleId = solution[i];
                vehicleAssignments.computeIfAbsent(vehicleId, k -> new ArrayList<>()).add(i);
            }

            // Identify and fix overloaded vehicles
            for (Map.Entry<Integer, List<Integer>> entry : vehicleAssignments.entrySet()) {
                int vehicleId = entry.getKey();
                List<Integer> assignedTrips = entry.getValue();

                if (assignedTrips.size() > VEHICLE_CAPACITY) {
                    // Sort trips by arrival time at the station
                    assignedTrips.sort((a, b) -> {
                        TripItemForOptimization tripA = subTrips.get(a);
                        TripItemForOptimization tripB = subTrips.get(b);
                        Vertiport station = vehicleOriginStationMap.get(Id.create(String.valueOf(vehicleId), DvrpVehicle.class));
                        double arrivalTimeA = tripA.departureTime + tripA.originNeighborVertiportCandidatesTimeAndDistance.get(station).get("travelTime");
                        double arrivalTimeB = tripB.departureTime + tripB.originNeighborVertiportCandidatesTimeAndDistance.get(station).get("travelTime");
                        return Double.compare(arrivalTimeA, arrivalTimeB);
                    });

                    // Reassign extra trips
                    for (int i = VEHICLE_CAPACITY; i < assignedTrips.size(); i++) {
                        int tripIndex = assignedTrips.get(i);
                        TripItemForOptimization trip = subTrips.get(tripIndex);

                        // Check if all available vehicles are at capacity
                        List<UAMVehicle> availableVehicles = new ArrayList<>(tripVehicleMap.get(trip.tripID));
                        boolean allVehiclesAtCapacity = availableVehicles.stream()
                                .allMatch(v -> vehicleAssignments.getOrDefault(Integer.parseInt(v.getId().toString()), Collections.emptyList()).size() >= VEHICLE_CAPACITY);

                        if (allVehiclesAtCapacity) {
                            // Create a new vehicle for this trip
                            UAMVehicle newVehicle = feedDataForVehicleCreation(trip, false);
                            availableVehicles.add(newVehicle);
                            tripVehicleMap.put(trip.tripID, availableVehicles);

                            // Update tripVehicleMap for the new vehicle
                            updateTripVehicleMapForNewVehicle(newVehicle);

                            // Assign the trip to the new vehicle
                            solution[tripIndex] = Integer.parseInt(newVehicle.getId().toString());
                        } else {
                            // Assign to an available vehicle that's not at capacity
                            assignAvailableVehicle(tripIndex, solution);
                        }
                    }
                }
            }

            // Check if the solution is now feasible
            isSolutionFeasible = isFeasible(solution, false);

            // If the solution is still infeasible
            // create new vehicles for all overloaded trips
            if (!isSolutionFeasible) {
                forceCreateNewVehicles(solution);
            }

            iterationCount++;
        }

/*        if (!isSolutionFeasible) {
            throw new IllegalArgumentException("Error: Could not find a feasible solution after " + MAX_ITERATIONS + " iterations.");
        }*/
        return solution;
    }
    private void forceCreateNewVehicles(int[] solution) {
        Map<Integer, List<Integer>> vehicleAssignments = new HashMap<>();

        // Group trips by assigned vehicle
        for (int i = 0; i < solution.length; i++) {
            int vehicleId = solution[i];
            vehicleAssignments.computeIfAbsent(vehicleId, k -> new ArrayList<>()).add(i);
        }

        for (Map.Entry<Integer, List<Integer>> entry : vehicleAssignments.entrySet()) {
            List<Integer> assignedTrips = entry.getValue();
            if (assignedTrips.size() > VEHICLE_CAPACITY) {
                for (int i = VEHICLE_CAPACITY; i < assignedTrips.size(); i++) {
                    int tripIndex = assignedTrips.get(i);
                    TripItemForOptimization trip = subTrips.get(tripIndex);
                    UAMVehicle newVehicle = feedDataForVehicleCreation(trip, false);
                    List<UAMVehicle> availableVehicles = new ArrayList<>(tripVehicleMap.getOrDefault(trip.tripID, new ArrayList<>()));
                    availableVehicles.add(newVehicle);
                    tripVehicleMap.put(trip.tripID, availableVehicles);

                    // Update tripVehicleMap for the new vehicle
                    updateTripVehicleMapForNewVehicle(newVehicle);

                    solution[tripIndex] = Integer.parseInt(newVehicle.getId().toString());
                }
            }
        }
    }

    // Performance indicators ==========================================================================================
    // Method to calculate and print the performance indicators
    private void printPerformanceIndicators(int[] solution, SolutionIndicatorData indicatorData) {
        // Print the pooling rate (already calculated)
        log.info("Pooling rate: " + indicatorData.getPoolingRate());

        // Print vehicle capacity rates (already calculated)
        log.info("Vehicle Capacity Rates:");
        for (int capacity = 0; capacity <= VEHICLE_CAPACITY; capacity++) {
            double rate = indicatorData.getVehicleCapacityRates().getOrDefault(capacity, 0.0);
            int count = (int)(rate * indicatorData.getNumberOfUAMVehiclesUsed());
            log.info("Capacity " + capacity + ": " + count + " vehicles, Rate: " + rate);
        }

        // Print travel time statistics (already calculated)
        log.info("Average travel time change: " + indicatorData.getAverageTravelTimeChange());
        log.info("5th percentile of travel time change: " + indicatorData.getPercentile5thTravelTimeChange());
        log.info("95th percentile of travel time change: " + indicatorData.getPercentile95thTravelTimeChange());

        // Print flight distance statistics (already calculated)
        log.info("Average flight distance change: " + indicatorData.getAverageFlightDistanceChange());
        log.info("5th percentile of flight distance change: " + indicatorData.getPercentile5thFlightDistanceChange());
        log.info("95th percentile of flight distance change: " + indicatorData.getPercentile95thFlightDistanceChange());

        // Print departure redirection rates (already calculated)
        log.info("Average departure redirection rate: " + indicatorData.getAverageDepartureRedirectionRate());
        log.info("5th percentile of departure redirection rate: " + indicatorData.getPercentile5thDepartureRedirectionRate());
        log.info("95th percentile of departure redirection rate: " + indicatorData.getPercentile95thDepartureRedirectionRate());

        // Print arrival redirection rates (already calculated)
        log.info("Average arrival redirection rate: " + indicatorData.getAverageArrivalRedirectionRate());
        log.info("5th percentile of arrival redirection rate: " + indicatorData.getPercentile5thArrivalRedirectionRate());
        log.info("95th percentile of arrival redirection rate: " + indicatorData.getPercentile95thArrivalRedirectionRate());

        // Print threshold statistics (already calculated)
        log.info("Share of shared rides with travel time changes exceeding " + SHARED_RIDE_TRAVEL_TIME_CHANGE_THRESHOLD +
                ": " + indicatorData.getSharedRidesExceedingThresholdRate());
        log.info("Total share of shared rides with travel time changes exceeding " + SHARED_RIDE_TRAVEL_TIME_CHANGE_THRESHOLD +
                ": " + indicatorData.getTotalSharedRidesExceedingThresholdRate());

        // Print UAM vehicle statistics
        log.info("Total UAM vehicle meter: " + indicatorData.getUamVehicleMeter());
        log.info("Number of UAM vehicles used: " + indicatorData.getNumberOfUAMVehiclesUsed());

        // Print deadheading and fleet size changes
        log.info("Deadheading flight distance change: " + (indicatorData.getDeadHeadingFlightDistance() - getNonPooledDeadheadingDistance()));
        log.info("Fleet size change: " + (indicatorData.getFleetSize() - getNonPooledFleetSize()));
    }
    // Method to print statistics to a CSV file
    private void printStatisticsToCsv(int[] solution, SolutionIndicatorData indicatorData, String fileName) {
        // Create a map to count the number of trips assigned to each vehicle
        Map<Integer, Integer> vehicleTripCount = new HashMap<>();
        for (int vehicleId : solution) {
            vehicleTripCount.put(vehicleId, vehicleTripCount.getOrDefault(vehicleId, 0) + 1);
        }

        try (FileWriter writer = new FileWriter(fileName)) {
            writer.append("TripId,AssignedVehicleId,AccessStationId,EgressStationId,TotalTravelTime,TravelTimeChange,FlightDistanceChange,DepartureRedirectionRate,ArrivalRedirectionRate,VehicleTripCount\n");
            for (int i = 0; i < solution.length; i++) {
                TripItemForOptimization trip = subTrips.get(i);
                String tripId = trip.tripID;
                int assignedVehicleId = solution[i];
                double travelTimeChange = indicatorData.getTravelTimeChanges().get(tripId);
                double flightDistanceChange = indicatorData.getFlightDistanceChanges().get(tripId);
                double departureRedirectionRate = indicatorData.getDepartureRedirectionRates().get(tripId);
                double arrivalRedirectionRate = indicatorData.getArrivalRedirectionRates().get(tripId);
                String assignedAccessStation = indicatorData.getAssignedAccessStations().get(tripId);
                String assignedEgressStation = indicatorData.getAssignedEgressStations().get(tripId);
                double totalTravelTime = indicatorData.getTotalTravelTimes().get(tripId);

                // Get the number of trips assigned to the same vehicle
                int tripCountForVehicle = vehicleTripCount.getOrDefault(assignedVehicleId, 0);

                writer.append(String.format("%s,%d,%s,%s,%.2f,%.2f,%.2f,%.2f,%.2f,%d\n",
                        tripId, assignedVehicleId, assignedAccessStation, assignedEgressStation, totalTravelTime, travelTimeChange, flightDistanceChange, departureRedirectionRate, arrivalRedirectionRate, tripCountForVehicle));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public class SolutionIndicatorData {
        private final int[] solution;
        private final Map<String, Double> travelTimeChanges = new HashMap<>();
        private final Map<String, Double> flightDistanceChanges = new HashMap<>();
        private final Map<String, Double> departureRedirectionRates = new HashMap<>();
        private final Map<String, Double> arrivalRedirectionRates = new HashMap<>();
        private final Map<String, Double> totalTravelTimes = new HashMap<>();
        private final Map<String, String> assignedAccessStations = new HashMap<>();
        private final Map<String, String> assignedEgressStations = new HashMap<>();

        private double[] fitness;
        private double poolingRate;
        private Map<Integer, Double> vehicleCapacityRates = new HashMap<>();
        private double sharedRidesExceedingThresholdRate;
        private double totalSharedRidesExceedingThresholdRate;
        private double averageTravelTimeChange;
        private double percentile5thTravelTimeChange;
        private double percentile95thTravelTimeChange;
        private double averageFlightDistanceChange;
        private double percentile5thFlightDistanceChange;
        private double percentile95thFlightDistanceChange;
        private double averageDepartureRedirectionRate;
        private double percentile5thDepartureRedirectionRate;
        private double percentile95thDepartureRedirectionRate;
        private double averageArrivalRedirectionRate;
        private double percentile5thArrivalRedirectionRate;
        private double percentile95thArrivalRedirectionRate;

        private double averageTotalTravelTime;
        private double percentile5thTotalTravelTime;
        private double percentile95thTotalTravelTime;

        private int numberOfUAMVehiclesUsed;
        private double uamVehicleMeter;

        private double deadHeadingFlightDistance;
        private int fleetSize;
        private Map<Integer, List<TripItemForOptimization>> vehicleAssignments;

        public SolutionIndicatorData(int[] solution) {
            this.solution = solution;
        }

        // Getters and setters for all fields
        public int[] getSolution() { return solution; }
        private Map<String, Double> getTravelTimeChanges() { return travelTimeChanges; }
        private Map<String, Double> getFlightDistanceChanges() { return flightDistanceChanges; }
        private Map<String, Double> getDepartureRedirectionRates() { return departureRedirectionRates; }
        private Map<String, Double> getArrivalRedirectionRates() { return arrivalRedirectionRates; }
        private Map<String, Double> getTotalTravelTimes() { return totalTravelTimes; }
        private  Map<String, String> getAssignedAccessStations() { return assignedAccessStations; }
        private  Map<String, String> getAssignedEgressStations() { return assignedEgressStations; }
        public double[] getFitness() { return fitness; }
        public void setFitness(double[] fitness) { this.fitness = fitness; }
        public double getPoolingRate() { return poolingRate; }
        public void setPoolingRate(double poolingRate) { this.poolingRate = poolingRate; }
        public Map<Integer, Double> getVehicleCapacityRates() { return vehicleCapacityRates; }
        //public void setVehicleCapacityRates(Map<Integer, Double> vehicleCapacityRates) { this.vehicleCapacityRates = vehicleCapacityRates; }
        public double getSharedRidesExceedingThresholdRate() { return sharedRidesExceedingThresholdRate; }
        public void setSharedRidesExceedingThresholdRate(double sharedRidesExceedingThresholdRate) { this.sharedRidesExceedingThresholdRate = sharedRidesExceedingThresholdRate; }
        public double getTotalSharedRidesExceedingThresholdRate() { return totalSharedRidesExceedingThresholdRate; }
        public void setTotalSharedRidesExceedingThresholdRate(double totalSharedRidesExceedingThresholdRate) { this.totalSharedRidesExceedingThresholdRate = totalSharedRidesExceedingThresholdRate; }
        public double getAverageTravelTimeChange() { return averageTravelTimeChange; }
        public void setAverageTravelTimeChange(double averageTravelTimeChange) { this.averageTravelTimeChange = averageTravelTimeChange; }
        public double getPercentile5thTravelTimeChange() { return percentile5thTravelTimeChange; }
        public void setPercentile5thTravelTimeChange(double percentile5thTravelTimeChange) { this.percentile5thTravelTimeChange = percentile5thTravelTimeChange; }
        public double getPercentile95thTravelTimeChange() { return percentile95thTravelTimeChange; }
        public void setPercentile95thTravelTimeChange(double percentile95thTravelTimeChange) { this.percentile95thTravelTimeChange = percentile95thTravelTimeChange; }
        public double getAverageFlightDistanceChange() { return averageFlightDistanceChange; }
        public void setAverageFlightDistanceChange(double averageFlightDistanceChange) { this.averageFlightDistanceChange = averageFlightDistanceChange; }
        public double getPercentile5thFlightDistanceChange() { return percentile5thFlightDistanceChange; }
        public void setPercentile5thFlightDistanceChange(double percentile5thFlightDistanceChange) { this.percentile5thFlightDistanceChange = percentile5thFlightDistanceChange; }
        public double getPercentile95thFlightDistanceChange() { return percentile95thFlightDistanceChange; }
        public void setPercentile95thFlightDistanceChange(double percentile95thFlightDistanceChange) { this.percentile95thFlightDistanceChange = percentile95thFlightDistanceChange; }
        public double getAverageDepartureRedirectionRate() { return averageDepartureRedirectionRate; }
        public void setAverageDepartureRedirectionRate(double averageDepartureRedirectionRate) { this.averageDepartureRedirectionRate = averageDepartureRedirectionRate; }
        public double getPercentile5thDepartureRedirectionRate() { return percentile5thDepartureRedirectionRate; }
        public void setPercentile5thDepartureRedirectionRate(double percentile5thDepartureRedirectionRate) { this.percentile5thDepartureRedirectionRate = percentile5thDepartureRedirectionRate; }
        public double getPercentile95thDepartureRedirectionRate() { return percentile95thDepartureRedirectionRate; }
        public void setPercentile95thDepartureRedirectionRate(double percentile95thDepartureRedirectionRate) { this.percentile95thDepartureRedirectionRate = percentile95thDepartureRedirectionRate; }
        public double getAverageArrivalRedirectionRate() { return averageArrivalRedirectionRate; }
        public void setAverageArrivalRedirectionRate(double averageArrivalRedirectionRate) { this.averageArrivalRedirectionRate = averageArrivalRedirectionRate; }
        public double getPercentile5thArrivalRedirectionRate() { return percentile5thArrivalRedirectionRate; }
        public void setPercentile5thArrivalRedirectionRate(double percentile5thArrivalRedirectionRate) { this.percentile5thArrivalRedirectionRate = percentile5thArrivalRedirectionRate; }
        public double getPercentile95thArrivalRedirectionRate() { return percentile95thArrivalRedirectionRate; }
        public void setPercentile95thArrivalRedirectionRate(double percentile95thArrivalRedirectionRate) { this.percentile95thArrivalRedirectionRate = percentile95thArrivalRedirectionRate; }

        public double getAverageTotalTravelTime() { return averageTotalTravelTime; }
        public void setAverageTotalTravelTime(double averageTotalTravelTime) { this.averageTotalTravelTime = averageTotalTravelTime; }
        public double getPercentile5thTotalTravelTime() { return percentile5thTotalTravelTime; }
        public void setPercentile5thTotalTravelTime(double percentile5thTotalTravelTime) { this.percentile5thTotalTravelTime = percentile5thTotalTravelTime; }
        public double getPercentile95thTotalTravelTime() { return percentile95thTotalTravelTime; }
        public void setPercentile95thTotalTravelTime(double percentile95thTotalTravelTime) { this.percentile95thTotalTravelTime = percentile95thTotalTravelTime; }

        // setArrivalRedirectionRate, setDepartureRedirectionRate, setTotalTravelTime, setAssignedAccessStation, setAssignedEgressStation
        public void setTravelTimeChanges(String tripId, Double travelTimeChanges) {
            this.travelTimeChanges.put(tripId, travelTimeChanges);
        }
        public void setFlightDistanceChanges(String tripId, Double flightDistanceChanges) {
            this.flightDistanceChanges.put(tripId, flightDistanceChanges);
        }
        public void setDepartureRedirectionRate(String tripId, Double departureRedirectionRates) {
            this.departureRedirectionRates.put(tripId, departureRedirectionRates);
        }
        public void setArrivalRedirectionRate(String tripId, Double arrivalRedirectionRates) {
            this.arrivalRedirectionRates.put(tripId, arrivalRedirectionRates);
        }
        public void setTotalTravelTime(String tripId, Double totalTravelTimes) {
            this.totalTravelTimes.put(tripId, totalTravelTimes);
        }
        public void setAssignedAccessStation(String tripId, String assignedAccessStations) {
            this.assignedAccessStations.put(tripId, assignedAccessStations);
        }
        public void setAssignedEgressStation(String tripId, String assignedEgressStations) {
            this.assignedEgressStations.put(tripId, assignedEgressStations);
        }

        public double getUamVehicleMeter() { return uamVehicleMeter; }
        public void setUamVehicleMeter(double uamVehicleKilometer) { this.uamVehicleMeter = uamVehicleKilometer; }
        public int getNumberOfUAMVehiclesUsed() { return numberOfUAMVehiclesUsed; }
        public void setNumberOfUAMVehiclesUsed(int numberOfUAMVehiclesUsed) { this.numberOfUAMVehiclesUsed = numberOfUAMVehiclesUsed; }

        public void setDeadheadingFlightDistance(double deadHeadingDistance) {
            this.deadHeadingFlightDistance = deadHeadingDistance;
        }
        public double getDeadHeadingFlightDistance() {
            return this.deadHeadingFlightDistance;
        }
        public void setFleetSize(int fleetSize) {
            this.fleetSize = fleetSize;
        }
        public int getFleetSize() {
            return this.fleetSize;
        }
        public void setVehicleAssignments(Map<Integer, List<TripItemForOptimization>> vehicleAssignments) {
            this.vehicleAssignments = vehicleAssignments;
        }
        public Map<Integer, List<TripItemForOptimization>> getVehicleAssignments() {
            return this.vehicleAssignments;
        }
    }
    private SolutionFitnessPair calculatePopulationIndicators(List<SolutionFitnessPair> population) {
        List<SolutionIndicatorData> indicatorDataList = new ArrayList<>();
        SolutionFitnessPair bestSolution = null;
        double bestFitness = Double.NEGATIVE_INFINITY;

        for (SolutionFitnessPair solutionPair : population) {
            SolutionIndicatorData indicatorData = new SolutionIndicatorData(solutionPair.getSolution());
            calculateFitness(solutionPair.getSolution(), indicatorData, true);
            calculateAdditionalIndicators(indicatorData);

            // Store the indicatorData in the solution pair
            solutionPair.setIndicatorData(indicatorData);
            indicatorDataList.add(indicatorData);

            // Track the best solution based on the first fitness objective
            if (bestSolution == null || solutionPair.getFitness()[0] > bestFitness) {
                bestSolution = solutionPair;
                bestFitness = solutionPair.getFitness()[0];
            }
        }

        // Write indicators to CSV
        writeIndicatorsToCsv(indicatorDataList, outputSubFolder + "last_iteration_solutions_indicators.csv");

        return bestSolution;
    }
    private void calculateAdditionalIndicators(SolutionIndicatorData indicatorData) {
        List<Double> travelTimeChanges = new ArrayList<>(indicatorData.getTravelTimeChanges().values());
        List<Double> flightDistanceChanges = new ArrayList<>(indicatorData.getFlightDistanceChanges().values());
        List<Double> departureRedirectionRates = new ArrayList<>(indicatorData.getDepartureRedirectionRates().values());
        List<Double> arrivalRedirectionRates = new ArrayList<>(indicatorData.getArrivalRedirectionRates().values());
        List<Double> totalTravelTimes = new ArrayList<>(indicatorData.getTotalTravelTimes().values());

        // Sort lists for percentile calculations
        Collections.sort(travelTimeChanges);
        Collections.sort(flightDistanceChanges);
        Collections.sort(departureRedirectionRates);
        Collections.sort(arrivalRedirectionRates);
        Collections.sort(totalTravelTimes);

        // Calculate averages
        indicatorData.setAverageTravelTimeChange(calculateAverage(travelTimeChanges));
        indicatorData.setAverageFlightDistanceChange(calculateAverage(flightDistanceChanges));
        indicatorData.setAverageDepartureRedirectionRate(calculateAverage(departureRedirectionRates));
        indicatorData.setAverageArrivalRedirectionRate(calculateAverage(arrivalRedirectionRates));
        indicatorData.setAverageTotalTravelTime(calculateAverage(totalTravelTimes));

        // Calculate percentiles
        indicatorData.setPercentile5thTravelTimeChange(calculatePercentile(travelTimeChanges, 5));
        indicatorData.setPercentile95thTravelTimeChange(calculatePercentile(travelTimeChanges, 95));
        indicatorData.setPercentile5thFlightDistanceChange(calculatePercentile(flightDistanceChanges, 5));
        indicatorData.setPercentile95thFlightDistanceChange(calculatePercentile(flightDistanceChanges, 95));
        indicatorData.setPercentile5thDepartureRedirectionRate(calculatePercentile(departureRedirectionRates, 5));
        indicatorData.setPercentile95thDepartureRedirectionRate(calculatePercentile(departureRedirectionRates, 95));
        indicatorData.setPercentile5thArrivalRedirectionRate(calculatePercentile(arrivalRedirectionRates, 5));
        indicatorData.setPercentile95thArrivalRedirectionRate(calculatePercentile(arrivalRedirectionRates, 95));
        indicatorData.setPercentile5thTotalTravelTime(calculatePercentile(totalTravelTimes, 5));
        indicatorData.setPercentile95thTotalTravelTime(calculatePercentile(totalTravelTimes, 95));

        // Calculate deadheading distance and fleet size using shareability network
        UAMOptimizationController optimizer = new UAMOptimizationController(
                indicatorData.getVehicleAssignments(),  // Your vehicle assignments map
                MAX_DETOUR_RATIO,                // maxDetourRatio
                VEHICLE_CAPACITY,   // maxPassengersPerVehicle
                MAX_CONNECTION_TIME_MINUTES,                 // maxConnectionTimeMinutes
                VEHICLE_CRUISE_SPEED, // flightSpeedMetersPerSecond
                vehicleOriginStationMap,
                vehicleDestinationStationMap
        );
        OptimizationResult result = optimizer.optimize();
        double deadheadingDistance = result.getTotalDeadheadingFlightDistance();
        indicatorData.setDeadheadingFlightDistance(deadheadingDistance);
        indicatorData.setFleetSize(result.getFleetSize());
    }
    private double calculateAverage(List<Double> values) {
        return values.isEmpty() ? Double.NaN : values.stream().mapToDouble(Double::doubleValue).average().orElse(Double.NaN);
    }
    private double calculatePercentile(List<Double> sortedValues, int percentile) {
        if (sortedValues.isEmpty()) {
            return Double.NaN;
        }
        int index = (int) Math.ceil(percentile / 100.0 * sortedValues.size()) - 1;
        return sortedValues.get(Math.max(0, Math.min(sortedValues.size() - 1, index)));
    }
    private void writeIndicatorsToCsv(List<SolutionIndicatorData> indicatorDataList, String fileName) {
        try (FileWriter writer = new FileWriter(fileName)) {
            // Write header
            writer.append("TotalFitness,TotalFlightDistanceChange,TotalTravelTimeChange,TotalCapacityViolationPenalty,PoolingRate,Capacity0Rate,Capacity1Rate,Capacity2Rate,Capacity3Rate,Capacity4Rate,SharedRidesExceedingThresholdRate,TotalSharedRidesExceedingThresholdRate,AvgTravelTimeChange,5thPercentileTravelTimeChange,95thPercentileTravelTimeChange,AvgFlightDistanceChange,5thPercentileFlightDistanceChange,95thPercentileFlightDistanceChange,AvgDepartureRedirectionRate,5thPercentileDepartureRedirectionRate,95thPercentileDepartureRedirectionRate,AvgArrivalRedirectionRate,5thPercentileArrivalRedirectionRate,95thPercentileArrivalRedirectionRate,AvgTotalTravelTime,5thPercentileTotalTravelTime,95thPercentileTotalTravelTime,TotalVehicleMeter,NumberOfVehiclesUsed,DeadheadingFlightDistanceChange,FleetSizeChange\n");

            // Write data for each solution
            for (SolutionIndicatorData data : indicatorDataList) {
                writer.append(String.format("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%f,%d\n",
                        data.getFitness()[0], REVERT_SIGN * data.getFitness()[1], REVERT_SIGN * data.getFitness()[2], data.getFitness()[3],
                        data.getPoolingRate(),
                        data.getVehicleCapacityRates().getOrDefault(0, 0.0),
                        data.getVehicleCapacityRates().getOrDefault(1, 0.0),
                        data.getVehicleCapacityRates().getOrDefault(2, 0.0),
                        data.getVehicleCapacityRates().getOrDefault(3, 0.0),
                        data.getVehicleCapacityRates().getOrDefault(4, 0.0),
                        data.getSharedRidesExceedingThresholdRate(),
                        data.getTotalSharedRidesExceedingThresholdRate(),
                        data.getAverageTravelTimeChange(),
                        data.getPercentile5thTravelTimeChange(),
                        data.getPercentile95thTravelTimeChange(),
                        data.getAverageFlightDistanceChange(),
                        data.getPercentile5thFlightDistanceChange(),
                        data.getPercentile95thFlightDistanceChange(),
                        data.getAverageDepartureRedirectionRate(),
                        data.getPercentile5thDepartureRedirectionRate(),
                        data.getPercentile95thDepartureRedirectionRate(),
                        data.getAverageArrivalRedirectionRate(),
                        data.getPercentile5thArrivalRedirectionRate(),
                        data.getPercentile95thArrivalRedirectionRate(),
                        data.getAverageTotalTravelTime(),
                        data.getPercentile5thTotalTravelTime(),
                        data.getPercentile95thTotalTravelTime(),
                        data.getUamVehicleMeter(),
                        data.getNumberOfUAMVehiclesUsed(), // This is actually the number of UAM vehicle-operations
                        data.getDeadHeadingFlightDistance()-getNonPooledDeadheadingDistance(),
                        data.getFleetSize() - getNonPooledFleetSize()
                ));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    // Initial data extraction methods =================================================================================
/*    private ArrayList<UAMTrip> extractSubTrips(List<UAMTrip> uamTrips) {
        // extract sub trips from uamTrips based on the departure time of trips falling between buffer start and end time
        return uamTrips.stream()
                .filter(trip -> trip.getDepartureTime() >= BUFFER_START_TIME && trip.getDepartureTime() < BUFFER_END_TIME)
                .collect(Collectors.toCollection(ArrayList::new));
    }*/

    // Method to create UAM vehicles and assign them to stations in the initialization phase (could also be used in later stage)
    private void saveStationVehicleNumber(List<TripItemForOptimization> subTrips) {
        // save the station's vehicle number for the current time based on the UAMTrips' origin and destination station
        for (TripItemForOptimization subTrip : subTrips) {
            feedDataForVehicleCreation(subTrip, true);
        }
    }

    private UAMVehicle feedDataForVehicleCreation(TripItemForOptimization subTrip, boolean isAddingVehicleBeforeInitialization) {
        Vertiport nearestOriginStation = findNearestStation(subTrip, vertiportsMap, true);
        Vertiport nearestDestinationStation = findNearestStation(subTrip, vertiportsMap, false);

        if (nearestOriginStation == null || nearestDestinationStation == null) {
            log.error("Found null station for trip: " + subTrip.tripID);
        }

        UAMVehicle vehicle = createVehicle(nearestOriginStation);

        //vehicles.put(vehicle.getId(), vehicle);
        vehicleOriginStationMap.put(vehicle.getId(), nearestOriginStation);
        vehicleDestinationStationMap.put(vehicle.getId(), nearestDestinationStation);
        //vehicleOccupancyMap.put(vehicle, VEHICLE_CAPACITY);

        if (isAddingVehicleBeforeInitialization){
            // Get the station ID
            Id<Vertiport> nearestOriginStationId = Id.create(nearestOriginStation.ID, Vertiport.class);
            // Check if there is already a list for this station ID, if not, create one
            List<UAMVehicle> vehiclesAtStation = originStationVehicleMap.computeIfAbsent(nearestOriginStationId, k -> new ArrayList<>());
            // Add the new vehicle to the list
            vehiclesAtStation.add(vehicle);
            originStationVehicleMap.put(nearestOriginStationId, vehiclesAtStation);
        }
        return vehicle;
    }
    // vehicle creator function
    private UAMVehicle createVehicle(Vertiport uamStation) {
        /*        UAMVehicleType vehicleType = new UAMVehicleType(id, capacity, range, horizontalSpeed, verticalSpeed,
                boardingTime, deboardingTime, turnAroundTime, energyConsumptionVertical, energyConsumptionHorizontal,
                maximumCharge);*/
        //final Map<Id<UAMVehicleType>, UAMVehicleType> vehicleTypes = new HashMap<>();
        Id<UAMVehicleType> vehicleTypeId = Id.create("poolingVehicle", UAMVehicleType.class);
        UAMVehicleType vehicleType = new UAMVehicleType(vehicleTypeId, 0, 0, 0, 0,
                0, 0, 0);
        //vehicleTypes.put(vehicleTypeId, vehicleType);

        // Create a builder instance
        ImmutableDvrpVehicleSpecification.Builder builder = ImmutableDvrpVehicleSpecification.newBuilder();
        // Set the properties of the vehicle
        builder.id(Id.create(String.valueOf(FIRST_UAM_VEHICLE_ID++), DvrpVehicle.class));
        builder.startLinkId(Id.create("0001", Link.class));
        builder.capacity(VEHICLE_CAPACITY);
        builder.serviceBeginTime(BUFFER_START_TIME);
        builder.serviceEndTime(END_SERVICE_TIME_OF_THE_DAY);
        // Build the vehicle specification
        ImmutableDvrpVehicleSpecification vehicleSpecification = builder.build();

        return new UAMVehicle(vehicleSpecification,
                network.getLinks().get(Id.create("0001", Link.class)), // TODO: need do it more elegantly
                Id.create(uamStation.ID, UAMStation.class),
                vehicleType);
    }

    private static Vertiport findNearestStation(TripItemForOptimization trip, HashMap<Integer, Vertiport> vertiportsMap, boolean accessLeg) {
        Vertiport nearestStation = null;
        double shortestDistance = Double.MAX_VALUE;
        for (Vertiport station : vertiportsMap.values()) {
            if (station == null) {
                log.error("Encountered null station in stations map");
                continue; // Skip null stations
            }
            double distance = Double.MAX_VALUE;
            if (accessLeg) {
                if (trip.originNeighborVertiportCandidatesTimeAndDistance.containsKey(station)) {
                    distance = trip.originNeighborVertiportCandidatesTimeAndDistance.get(station).get("distance");
                }
            } else {
                if (trip.destinationNeighborVertiportCandidatesTimeAndDistance.containsKey(station)) {
                    distance = trip.destinationNeighborVertiportCandidatesTimeAndDistance.get(station).get("distance");
                }
            }
            if (distance < shortestDistance) {
                nearestStation = station;
                shortestDistance = distance;
            }
        }
        if (nearestStation == null) {
            log.warn("No nearest station found for trip: " + trip.tripID);
        }
        return nearestStation;
    }
    private static Vertiport findEuclideanNearestStation(TripItemForOptimization trip, HashMap<Integer, Vertiport> vertiportsMap, boolean accessLeg) {
        Vertiport nearestStation = null;
        double shortestDistance = Double.MAX_VALUE;
        for (Vertiport station : vertiportsMap.values()) {
            if (station == null) {
                log.error("Encountered null station in stations map");
                continue; // Skip null stations
            }
            double distance = Double.MAX_VALUE;
            if (accessLeg) {
                    distance = calculateEuclideanDistance(trip.origin, station.coord);
            } else {
                distance = calculateEuclideanDistance(trip.destination, station.coord);
            }
            if (distance < shortestDistance) {
                nearestStation = station;
                shortestDistance = distance;
            }
        }
        if (nearestStation == null) {
            log.warn("No nearest station found for trip: " + trip.tripID);
        }
        return nearestStation;
    }

}
