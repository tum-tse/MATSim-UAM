package net.bhl.matsim.uam.optimization.pooling.sn;

import net.bhl.matsim.uam.optimization.pooling.MultiObjectiveNSGAII;
import org.matsim.api.core.v01.Coord;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.*;
import java.util.stream.Collectors;

import static net.bhl.matsim.uam.optimization.pooling.sn.EVTOLBatteryManager.DEFAULT_CHARGING_RATE_KWH_PER_SECOND;

/**
 * Enhanced battery-aware eVTOL optimization system with Monte Carlo analysis
 * Runs multiple simulations with different random seeds and averages results
 */
public class BatteryAwareOptimizationExample {

    // Essential constants for extreme stress test
    private static final int VEHICLE_CAPACITY = MultiObjectiveNSGAII.VEHICLE_CAPACITY;
    private static final double VEHICLE_CRUISE_SPEED = MultiObjectiveNSGAII.VEHICLE_CRUISE_SPEED;
    private static final double MAX_DETOUR_RATIO = MultiObjectiveNSGAII.MAX_DETOUR_RATIO;
    private static final int MAX_CONNECTION_TIME_MINUTES = MultiObjectiveNSGAII.MAX_CONNECTION_TIME_MINUTES;

    // Monte Carlo simulation parameters
    private static final int DEFAULT_NUM_RUNS = 10; // Number of Monte Carlo runs
    private static final long BASE_SEED = 54321; // Base seed for reproducibility

    // Only the charging rates used in extreme stress test
    private static final double[] CRITICAL_CHARGING_RATES = {
            0 / 60.0,       // No charging (0C)
            4.17 / 60.0,    // Slow charging (2C)
            6.24 / 60.0,    // Default charging rate (3C)
            8.34 / 60.0,    // Fast charging (4C)
            10.42 / 60.0,   // Very fast charging (5C)
            16.68 / 60.0    // Maximum charging (8C)
    };

    private static final String[] CRITICAL_RATE_NAMES = {
            "No charging", "Slow charging", "Default Rate", "Fast charging", "Very fast charging", "Maximum charging"
    };

    public static void main(String[] args) {
        System.out.println("=".repeat(80));
        System.out.println("eVTOL MONTE CARLO ANALYSIS");
        System.out.println("=".repeat(80));

        int numRuns = DEFAULT_NUM_RUNS;
        if (args.length > 0) {
            try {
                numRuns = Integer.parseInt(args[0]);
                System.out.println("Using " + numRuns + " Monte Carlo runs");
            } catch (NumberFormatException e) {
                System.out.println("Invalid number of runs specified, using default: " + DEFAULT_NUM_RUNS);
            }
        }

        runMonteCarloAnalysis(numRuns);

        System.out.println("\n" + "=".repeat(80));
        System.out.println("MONTE CARLO ANALYSIS COMPLETED!");
        System.out.println("=".repeat(80));
    }

    /**
     * Run Monte Carlo analysis with multiple random seeds
     */
    private static void runMonteCarloAnalysis(int numRuns) {
        System.out.println("=".repeat(80));
        System.out.println("MONTE CARLO ANALYSIS: " + numRuns + " SIMULATION RUNS");
        System.out.println("=".repeat(80));

        // Use thread pool for parallel execution
        int numThreads = Math.min(numRuns, Runtime.getRuntime().availableProcessors());
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        try {
            // Submit all runs for parallel execution
            List<Future<SimulationRunResult>> futures = new ArrayList<>();

            for (int run = 0; run < numRuns; run++) {
                final int runIndex = run;
                futures.add(executor.submit(() -> runSingleSimulation(runIndex)));
            }

            // Collect results
            List<SimulationRunResult> results = new ArrayList<>();
            for (int i = 0; i < futures.size(); i++) {
                try {
                    SimulationRunResult result = futures.get(i).get(30, TimeUnit.MINUTES);
                    results.add(result);
                    System.out.printf("Completed run %d/%d\n", i + 1, numRuns);
                } catch (Exception e) {
                    System.err.printf("Run %d failed: %s\n", i + 1, e.getMessage());
                }
            }

            // Analyze and present results
            analyzeMonteCarloResults(results);

        } finally {
            executor.shutdown();
        }
    }

    /**
     * Run a single simulation with a unique random seed
     */
    private static SimulationRunResult runSingleSimulation(int runIndex) {
        long seed = BASE_SEED + runIndex;
        System.out.printf("Starting run %d with seed %d\n", runIndex + 1, seed);

        // Create trips with unique seed
        List<VehicleTrip> trips = createExtremeStressTestTrips(seed);

        // Test all charging rates
        List<ChargingRateResult> chargingResults = new ArrayList<>();

        for (int i = 0; i < CRITICAL_CHARGING_RATES.length; i++) {
            try {
                long startTime = System.currentTimeMillis();

                UAMOptimizationController optimizer = new UAMOptimizationController(
                        trips,
                        MAX_DETOUR_RATIO,
                        VEHICLE_CAPACITY,
                        MAX_CONNECTION_TIME_MINUTES,
                        VEHICLE_CRUISE_SPEED,
                        CRITICAL_CHARGING_RATES[i],
                        false
                );

                OptimizationResult result = optimizer.optimize();
                long duration = System.currentTimeMillis() - startTime;

                // Extract performance metrics
                OptimizationResult.PerformanceMetrics metrics = result.getPerformanceMetrics();

                ChargingRateResult chargingResult = new ChargingRateResult(
                        CRITICAL_RATE_NAMES[i],
                        metrics.fleetSize,
                        metrics.vtolOperations,
                        metrics.totalEnergyUsed,
                        metrics.averageEnergyPerVehicle,
                        metrics.totalChargingTime,
                        metrics.averageChargingTimePerVehicle,
                        duration / 1000.0,
                        true
                );

                chargingResults.add(chargingResult);

            } catch (Exception e) {
                // Add failed result
                ChargingRateResult failedResult = new ChargingRateResult(
                        CRITICAL_RATE_NAMES[i],
                        0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, false
                );
                chargingResults.add(failedResult);
                System.err.printf("Run %d, charging rate %s failed: %s\n",
                        runIndex + 1, CRITICAL_RATE_NAMES[i], e.getMessage());
            }
        }

        return new SimulationRunResult(runIndex + 1, seed, chargingResults);
    }

    /**
     * Analyze Monte Carlo results and present averaged statistics
     */
    private static void analyzeMonteCarloResults(List<SimulationRunResult> results) {
        if (results.isEmpty()) {
            System.out.println("No successful simulation runs to analyze!");
            return;
        }

        System.out.println("\n" + "=".repeat(80));
        System.out.printf("MONTE CARLO RESULTS ANALYSIS (%d successful runs)\n", results.size());
        System.out.println("=".repeat(80));

        // Calculate statistics for each charging rate
        for (int rateIndex = 0; rateIndex < CRITICAL_RATE_NAMES.length; rateIndex++) {
            final int currentRateIndex = rateIndex;

            // Extract results for this charging rate from all runs
            List<ChargingRateResult> rateResults = results.stream()
                    .map(run -> run.chargingResults.get(currentRateIndex))
                    .filter(result -> result.success)
                    .collect(Collectors.toList());

            if (rateResults.isEmpty()) {
                System.out.printf("\n%-18s: NO SUCCESSFUL RUNS\n", CRITICAL_RATE_NAMES[rateIndex]);
                continue;
            }

            // Calculate statistics
            ChargingeRateStatistics stats = calculateChargingRateStatistics(rateResults);

            System.out.printf("\n%-18s: (%d successful runs)\n", CRITICAL_RATE_NAMES[rateIndex], rateResults.size());
            System.out.printf("  Fleet Size:          %.1f ± %.1f (range: %.0f - %.0f)\n",
                    stats.avgFleetSize, stats.stdFleetSize, stats.minFleetSize, stats.maxFleetSize);
            System.out.printf("  VTOL Operations:     %.1f ± %.1f (range: %.0f - %.0f)\n",
                    stats.avgVtolOps, stats.stdVtolOps, stats.minVtolOps, stats.maxVtolOps);
            System.out.printf("  Total Energy (kWh):  %.2f ± %.2f (range: %.2f - %.2f)\n",
                    stats.avgTotalEnergy, stats.stdTotalEnergy, stats.minTotalEnergy, stats.maxTotalEnergy);
            System.out.printf("  Avg Energy/Vehicle:  %.2f ± %.2f kWh\n",
                    stats.avgEnergyPerVehicle, stats.stdEnergyPerVehicle);
            System.out.printf("  Total Charging (min): %.1f ± %.1f (range: %.1f - %.1f)\n",
                    stats.avgTotalCharging, stats.stdTotalCharging, stats.minTotalCharging, stats.maxTotalCharging);
            System.out.printf("  Avg Charging/Vehicle: %.1f ± %.1f min\n",
                    stats.avgChargingPerVehicle, stats.stdChargingPerVehicle);
            System.out.printf("  Execution Time:      %.1f ± %.1f s\n",
                    stats.avgExecutionTime, stats.stdExecutionTime);

            if (rateResults.size() < results.size()) {
                int failedRuns = results.size() - rateResults.size();
                System.out.printf("  Failed Runs:         %d (%.1f%%)\n",
                        failedRuns, (failedRuns * 100.0) / results.size());
            }
        }

        // Overall summary
        System.out.println("\n" + "=".repeat(80));
        System.out.println("SUMMARY STATISTICS");
        System.out.println("=".repeat(80));

        System.out.printf("Total simulation runs: %d\n", results.size());
        System.out.printf("Trip generation seeds: %d - %d\n",
                results.get(0).seed, results.get(results.size() - 1).seed);

        // Find best performing charging rate on average
        double bestAvgFleetSize = Double.MAX_VALUE;
        String bestChargingRate = "";

        for (int rateIndex = 0; rateIndex < CRITICAL_RATE_NAMES.length; rateIndex++) {
            final int currentRateIndex = rateIndex;

            List<ChargingRateResult> rateResults = results.stream()
                    .map(run -> run.chargingResults.get(currentRateIndex))
                    .filter(result -> result.success)
                    .collect(Collectors.toList());

            if (!rateResults.isEmpty()) {
                double avgFleetSize = rateResults.stream()
                        .mapToDouble(r -> r.fleetSize)
                        .average()
                        .orElse(Double.MAX_VALUE);

                if (avgFleetSize < bestAvgFleetSize) {
                    bestAvgFleetSize = avgFleetSize;
                    bestChargingRate = CRITICAL_RATE_NAMES[rateIndex];
                }
            }
        }

        if (!bestChargingRate.isEmpty()) {
            System.out.printf("Best performing charging strategy: %s (avg fleet size: %.1f)\n",
                    bestChargingRate, bestAvgFleetSize);
        }
    }

    /**
     * Calculate statistics for a charging rate across multiple runs
     */
    private static ChargingeRateStatistics calculateChargingRateStatistics(List<ChargingRateResult> results) {
        double avgFleetSize = results.stream().mapToDouble(r -> r.fleetSize).average().orElse(0);
        double avgVtolOps = results.stream().mapToDouble(r -> r.vtolOperations).average().orElse(0);
        double avgTotalEnergy = results.stream().mapToDouble(r -> r.totalEnergyUsed).average().orElse(0);
        double avgEnergyPerVehicle = results.stream().mapToDouble(r -> r.averageEnergyPerVehicle).average().orElse(0);
        double avgTotalCharging = results.stream().mapToDouble(r -> r.totalChargingTime).average().orElse(0);
        double avgChargingPerVehicle = results.stream().mapToDouble(r -> r.averageChargingTimePerVehicle).average().orElse(0);
        double avgExecutionTime = results.stream().mapToDouble(r -> r.executionTimeSeconds).average().orElse(0);

        // Calculate standard deviations
        double stdFleetSize = calculateStandardDeviation(results.stream().mapToDouble(r -> r.fleetSize).toArray(), avgFleetSize);
        double stdVtolOps = calculateStandardDeviation(results.stream().mapToDouble(r -> r.vtolOperations).toArray(), avgVtolOps);
        double stdTotalEnergy = calculateStandardDeviation(results.stream().mapToDouble(r -> r.totalEnergyUsed).toArray(), avgTotalEnergy);
        double stdEnergyPerVehicle = calculateStandardDeviation(results.stream().mapToDouble(r -> r.averageEnergyPerVehicle).toArray(), avgEnergyPerVehicle);
        double stdTotalCharging = calculateStandardDeviation(results.stream().mapToDouble(r -> r.totalChargingTime).toArray(), avgTotalCharging);
        double stdChargingPerVehicle = calculateStandardDeviation(results.stream().mapToDouble(r -> r.averageChargingTimePerVehicle).toArray(), avgChargingPerVehicle);
        double stdExecutionTime = calculateStandardDeviation(results.stream().mapToDouble(r -> r.executionTimeSeconds).toArray(), avgExecutionTime);

        // Calculate min/max
        double minFleetSize = results.stream().mapToDouble(r -> r.fleetSize).min().orElse(0);
        double maxFleetSize = results.stream().mapToDouble(r -> r.fleetSize).max().orElse(0);
        double minVtolOps = results.stream().mapToDouble(r -> r.vtolOperations).min().orElse(0);
        double maxVtolOps = results.stream().mapToDouble(r -> r.vtolOperations).max().orElse(0);
        double minTotalEnergy = results.stream().mapToDouble(r -> r.totalEnergyUsed).min().orElse(0);
        double maxTotalEnergy = results.stream().mapToDouble(r -> r.totalEnergyUsed).max().orElse(0);
        double minTotalCharging = results.stream().mapToDouble(r -> r.totalChargingTime).min().orElse(0);
        double maxTotalCharging = results.stream().mapToDouble(r -> r.totalChargingTime).max().orElse(0);

        return new ChargingeRateStatistics(
                avgFleetSize, stdFleetSize, minFleetSize, maxFleetSize,
                avgVtolOps, stdVtolOps, minVtolOps, maxVtolOps,
                avgTotalEnergy, stdTotalEnergy, minTotalEnergy, maxTotalEnergy,
                avgEnergyPerVehicle, stdEnergyPerVehicle,
                avgTotalCharging, stdTotalCharging, minTotalCharging, maxTotalCharging,
                avgChargingPerVehicle, stdChargingPerVehicle,
                avgExecutionTime, stdExecutionTime
        );
    }

    /**
     * Calculate standard deviation
     */
    private static double calculateStandardDeviation(double[] values, double mean) {
        if (values.length <= 1) return 0.0;

        double sumSquaredDiffs = 0.0;
        for (double value : values) {
            sumSquaredDiffs += Math.pow(value - mean, 2);
        }
        return Math.sqrt(sumSquaredDiffs / (values.length - 1));
    }

    /**
     * Create extreme stress test with specified random seed
     */
    private static List<VehicleTrip> createExtremeStressTestTrips(long seed) {
        List<VehicleTrip> trips = new ArrayList<>();
        Random random = new Random(seed);

        int baseTime = 28800; // 8:00 AM
        int timeWindow = 3600;  // 60 minutes (extreme density)
        int tripCount = 100;

        for (int i = 1; i <= tripCount; i++) {
            // Create high-density urban scenario with shorter distances
            double centerX = 75000; // City center
            double centerY = 75000;
            double urbanRadius = 10000; // 10km radius urban area

            // Origins clustered around urban centers
            double originAngle = random.nextDouble() * 2 * Math.PI;
            double originRadius = random.nextGaussian() * urbanRadius * 0.3 + urbanRadius * 0.5;
            Coord origin = new Coord(
                    centerX + originRadius * Math.cos(originAngle),
                    centerY + originRadius * Math.sin(originAngle)
            );

            // Destinations with preference for city center and airports
            double destAngle = random.nextDouble() * 2 * Math.PI;
            double destRadius = random.nextGaussian() * urbanRadius * 0.4 + urbanRadius * 0.6;
            Coord destination = new Coord(
                    centerX + destRadius * Math.cos(destAngle),
                    centerY + destRadius * Math.sin(destAngle)
            );

            // Peak hour departure times (concentrated)
            int departureTime = baseTime + random.nextInt(timeWindow);

            // Calculate trip duration based on distance
            double distance = Math.sqrt(Math.pow(destination.getX() - origin.getX(), 2) +
                    Math.pow(destination.getY() - origin.getY(), 2));
            int tripDuration = (int) (distance / VEHICLE_CRUISE_SPEED) + 300 + random.nextInt(300);

            // Higher passenger loads in extreme scenario
            //int passengers = random.nextDouble() < 0.7 ? VEHICLE_CAPACITY - 1 + random.nextInt(2) : 1;

            trips.add(new VehicleTrip("EX" + i,
                    origin,
                    destination,
                    departureTime,
                    departureTime + tripDuration,
                    VEHICLE_CAPACITY));
        }

        return trips;
    }

    /**
     * Data classes for storing results
     */
    private static class SimulationRunResult {
        final int runNumber;
        final long seed;
        final List<ChargingRateResult> chargingResults;

        public SimulationRunResult(int runNumber, long seed, List<ChargingRateResult> chargingResults) {
            this.runNumber = runNumber;
            this.seed = seed;
            this.chargingResults = chargingResults;
        }
    }

    private static class ChargingRateResult {
        final String chargingRateName;
        final int fleetSize;
        final int vtolOperations;
        final double totalEnergyUsed;
        final double averageEnergyPerVehicle;
        final double totalChargingTime;
        final double averageChargingTimePerVehicle;
        final double executionTimeSeconds;
        final boolean success;

        public ChargingRateResult(String chargingRateName, int fleetSize, int vtolOperations,
                                  double totalEnergyUsed, double averageEnergyPerVehicle,
                                  double totalChargingTime, double averageChargingTimePerVehicle,
                                  double executionTimeSeconds, boolean success) {
            this.chargingRateName = chargingRateName;
            this.fleetSize = fleetSize;
            this.vtolOperations = vtolOperations;
            this.totalEnergyUsed = totalEnergyUsed;
            this.averageEnergyPerVehicle = averageEnergyPerVehicle;
            this.totalChargingTime = totalChargingTime;
            this.averageChargingTimePerVehicle = averageChargingTimePerVehicle;
            this.executionTimeSeconds = executionTimeSeconds;
            this.success = success;
        }
    }

    private static class ChargingeRateStatistics {
        final double avgFleetSize, stdFleetSize, minFleetSize, maxFleetSize;
        final double avgVtolOps, stdVtolOps, minVtolOps, maxVtolOps;
        final double avgTotalEnergy, stdTotalEnergy, minTotalEnergy, maxTotalEnergy;
        final double avgEnergyPerVehicle, stdEnergyPerVehicle;
        final double avgTotalCharging, stdTotalCharging, minTotalCharging, maxTotalCharging;
        final double avgChargingPerVehicle, stdChargingPerVehicle;
        final double avgExecutionTime, stdExecutionTime;

        public ChargingeRateStatistics(double avgFleetSize, double stdFleetSize, double minFleetSize, double maxFleetSize,
                                       double avgVtolOps, double stdVtolOps, double minVtolOps, double maxVtolOps,
                                       double avgTotalEnergy, double stdTotalEnergy, double minTotalEnergy, double maxTotalEnergy,
                                       double avgEnergyPerVehicle, double stdEnergyPerVehicle,
                                       double avgTotalCharging, double stdTotalCharging, double minTotalCharging, double maxTotalCharging,
                                       double avgChargingPerVehicle, double stdChargingPerVehicle,
                                       double avgExecutionTime, double stdExecutionTime) {
            this.avgFleetSize = avgFleetSize;
            this.stdFleetSize = stdFleetSize;
            this.minFleetSize = minFleetSize;
            this.maxFleetSize = maxFleetSize;
            this.avgVtolOps = avgVtolOps;
            this.stdVtolOps = stdVtolOps;
            this.minVtolOps = minVtolOps;
            this.maxVtolOps = maxVtolOps;
            this.avgTotalEnergy = avgTotalEnergy;
            this.stdTotalEnergy = stdTotalEnergy;
            this.minTotalEnergy = minTotalEnergy;
            this.maxTotalEnergy = maxTotalEnergy;
            this.avgEnergyPerVehicle = avgEnergyPerVehicle;
            this.stdEnergyPerVehicle = stdEnergyPerVehicle;
            this.avgTotalCharging = avgTotalCharging;
            this.stdTotalCharging = stdTotalCharging;
            this.minTotalCharging = minTotalCharging;
            this.maxTotalCharging = maxTotalCharging;
            this.avgChargingPerVehicle = avgChargingPerVehicle;
            this.stdChargingPerVehicle = stdChargingPerVehicle;
            this.avgExecutionTime = avgExecutionTime;
            this.stdExecutionTime = stdExecutionTime;
        }
    }
}