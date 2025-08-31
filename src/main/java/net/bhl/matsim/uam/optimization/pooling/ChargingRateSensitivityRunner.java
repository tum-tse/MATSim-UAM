package net.bhl.matsim.uam.optimization.pooling;

import java.io.IOException;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.List;
import java.util.ArrayList;

import static net.bhl.matsim.uam.optimization.pooling.GridSearch.TIMEOUT_MINUTES;
import static net.bhl.matsim.uam.optimization.pooling.MultiObjectiveNSGAII.setFilePaths;

/**
 * Runner class for conducting sensitivity analysis on charging rates (CHARGING_RATE_KWH_PER_SECOND parameter)
 * with fixed pooling parameters: pooling time window = 3 minutes, origin search radius = 6000m, destination search radius = 6000m
 */
public class ChargingRateSensitivityRunner {
    
    // Fixed pooling parameters as requested (obtained from SensitivityConfig)
    private static final boolean ENABLE_LOCAL_SEARCH = true;
    private static final boolean ENABLE_PRINT_RESULTS = true;
    
    // Charging rate values to test (kWh/second)
    // Based on the GuideForSensitivityAnalysis.md document
    private static final double[] CHARGING_RATES_KWH_PER_SECOND = {
        2.08 / 60.0,   // Very slow charging (1C) - 2.08 kWh/min = 0.0347 kWh/s
        4.17 / 60.0,   // Slow charging (2C) - 4.17 kWh/min = 0.0695 kWh/s
        6.24 / 60.0,   // Default charging rate (3C) - 6.24 kWh/min = 0.104 kWh/s
        8.34 / 60.0,   // Fast charging (4C) - 8.34 kWh/min = 0.139 kWh/s
        10.42 / 60.0,  // Very fast charging (5C) - 10.42 kWh/min = 0.1737 kWh/s
        12.48 / 60.0,  // Ultra fast charging (6C) - 12.48 kWh/min = 0.208 kWh/s
        14.56 / 60.0,  // Extreme charging (7C) - 14.56 kWh/min = 0.243 kWh/s
        16.68 / 60.0   // Maximum charging (8C) - 16.68 kWh/min = 0.278 kWh/s
    };
    
    private static final String[] CHARGING_RATE_LABELS = {
        "1c", "2c", "3c", "4c",
        "5c", "6c", "7c", "8c"
    };
    
    public static void main(String[] args) throws IOException, InterruptedException {
        // Initialize and run the optimization
        MultiObjectiveNSGAII.initialization(args);

        if (args.length < 5) {
            System.out.println("Usage: ChargingRateSensitivityRunner <Trip_Item> <Config> <Vertiport_Unit_Candidate> <Scenario_Configuration> <Result_Output>");
            System.exit(1);
        }
        
        String tripItemFile = args[0];
        String configFile = args[1];
        String vertiportFile = args[2];
        String scenarioFile = args[3];
        String baseOutputDir = args[4];

        // Create directory for this sensitivity analysis
        String sensitivityOutputDir = baseOutputDir + "/charging_rate/";
        MultiObjectiveNSGAII.createFolder(sensitivityOutputDir);
        setFilePaths(args[0], args[1], args[2], args[3], sensitivityOutputDir);
        
        // Get fixed parameters from SensitivityConfig
        SensitivityConfig defaultConfig = new SensitivityConfig();
        
        System.out.println("Starting Charging Rate Sensitivity Analysis...");
        System.out.println("Fixed parameters:");
        System.out.println("  Pooling time window: " + defaultConfig.getPoolingTimeWindow() + " minutes");
        System.out.println("  Origin search radius: " + defaultConfig.getOriginSearchRadius() + " meters");
        System.out.println("  Destination search radius: " + defaultConfig.getDestinationSearchRadius() + " meters");
        System.out.println("  Number of Monte Carlo simulations: " + defaultConfig.getNumSimulations());
        System.out.println("Output directory: " + sensitivityOutputDir);
        System.out.println("Testing charging rates:");
        for (int i = 0; i < CHARGING_RATES_KWH_PER_SECOND.length; i++) {
            System.out.printf("  %s: %.4f kWh/s (%.2f kWh/min)%n", 
                CHARGING_RATE_LABELS[i], 
                CHARGING_RATES_KWH_PER_SECOND[i], 
                CHARGING_RATES_KWH_PER_SECOND[i] * 60);
        }
        
        // Run experiments in parallel
        runParallelExperiments(tripItemFile, configFile, vertiportFile, scenarioFile, sensitivityOutputDir);
        
        System.out.println("Charging Rate Sensitivity Analysis completed!");
        System.out.println("Results saved in: " + sensitivityOutputDir);
    }
    
    private static void runParallelExperiments(String tripItemFile, String configFile, 
                                             String vertiportFile, String scenarioFile, String outputDir) {
        int numThreads = Math.min(CHARGING_RATES_KWH_PER_SECOND.length, Runtime.getRuntime().availableProcessors());
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);
        List<Future<?>> futures = new ArrayList<>();
        
        for (int i = 0; i < CHARGING_RATES_KWH_PER_SECOND.length; i++) {
            final double chargingRate = CHARGING_RATES_KWH_PER_SECOND[i];
            final String label = CHARGING_RATE_LABELS[i];
            
            Future<?> future = executor.submit(() -> runSingleExperiment(
                tripItemFile, configFile, vertiportFile, scenarioFile, outputDir, chargingRate, label));
            futures.add(future);
        }
        
        // Wait for all experiments to complete
        for (Future<?> future : futures) {
            try {
                future.get(); // This will block until the task is complete
            } catch (Exception e) {
                System.err.println("Error in experiment execution: " + e.getMessage());
                e.printStackTrace();
            }
        }
        
        executor.shutdown();
        try {
            if (!executor.awaitTermination(TIMEOUT_MINUTES, TimeUnit.MINUTES)) {
                executor.shutdownNow();
            }
        } catch (InterruptedException e) {
            executor.shutdownNow();
            Thread.currentThread().interrupt();
        }
    }
    
    private static void runSingleExperiment(String tripItemFile, String configFile, 
                                          String vertiportFile, String scenarioFile, 
                                          String baseOutputDir, double chargingRate, String label) {
        try {
            // Create SensitivityConfig object and convert to string format
            SensitivityConfig config = SensitivityConfig.forChargingRateAnalysis(chargingRate);
            String sensitivityConfig = config.getNumSimulations() + "," + config.getChargingRateKwhPerSecond();
            
            // Build arguments array for MultiObjectiveNSGAII using SensitivityConfig values
            String[] optimizationArgs = {
                tripItemFile,                    // Trip_Item
                configFile,                      // Config
                vertiportFile,                   // Vertiport_Unit_Candidate
                scenarioFile,                    // Scenario_Configuration
                baseOutputDir,                   // Result_Output
                String.valueOf(config.getPoolingTimeWindow()),  // BUFFER_END_TIME (convert seconds to minutes)
                String.valueOf(config.getOriginSearchRadius()), // SEARCH_RADIUS_ORIGIN
                String.valueOf(config.getDestinationSearchRadius()), // SEARCH_RADIUS_DESTINATION
                String.valueOf(ENABLE_LOCAL_SEARCH),  // ENABLE_LOCAL_SEARCH
                String.valueOf(ENABLE_PRINT_RESULTS), // ENABLE_PRINT_RESULTS
                label + "/",                               // OUTPUT_SUBFOLDER
                sensitivityConfig                     // SENSITIVITY_CONFIG
            };
            
            System.out.printf("Running experiment with charging rate: %.4f kWh/s (%.2f kWh/min)...%n", 
                             chargingRate, chargingRate * 60);

            double[] results = MultiObjectiveNSGAII.callAlgorithm(optimizationArgs);
            
            System.out.printf("Completed experiment with charging rate %.4f kWh/s. Best fitness: %.6f%n", 
                             chargingRate, results[3]);
            
        } catch (Exception e) {
            System.err.printf("Error running experiment with charging rate %.4f kWh/s: %s%n", 
                             chargingRate, e.getMessage());
            e.printStackTrace();
        }
    }
}