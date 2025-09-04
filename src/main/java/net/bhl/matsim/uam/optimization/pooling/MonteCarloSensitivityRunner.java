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
 * Runner class for conducting sensitivity analysis on Monte Carlo simulations (NUM_SIMULATIONS parameter)
 * with fixed pooling parameters: pooling time window = 3 minutes, origin search radius = 6000m, destination search radius = 6000m
 */
public class MonteCarloSensitivityRunner {
    
    // Fixed pooling parameters as requested (obtained from SensitivityConfig)
    private static final boolean ENABLE_LOCAL_SEARCH = true;
    private static final boolean ENABLE_PRINT_RESULTS = true;
    
    // Monte Carlo simulation values to test
    private static final int[] NUM_SIMULATIONS_VALUES = {25, 50, 100, 200, 400, 800, 1600, 3200, 6400, 12800};
    
    public static void main(String[] args) throws IOException, InterruptedException {
        // Initialize and run the optimization
        MultiObjectiveNSGAII.initialization(args);

        if (args.length < 6) {
            System.out.println("Usage: MonteCarloSensitivityRunner <Trip_Item> <Config> <Vertiport_Unit_Candidate> <Scenario_Configuration> <Result_Output> <Sensitivity_Config_File>");
            System.exit(1);
        }
        
        String tripItemFile = args[0];
        String configFile = args[1];
        String vertiportFile = args[2];
        String scenarioFile = args[3];
        String baseOutputDir = args[4];
        String sensitivityConfigFile = args[5];

        // Create directory for this sensitivity analysis
        String sensitivityOutputDir = baseOutputDir + "monte_carlo/";
        MultiObjectiveNSGAII.createFolder(sensitivityOutputDir);
        setFilePaths(args[0], args[1], args[2], args[3], sensitivityOutputDir);
        
        // Load base configuration from file
        SensitivityConfig baseConfig;
        try {
            baseConfig = SensitivityConfig.fromFile(sensitivityConfigFile);
        } catch (IOException e) {
            System.err.println("Failed to load sensitivity config from file: " + e.getMessage());
            System.exit(1);
            return;
        }
        
        System.out.println("Starting Monte Carlo Sensitivity Analysis...");
        System.out.println("Fixed parameters:");
        System.out.println("Output directory: " + sensitivityOutputDir);
        System.out.println("Testing Monte Carlo simulation values:");
        for (int numSims : NUM_SIMULATIONS_VALUES) {
            System.out.println("  " + numSims + " simulations");
        }
        
        // Run experiments in parallel
        runParallelExperiments(tripItemFile, configFile, vertiportFile, scenarioFile, sensitivityOutputDir, baseConfig);
        
        System.out.println("Monte Carlo Sensitivity Analysis completed!");
        System.out.println("Results saved in: " + sensitivityOutputDir);
    }
    
    private static void runParallelExperiments(String tripItemFile, String configFile, 
                                             String vertiportFile, String scenarioFile, String outputDir, SensitivityConfig baseConfig) {
        int numThreads = Math.min(NUM_SIMULATIONS_VALUES.length, Runtime.getRuntime().availableProcessors());
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);
        List<Future<?>> futures = new ArrayList<>();
        
        for (int numSims : NUM_SIMULATIONS_VALUES) {
            Future<?> future = executor.submit(() -> runSingleExperiment(
                tripItemFile, configFile, vertiportFile, scenarioFile, outputDir, numSims, baseConfig));
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
                                          String baseOutputDir, int numSimulations, SensitivityConfig baseConfig) {
        try {
            // Create a temporary config file for this number of simulations
            String tempConfigPath = createTempConfigFile(baseConfig, numSimulations, baseOutputDir);
            
            // Build arguments array for MultiObjectiveNSGAII using SensitivityConfig values
            String[] optimizationArgs = {
                tripItemFile,                    // Trip_Item
                configFile,                      // Config
                vertiportFile,                   // Vertiport_Unit_Candidate
                scenarioFile,                    // Scenario_Configuration
                baseOutputDir,                   // Result_Output
                String.valueOf(baseConfig.getPoolingTimeWindow()),  // BUFFER_END_TIME (convert seconds to minutes)
                String.valueOf(baseConfig.getOriginSearchRadius()), // SEARCH_RADIUS_ORIGIN
                String.valueOf(baseConfig.getDestinationSearchRadius()), // SEARCH_RADIUS_DESTINATION
                String.valueOf(ENABLE_LOCAL_SEARCH),  // ENABLE_LOCAL_SEARCH
                String.valueOf(ENABLE_PRINT_RESULTS), // ENABLE_PRINT_RESULTS
                numSimulations + "/",            // OUTPUT_SUBFOLDER
                tempConfigPath                   // SENSITIVITY_CONFIG (path to temp config file)
            };
            
            System.out.println("Running experiment with " + numSimulations + " Monte Carlo simulations...");

            double[] results = MultiObjectiveNSGAII.callAlgorithm(optimizationArgs);
            
            System.out.println("Completed experiment with " + numSimulations + 
                             " simulations. Best fitness: " + results[3]);
                             
/*            // Clean up temp file
            try {
                java.nio.file.Files.deleteIfExists(java.nio.file.Paths.get(tempConfigPath));
            } catch (IOException e) {
                System.err.println("Warning: Could not delete temporary config file: " + tempConfigPath);
            }*/
            
        } catch (Exception e) {
            System.err.println("Error running experiment with " + numSimulations + 
                             " simulations: " + e.getMessage());
            e.printStackTrace();
        }
    }
    
    private static String createTempConfigFile(SensitivityConfig baseConfig, int numSimulations, String baseOutputDir) throws IOException {
        //String tempDir = System.getProperty("java.io.tmpdir");
        String tempConfigPath = baseOutputDir + "sensitivity_config_mc_" + numSimulations + "_" + System.currentTimeMillis() + ".properties";
        
        try (java.io.FileWriter writer = new java.io.FileWriter(tempConfigPath)) {
            writer.write("# Temporary sensitivity config for Monte Carlo simulations: " + numSimulations + "\n");
            writer.write("numSimulations=" + numSimulations + "\n");
            writer.write("chargingRateKwhPerSecond=" + baseConfig.getChargingRateKwhPerSecond() + "\n");
        }
        
        return tempConfigPath;
    }
}