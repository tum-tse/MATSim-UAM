package net.bhl.matsim.uam.optimization.pooling;

import net.bhl.matsim.uam.analysis.traveltimes.utils.ThreadCounter;
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.Attribute;

import java.util.ArrayList;
import java.util.concurrent.*;
import java.util.List;
import java.util.logging.*;

import static net.bhl.matsim.uam.optimization.pooling.MultiObjectiveNSGAII.setFilePaths;

public class GridSearch {
    private static final int numProcessors = Runtime.getRuntime().availableProcessors();
    private static final int bufferDivider = 1;
    private static final Logger logger = Logger.getLogger(GridSearch.class.getName());
    public static final int TIMEOUT_MINUTES = 60*24*10;

    public static void main(String[] args) throws Exception {
        // Configure logger
        Logger logger = Logger.getLogger(GridSearch.class.getName());
        // Remove any existing handlers
        for (Handler handler : logger.getHandlers()) {
            logger.removeHandler(handler);
        }
        // Prevent parent loggers from also logging the messages
        logger.setUseParentHandlers(false);
        // Add our handler
        ConsoleHandler handler = new ConsoleHandler();
        handler.setLevel(Level.INFO);
        logger.setLevel(Level.INFO);
        logger.addHandler(handler);

        // Define the attributes
        ArrayList<Attribute> attributes = new ArrayList<>();
        attributes.add(new Attribute("poolingTimeWindow"));
        attributes.add(new Attribute("searchRadiusOrigin"));
        attributes.add(new Attribute("searchRadiusDestination"));
        attributes.add(new Attribute("FitnessScore"));

        // Create the dataset with these attributes
        Instances dataset = new Instances("OptimizationData", attributes, 0);
        dataset.setClassIndex(dataset.numAttributes() - 1);

        // Set up the output directory
        String outputSubFolder = args[4].endsWith("/") ? args[4] + "grid_search" : args[4] + "/" + "grid_search";
        MultiObjectiveNSGAII.createFolder(outputSubFolder);
        setFilePaths(args[0], args[1], args[2], args[3], outputSubFolder);

        // Create a thread pool with fixed number of threads
        int numThreads = Runtime.getRuntime().availableProcessors() / bufferDivider;
        ThreadPoolExecutor executor = (ThreadPoolExecutor) Executors.newFixedThreadPool(numThreads);

        // Use a completed tasks counter to track progress
        final CountDownLatch completionLatch = new CountDownLatch(5 * 10 * 10); // Total tasks
        ThreadCounter threadCounter = new ThreadCounter();

        try {
            logger.info("Starting grid search with parameters:");
            logger.info("Pooling time window range: 1.0 to 5.0");
            logger.info("Search radius origin range: 1000 to 10000");
            logger.info("Search radius destination range: 1000 to 10000");
            logger.info("Total combinations to evaluate: " + (5 * 10 * 10));
            logger.info("Using " + numThreads + " threads");

            // Store parameter sets and results for tracking
            List<double[]> parameterSets = new ArrayList<>();
            List<Future<double[]>> futures = new ArrayList<>();

            // Submit all tasks
            for (double ptw = 1; ptw <= 5.0; ptw += 1.0) {
                for (double sro = 1000; sro <= 10000; sro += 1000) {
                    for (double srd = 1000; srd <= 10000; srd += 1000) {
                        final double finalPtw = ptw;
                        final double finalSro = sro;
                        final double finalSrd = srd;

                        // Store parameter set for reference
                        parameterSets.add(new double[]{finalPtw, finalSro, finalSrd});

                        // Wait until thread is available in the pool
                        while (threadCounter.getProcesses() >= numProcessors/bufferDivider - 1)
                            Thread.sleep(200);

                        logger.info("Submitting task: PTW=" + finalPtw + ", SRO=" + finalSro + ", SRD=" + finalSrd);

                        // Submit task to thread pool
                        Future<double[]> future = executor.submit(new Callable<double[]>() {
                            @Override
                            public double[] call() throws Exception {
                                threadCounter.register(); // Register at the start of the task
                                try {
                                    logger.info("Starting execution: PTW=" + finalPtw + ", SRO=" + finalSro + ", SRD=" + finalSrd);

                                    String[] multiObjectiveArgs = {
                                            String.valueOf("" ), // INPUT_FILE
                                            String.valueOf("" ), // INPUT_FILE
                                            String.valueOf("" ), // INPUT_FILE
                                            String.valueOf("" ), // INPUT_FILE
                                            String.valueOf("" ), // OUTPUT_DIRECTORY
                                            String.valueOf(finalPtw), // BUFFER_END_TIME
                                            String.valueOf(finalSro), // SEARCH_RADIUS_ORIGIN
                                            String.valueOf(finalSrd), // SEARCH_RADIUS_DESTINATION
                                            String.valueOf(true),  // ENABLE_LOCAL_SEARCH
                                            String.valueOf(true),  // ENABLE_PRINT_RESULTS
                                            String.valueOf(finalPtw + "_" + finalSro + "_" + finalSrd + "/") // OUTPUT_SUB_DIRECTORY
                                    };
                                    double[] result = MultiObjectiveNSGAII.callAlgorithm(multiObjectiveArgs);
                                    logger.info("Completed execution: PTW=" + finalPtw + ", SRO=" + finalSro + ", SRD=" + finalSrd);
                                    return result;
                                } catch (Exception e) {
                                    logger.log(Level.SEVERE, "Task failed for ptw=" + finalPtw + ", sro=" + finalSro + ", srd=" + finalSrd, e);
                                    return null;
                                } finally {
                                    threadCounter.deregister(); // Deregister at the end of the task, even if an exception occurs
                                    completionLatch.countDown(); // Decrement completion counter
                                    logger.info("Remaining tasks: " + completionLatch.getCount());
                                }
                            }
                        });

                        futures.add(future);
                    }
                }
            }

            // Signal no more tasks will be submitted
            logger.info("All tasks submitted. No more tasks will be submitted. Waiting for completion...");
            executor.shutdown();

            // Wait for tasks to complete with a timeout
            logger.info("Waiting for all tasks to complete...");
            boolean completed = completionLatch.await(TIMEOUT_MINUTES, TimeUnit.MINUTES);

            if (!completed) {
                logger.warning("Timeout reached. Not all tasks have completed.");

                // Print information about incomplete tasks
                int activeCount = executor.getActiveCount();
                int queueSize = executor.getQueue().size();
                logger.info("Active threads: " + activeCount);
                logger.info("Tasks still in queue: " + queueSize);

                // Continue waiting until all tasks complete
                logger.info("Continuing to wait until all tasks finish...");

                boolean allDone = false;
                while (!allDone) {
                    allDone = true;
                    int remainingActive = 0;

                    for (Future<double[]> future : futures) {
                        if (!future.isDone()) {
                            allDone = false;
                            remainingActive++;
                        }
                    }

                    if (!allDone) {
                        logger.info("Still waiting for " + remainingActive + " tasks to complete...");
                        Thread.sleep(600000); // Check every 10 minute
                    }
                }
            }

            logger.info("All tasks have completed. Processing results...");

            // Process results
            int completedTasks = 0;
            int failedTasks = 0;
            double bestScore = Double.NEGATIVE_INFINITY;
            double[] bestParams = null;

            for (int i = 0; i < futures.size(); i++) {
                Future<double[]> future = futures.get(i);
                double[] params = parameterSets.get(i);

                try {
                    double[] result = future.get(); // Should return immediately if done
                    if (result != null) {
                        double score = result[3]; // Assume fitness score is at index 3

                        // Track best score
                        if (score > bestScore) {
                            bestScore = score;
                            bestParams = new double[]{params[0], params[1], params[2], score};
                        }

                        completedTasks++;
                        logger.info(String.format("Result for PTW=%.1f, SRO=%.1f, SRD=%.1f: Score=%.4f",
                                params[0], params[1], params[2], score));
                    } else {
                        failedTasks++;
                        logger.warning(String.format("Null result for PTW=%.1f, SRO=%.1f, SRD=%.1f",
                                params[0], params[1], params[2]));
                    }
                } catch (Exception e) {
                    failedTasks++;
                    logger.log(Level.WARNING, String.format("Error getting result for PTW=%.1f, SRO=%.1f, SRD=%.1f",
                            params[0], params[1], params[2]), e);
                }
            }

            // Report results
            logger.info("Grid search completed.");
            logger.info("Total tasks: " + futures.size());
            logger.info("Completed tasks: " + completedTasks);
            logger.info("Failed tasks: " + failedTasks);

            if (bestParams != null) {
                logger.info("Best configuration found:");
                logger.info(String.format("Pooling Time Window: %.1f", bestParams[0]));
                logger.info(String.format("Search Radius Origin: %.1f", bestParams[1]));
                logger.info(String.format("Search Radius Destination: %.1f", bestParams[2]));
                logger.info(String.format("Fitness Score: %.4f", bestParams[3]));
            } else {
                logger.warning("No valid results obtained. Cannot determine best parameters.");
            }

        } catch (Exception e) {
            logger.log(Level.SEVERE, "An error occurred in the main execution", e);
        } finally {
            // Final check to ensure all threads are terminated
            if (!executor.isTerminated()) {
                logger.warning("Forcing termination of executor service.");
                List<Runnable> droppedTasks = executor.shutdownNow();
                logger.info(droppedTasks.size() + " tasks were never executed.");
            }
        }

        logger.info("GridSearch execution completed.");
    }
}