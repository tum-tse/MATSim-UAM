package net.bhl.matsim.uam.optimization.pooling;

import net.bhl.matsim.uam.analysis.traveltimes.utils.ThreadCounter;
import weka.classifiers.Evaluation;
import weka.classifiers.functions.SimpleLinearRegression;
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
    private static final int TIMEOUT_MINUTES = 6000;

    public static void main(String[] args) throws Exception {
        // Define the attributes
        ArrayList<Attribute> attributes = new ArrayList<>();
        attributes.add(new Attribute("poolingTimeWindow"));
        attributes.add(new Attribute("searchRadiusOrigin"));
        attributes.add(new Attribute("searchRadiusDestination"));
        attributes.add(new Attribute("FitnessScore"));

        // Create the dataset with these attributes
        Instances dataset = new Instances("OptimizationData", attributes, 0);
        dataset.setClassIndex(dataset.numAttributes() - 1);

        // Create a thread pool
        int numThreads = Runtime.getRuntime().availableProcessors() / bufferDivider;
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        List<Future<double[]>> futures = new ArrayList<>();

        //MultiObjectiveNSGAII.initialization(args);
        String outputSubFolder = args[4].endsWith("/") ? args[4] + "grid_search" : args[4] + "/" + "grid_search";
        MultiObjectiveNSGAII.createFolder(outputSubFolder);
        setFilePaths(args[0], args[1], args[2], args[3], outputSubFolder);
        ThreadCounter threadCounter = new ThreadCounter();

        try {
            // Example loop to optimize parameters
            for (double ptw = 1; ptw <= 5.0; ptw += 1.0) {
                for (double sro = 1000; sro <= 10000; sro += 1000) {
                    for (double srd = 1000; srd <= 10000; srd += 1000) {
                        final double finalPtw = ptw;
                        final double finalSro = sro;
                        final double finalSrd = srd;

                        while (threadCounter.getProcesses() >= numProcessors/bufferDivider - 1)
                            Thread.sleep(200);
                        // Submit task to thread pool
                        Future<double[]> future = executor.submit(new Callable<>() {
                            @Override
                            public double[] call() throws Exception {
                                threadCounter.register(); // Register at the start of the task
                                try {
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
                                    return MultiObjectiveNSGAII.callAlgorithm(multiObjectiveArgs);
                                } catch (Exception e) {
                                    logger.log(Level.SEVERE, "Task failed for ptw=" + finalPtw + ", sro=" + finalSro + ", srd=" + finalSrd, e);
                                    return null;
                                } finally {
                                    threadCounter.deregister(); // Deregister at the end of the task, even if an exception occurs
                                }
                            }
                        });

                        futures.add(future);
                    }
                }
            }

            executor.shutdown();
            if (!executor.awaitTermination(TIMEOUT_MINUTES, TimeUnit.MINUTES)) {
                logger.warning("Timeout occurred. Not all tasks completed.");
            }

            // Collect results and add to dataset
            for (Future<double[]> future : futures) {
                try {
                    double[] fitnessScore = future.get(6000, TimeUnit.MINUTES);
                    if (fitnessScore != null) {
                        // Add the instance to the dataset
                        // Note: You'll need to keep track of which parameters correspond to which future
                        // This might require additional bookkeeping
                        dataset.add(new DenseInstance(1.0, fitnessScore));
                    }
                } catch (Exception e) {
                    logger.log(Level.WARNING, "Error getting task result", e);
                }
            }
        } catch (Exception e) {
            logger.log(Level.SEVERE, "An error occurred in the main execution", e);
        } finally {
            executor.shutdownNow();
        }

        // Use Weka to find the best parameters using a simple linear regression
        //SimpleLinearRegression slr = new SimpleLinearRegression();
        //slr.buildClassifier(dataset);

        // Evaluate the model
        //Evaluation eval = new Evaluation(dataset);
        //eval.evaluateModel(slr, dataset);

        // Output the results
        //System.out.println("Best fitness score: " + eval.meanAbsoluteError());
        //System.out.println("Summary: " + eval.toSummaryString());
    }
}