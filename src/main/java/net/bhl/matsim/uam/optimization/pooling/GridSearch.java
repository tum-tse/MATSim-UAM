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

public class GridSearch {
    private static final int numProcessors = Runtime.getRuntime().availableProcessors();
    private static final int bufferDivider = 2;

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
        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        List<Future<double[]>> futures = new ArrayList<>();

        MultiObjectiveNSGAII.initialization(args);
        ThreadCounter threadCounter = new ThreadCounter();

        // Example loop to optimize parameters
        for (double ptw = 5.0; ptw <= 15.0; ptw += 1.0) {
            for (double sro = 500; sro <= 3000; sro += 500) {
                for (double srd = 500; srd <= 3000; srd += 500) {
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
                                    String.valueOf(false),  // ENABLE_LOCAL_SEARCH
                                    String.valueOf(true),  // ENABLE_PRINT_RESULTS
                                    String.valueOf(finalPtw + "_" + finalSro + "_" + finalSrd + "/") // OUTPUT_SUB_DIRECTORY
                            };
                            return MultiObjectiveNSGAII.callAlgorithm(multiObjectiveArgs);
                            } finally {
                                threadCounter.deregister(); // Deregister at the end of the task, even if an exception occurs
                            }
                        }
                    });

                    futures.add(future);
                }
            }
        }

        // Collect results and add to dataset
        for (Future<double[]> future : futures) {
            double[] fitnessScore = future.get();
            // Add the instance to the dataset
            // Note: You'll need to keep track of which parameters correspond to which future
            // This might require additional bookkeeping
            dataset.add(new DenseInstance(1.0, fitnessScore));
        }

        executor.shutdown();

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
