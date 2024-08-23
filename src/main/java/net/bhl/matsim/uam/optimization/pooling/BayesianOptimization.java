package net.bhl.matsim.uam.optimization.pooling;

import net.bhl.matsim.uam.analysis.traveltimes.utils.ThreadCounter;
import weka.classifiers.functions.GaussianProcesses;
import weka.core.DenseInstance;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Attribute;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;
import java.util.concurrent.*;
import java.util.List;
import java.util.concurrent.atomic.AtomicReference;
import java.util.logging.Logger;
import java.util.logging.Level;

import static net.bhl.matsim.uam.optimization.pooling.MultiObjectiveNSGAII.setFilePaths;

public class BayesianOptimization {

    private static final int numProcessors = Runtime.getRuntime().availableProcessors();
    private static final int bufferDivider = 1;
    private static final Logger LOGGER = Logger.getLogger(BayesianOptimization.class.getName());
    private static final int TIMEOUT_MINUTES = 6000;

    private GaussianProcesses gaussianProcess;
    private Instances dataset;
    private final long SEED = 4711; // MATSim default Random Seed
    private final Random rand = new Random(SEED);
    private final int NUM_THREADS = Runtime.getRuntime().availableProcessors();
    private final ExecutorService executorService = Executors.newFixedThreadPool(NUM_THREADS);

    // Initialize with empty dataset
    private final ParameterRange poolingTimeWindowRange;
    private final ParameterRange searchRadiusOriginRange;
    private final ParameterRange searchRadiusDestinationRange;

    public BayesianOptimization(ParameterRange poolingTimeWindowRange,
                                ParameterRange searchRadiusOriginRange,
                                ParameterRange searchRadiusDestinationRange) throws Exception {
        this.poolingTimeWindowRange = poolingTimeWindowRange;
        this.searchRadiusOriginRange = searchRadiusOriginRange;
        this.searchRadiusDestinationRange = searchRadiusDestinationRange;

        ArrayList<Attribute> attributes = new ArrayList<>();
        attributes.add(new Attribute("pooling_time_window"));
        attributes.add(new Attribute("search_radius_origin"));
        attributes.add(new Attribute("search_radius_destination"));
        attributes.add(new Attribute("performance_metric"));

        dataset = new Instances("OptimizationDataset", attributes, 0);
        dataset.setClassIndex(3);

        // Initialize GaussianProcesses but don't build it yet
        gaussianProcess = new GaussianProcesses();

        // Add initial data points
        addInitialDataPoints();
    }

    // Add a fixed number of initial data points
    private void addInitialDataPoints() throws Exception {
        int[][] initialPoints = {
                {poolingTimeWindowRange.getMinValue(), searchRadiusOriginRange.getMinValue(), searchRadiusDestinationRange.getMinValue()},
                {poolingTimeWindowRange.getMaxValue(), searchRadiusOriginRange.getMaxValue(), searchRadiusDestinationRange.getMaxValue()},
                {poolingTimeWindowRange.getMidValue(), searchRadiusOriginRange.getMidValue(), searchRadiusDestinationRange.getMidValue()}
        };

        int successfulPoints = 0;
        for (int[] point : initialPoints) {
            try {
                double performanceMetric = evaluateParameters(point[0], point[1], point[2]);
                addDataPoint(point[0], point[1], point[2], performanceMetric);
                successfulPoints++;
            } catch (IllegalStateException e) {
                LOGGER.warning("Failed to evaluate initial point: " + e.getMessage());
            }
        }

        if (successfulPoints == 0) {
            throw new IllegalStateException("Failed to evaluate any initial points. Check parameter ranges and constraints.");
        }
    }

    // Add a new data point to the dataset
    public void addDataPoint(int poolingTimeWindow, int searchRadiusOrigin, int searchRadiusDestination, double performanceMetric) throws Exception {
        double[] values = new double[]{poolingTimeWindow, searchRadiusOrigin, searchRadiusDestination, performanceMetric};
        Instance instance = new DenseInstance(1.0, values);
        dataset.add(instance);

        // Build the classifier after adding a new data point
        gaussianProcess.buildClassifier(dataset);
    }

    // Predict performance for a given set of parameters
    public double predictPerformance(int poolingTimeWindow, int searchRadiusOrigin, int searchRadiusDestination) throws Exception {
        if (dataset.numInstances() == 0) {
            throw new IllegalStateException("Cannot predict performance without any training data.");
        }
        double[] values = new double[]{poolingTimeWindow, searchRadiusOrigin, searchRadiusDestination, Double.NaN};
        Instance instance = new DenseInstance(1.0, values);
        instance.setDataset(dataset);
        return gaussianProcess.classifyInstance(instance);
    }

    // Deterministic acquisition function (Upper Confidence Bound)
    private double acquisitionFunction(int poolingTimeWindow, int searchRadiusOrigin, int searchRadiusDestination) throws Exception {
        double mean = predictPerformance(poolingTimeWindow, searchRadiusOrigin, searchRadiusDestination);
        double std = Math.sqrt(gaussianProcess.getStandardDeviation(new DenseInstance(1.0, new double[]{poolingTimeWindow, searchRadiusOrigin, searchRadiusDestination, Double.NaN})));
        return mean + 2 * std;  // UCB with kappa = 2
    }

    // Parallelized optimization method
    public int[] optimizeParameters(int iterations) throws Exception {
        AtomicReference<Double> bestPerformance = new AtomicReference<>(Double.NEGATIVE_INFINITY);
        AtomicReference<int[]> bestParams = new AtomicReference<>(new int[3]);

        ThreadCounter threadCounter = new ThreadCounter();

        for (int i = 0; i < iterations; i++) {
            List<Future<OptimizationResult>> futures = new ArrayList<>();

            for (int ptw = poolingTimeWindowRange.getMinValue(); ptw <= poolingTimeWindowRange.getMaxValue(); ptw++) {
                for (int sro = searchRadiusOriginRange.getMinValue(); sro <= searchRadiusOriginRange.getMaxValue(); sro += 100) {
                    for (int srd = searchRadiusDestinationRange.getMinValue(); srd <= searchRadiusDestinationRange.getMaxValue(); srd += 100) {
                        final int finalPtw = ptw;
                        final int finalSro = sro;
                        final int finalSrd = srd;

                        while (threadCounter.getProcesses() >= numProcessors/bufferDivider - 1)
                            Thread.sleep(200);

                        futures.add(executorService.submit(new Callable<OptimizationResult>() {
                            @Override
                            public OptimizationResult call() throws Exception {
                                threadCounter.register(); // Register at the start of the task
                                try {
                                    double acquisitionValue = acquisitionFunction(finalPtw, finalSro, finalSrd);
                                    return new OptimizationResult(finalPtw, finalSro, finalSrd, acquisitionValue);
                                } catch (Exception e) {
                                    LOGGER.warning("Failed to evaluate point (" + finalPtw + ", " + finalSro + ", " + finalSrd + "): " + e.getMessage());
                                    return null;
                                } finally {
                                threadCounter.deregister(); // Deregister at the end of the task, even if an exception occurs
                                }
                            }
                        }));
                    }
                }
            }

            OptimizationResult bestResult = null;
            double bestAcquisitionValue = Double.NEGATIVE_INFINITY;

            for (Future<OptimizationResult> future : futures) {
                OptimizationResult result = future.get();
                if (result != null && result.performance > bestAcquisitionValue) {
                    bestResult = result;
                    bestAcquisitionValue = result.performance;
                }
            }

            if (bestResult != null) {
                try {
                    double actualPerformance = evaluateParameters(bestResult.poolingTimeWindow, bestResult.searchRadiusOrigin, bestResult.searchRadiusDestination);
                    addDataPoint(bestResult.poolingTimeWindow, bestResult.searchRadiusOrigin, bestResult.searchRadiusDestination, actualPerformance);

                    if (actualPerformance > bestPerformance.get()) {
                        bestPerformance.set(actualPerformance);
                        bestParams.set(new int[]{bestResult.poolingTimeWindow, bestResult.searchRadiusOrigin, bestResult.searchRadiusDestination});
                    }
                    LOGGER.info("Iteration " + i + ": Best performance so far: " + bestPerformance.get());
                } catch (IllegalStateException e) {
                    LOGGER.warning("Failed to evaluate best point in iteration " + i + ": " + e.getMessage());
                }
            } else {
                LOGGER.warning("No valid point found in iteration " + i);
            }
        }

        executorService.shutdown();
        if (!executorService.awaitTermination(TIMEOUT_MINUTES, TimeUnit.MINUTES)) {
            LOGGER.warning("Timeout occurred. Not all tasks completed.");
        }

        return bestParams.get();
    }

    private double evaluateParameters(int poolingTimeWindow, int searchRadiusOrigin, int searchRadiusDestination) throws IOException, InterruptedException {
        String[] multiObjectiveArgs = new String[]{
                String.valueOf("" ), // INPUT_FILE
                String.valueOf("" ), // INPUT_FILE
                String.valueOf("" ), // INPUT_FILE
                String.valueOf("" ), // INPUT_FILE
                String.valueOf("" ), // OUTPUT_DIRECTORY
                String.valueOf(poolingTimeWindow),
                String.valueOf(searchRadiusOrigin),
                String.valueOf(searchRadiusDestination),
                String.valueOf(false),
                String.valueOf(true),
                String.valueOf(poolingTimeWindow + "_" + searchRadiusOrigin + "_" + searchRadiusDestination + "/") // OUTPUT_SUB_DIRECTORY
        };
        double[] actualPerformance = MultiObjectiveNSGAII.callAlgorithm(multiObjectiveArgs);
        return actualPerformance[3];
    }

    private static class OptimizationResult {
        final int poolingTimeWindow;
        final int searchRadiusOrigin;
        final int searchRadiusDestination;
        final double performance;

        OptimizationResult(int ptw, int sro, int srd, double perf) {
            poolingTimeWindow = ptw;
            searchRadiusOrigin = sro;
            searchRadiusDestination = srd;
            performance = perf;
        }
    }

    public static class ParameterRange {
        private final int min;
        private final int max;

        public ParameterRange(int min, int max) {
            this.min = min;
            this.max = max;
        }

        public int getMinValue() {
            return min;
        }

        public int getMaxValue() {
            return max;
        }

        public int getMidValue() {
            return min + (max - min) / 2;
        }
    }

    public static void main(String[] args) {
        try {
            ParameterRange poolingTimeWindowRange = new ParameterRange(1, 25);
            ParameterRange searchRadiusOriginRange = new ParameterRange(1000, 50000);
            ParameterRange searchRadiusDestinationRange = new ParameterRange(1000, 50000);

            //MultiObjectiveNSGAII.initialization(args);
            String outputSubFolder = args[4].endsWith("/") ? args[4] + "bayesian_optimization" : args[4] + "/" + "bayesian_optimization";
            MultiObjectiveNSGAII.createFolder(outputSubFolder);
            setFilePaths(args[0], args[1], args[2], args[3], outputSubFolder);

            BayesianOptimization optimization = new BayesianOptimization(
                    poolingTimeWindowRange, searchRadiusOriginRange, searchRadiusDestinationRange);

            int[] bestParams = optimization.optimizeParameters(20);
            System.out.println("Best Parameters: Pooling Time Window = " + bestParams[0] +
                    ", Search Radius Origin = " + bestParams[1] +
                    ", Search Radius Destination = " + bestParams[2]);
        } catch (Exception e) {
            LOGGER.log(Level.SEVERE, "An error occurred during optimization", e);
        }
    }
}
