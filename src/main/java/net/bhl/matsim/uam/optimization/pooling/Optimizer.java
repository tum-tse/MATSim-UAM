package net.bhl.matsim.uam.optimization.pooling;

import java.io.File;

public class Optimizer {
    public static void main(String[] args) throws Exception {
        // Check if output directory already has content
        String outputDir = args[4];
        if (!outputDir.endsWith("/")) {
            outputDir += "/";
        }

        File directory = new File(outputDir);
        if (directory.exists() && directory.isDirectory()) {
            File outputClusteredCandidatesFile = new File(directory, "output_clustered_candidates.csv");
            if (outputClusteredCandidatesFile.exists()) {
                System.out.println("Output directory " + outputDir + " already contains 'output_clustered_candidates.csv'. Terminating to avoid overwriting data.");
                return; // Exit the program
            }
        } else {
            System.out.println("Output directory " + outputDir + " does not exist. Creating a new one.");
        }

        // 1. Initialization
        MultiObjectiveNSGAII.initialization(args);

        // 2. Run the optimization
        GridSearch.main(args);
        //BayesianOptimization.main(args);
    }
}
