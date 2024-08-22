package net.bhl.matsim.uam.optimization.pooling;

public class Optimizer {
    public static void main(String[] args) throws Exception {
        // 1. Initialization
        MultiObjectiveNSGAII.initialization(args);

        // 2. Run the optimization
        GridSearch.main(args);
        BayesianOptimization.main(args);
    }
}
