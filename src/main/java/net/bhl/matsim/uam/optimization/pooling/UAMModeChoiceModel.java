package net.bhl.matsim.uam.optimization.pooling;

import java.util.Random;

/**
 * Class to handle UAM mode choice calculations using an ordinal logit model
 * Based on estimated coefficients from biogeme model
 */
public class UAMModeChoiceModel {
    // Model coefficients from estimation
    private static final double B_COPASSENGER = -0.0901;
    private static final double B_COST = -0.000815;
    private static final double B_SHARED_INVEHICLE_TIME = -0.00881;
    private static final double B_SHARED_REDIRECTION_TIME = 0.0316;
    private static final double B_UAM_INVEHICLE_TIME = 0.042;
    private static final double B_UAM_REDIRECTION_TIME = 0.0556;
    private static final double B_WAITING_TIME = -0.0606;

    // Threshold parameters
    private static final double TAU1 = -1.68;
    private static final double TAU1_DIFF_MINUS1 = 1.61;
    private static final double TAU1_DIFF_0 = 0.522;
    private static final double TAU1_DIFF_1 = 1.34;

    private final Random random;

    public UAMModeChoiceModel(Random random) {
        this.random = random;
    }

    /**
     * Calculates utility for shared UAM option
     */
    private double calculateSharedUAMUtility(double cost, double inVehicleTime,
                                             double redirectionTime, double waitingTime, int coPassengers) {
        return B_COST * cost +
                B_SHARED_INVEHICLE_TIME * inVehicleTime +
                B_SHARED_REDIRECTION_TIME * redirectionTime +
                B_WAITING_TIME * waitingTime +
                B_COPASSENGER * coPassengers;
    }

    /**
     * Calculates utility for non-shared UAM option
     */
    private double calculateNonSharedUAMUtility(double cost, double inVehicleTime,
                                                double redirectionTime, double waitingTime) {
        return B_COST * cost +
                B_UAM_INVEHICLE_TIME * inVehicleTime +
                B_UAM_REDIRECTION_TIME * redirectionTime +
                B_WAITING_TIME * waitingTime;
    }

    /**
     * Performs Monte Carlo simulation to determine acceptance probability
     * @return probability of accepting shared UAM
     */
    public double calculateAcceptanceProbability(
            double sharedCost, double sharedInVehicleTime, double sharedRedirectionTime,
            double sharedWaitingTime, int coPassengers,
            double nonSharedCost, double nonSharedInVehicleTime, double nonSharedRedirectionTime,
            double nonSharedWaitingTime,
            int numSimulations) {

        int acceptanceCount = 0;

        for (int i = 0; i < numSimulations; i++) {
            // Add Gumbel-distributed error terms
            double epsilon1 = drawGumbel();
            double epsilon2 = drawGumbel();

            double sharedUtility = calculateSharedUAMUtility(
                    sharedCost, sharedInVehicleTime, sharedRedirectionTime,
                    sharedWaitingTime, coPassengers) + epsilon1;

            double nonSharedUtility = calculateNonSharedUAMUtility(
                    nonSharedCost, nonSharedInVehicleTime, nonSharedRedirectionTime,
                    nonSharedWaitingTime) + epsilon2;

            if (sharedUtility > nonSharedUtility) {
                acceptanceCount++;
            }
        }

        return (double) acceptanceCount / numSimulations;
    }

    /**
     * Draws from standard Gumbel distribution using inverse transform method
     */
    private double drawGumbel() {
        double u = random.nextDouble();
        return -Math.log(-Math.log(u));
    }

    /**
     * Converts 5-category decision to binary choice probability
     * using estimated thresholds
     */
    public double convertToBinaryProbability(double utilityDifference) {
        // Calculate probabilities for each category using ordinal logit formulation
        double p1 = 1.0 / (1.0 + Math.exp(-(TAU1 - utilityDifference)));
        double p2 = 1.0 / (1.0 + Math.exp(-(TAU1 + TAU1_DIFF_MINUS1 - utilityDifference))) - p1;
        double p3 = 1.0 / (1.0 + Math.exp(-(TAU1 + TAU1_DIFF_MINUS1 + TAU1_DIFF_0 - utilityDifference))) - (p1 + p2);
        double p4 = 1.0 / (1.0 + Math.exp(-(TAU1 + TAU1_DIFF_MINUS1 + TAU1_DIFF_0 + TAU1_DIFF_1 - utilityDifference))) - (p1 + p2 + p3);
        double p5 = 1.0 - (p1 + p2 + p3 + p4);

        // Conservative approach: only "probably_shared" (p4) and "definitely_shared" (p5)
        // are considered as acceptance
        return p4 + p5;
    }
}
