package net.bhl.matsim.uam.optimization.pooling;

import java.util.Random;

/**
 * Class to handle UAM mode choice calculations using an ordinal logit model
 * Based on estimated coefficients from biogeme model
 */
public class UAMModeChoiceModel {
    // Model coefficients from estimation
    // UAM coefficients
    private static final double B_ASC_UAM = -4.66;
    private static final double B_UAM_COST = -0.0217;
    private static final double B_UAM_INVEHICLE_TIME = -0.153;

    private static final double B_OUT_OF_VEHICLE_TIME = -0.0813;

    // shared UAM coefficients
    private static final double B_SHARED_COST = -0.025;
    private static final double B_SHARED_INVEHICLE_TIME = -0.177;
    private static final double B_COPASSENGER = -0.571;


    private final Random random;
    private final int numberSimulations;

    public UAMModeChoiceModel(Random random, int numberSimulations) {
        this.random = random;
        this.numberSimulations = numberSimulations;
    }

    /**
     * Calculates utility for shared UAM option
     */
    private double calculateSharedUAMUtility(double cost, double inVehicleTime,
                                             double outOfVehicleTime, int coPassengers) {
        return B_SHARED_COST * cost +
                B_SHARED_INVEHICLE_TIME * inVehicleTime/60 +
                B_OUT_OF_VEHICLE_TIME * outOfVehicleTime/60 +
                B_COPASSENGER * coPassengers;
    }

    /**
     * Calculates utility for non-shared UAM option
     */
    private double calculateNonSharedUAMUtility(double cost, double inVehicleTime,
                                                double outOfVehicleTime) {
        return B_UAM_COST * cost +
                B_UAM_INVEHICLE_TIME * inVehicleTime/60 +
                B_OUT_OF_VEHICLE_TIME * outOfVehicleTime/60 +
                B_ASC_UAM;
    }

    /**
     * Performs Monte Carlo simulation to determine acceptance probability
     * @return probability of accepting shared UAM
     */
    public double calculateAcceptanceProbability(
            double sharedCost, double sharedInVehicleTime, double sharedRedirectionTime,
            double sharedWaitingTime, int coPassengers,
            double nonSharedCost, double nonSharedInVehicleTime, double nonSharedRedirectionTime) {

        int acceptanceCount = 0;

        for (int i = 0; i < numberSimulations; i++) {
            // Add Gumbel-distributed error terms
            double epsilon1 = drawGumbel();
            double epsilon2 = drawGumbel();

            double sharedUtility = calculateSharedUAMUtility(
                    sharedCost, sharedInVehicleTime, sharedRedirectionTime+sharedWaitingTime, coPassengers) + epsilon1;

            double nonSharedUtility = calculateNonSharedUAMUtility(
                    nonSharedCost, nonSharedInVehicleTime, nonSharedRedirectionTime) + epsilon2;

            if (sharedUtility > nonSharedUtility) {
                acceptanceCount++;
            }
        }

        return (double) acceptanceCount / numberSimulations;
    }

    /**
     * Draws from standard Gumbel distribution using inverse transform method
     */
    private double drawGumbel() {
        double u = random.nextDouble();
        return -Math.log(-Math.log(u));
    }
}
