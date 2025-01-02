package net.bhl.matsim.uam.optimization.pooling;

import java.util.HashMap;
import java.util.Map;

/**
 * Utility class for calculating shared UAM costs based on number of co-passengers
 * and redirection time.
 */
public class SharedUAMCostCalculator {

    /**
     * Calculate beta factor based on number of co-passengers
     * @param numCoPassengers Number of co-passengers
     * @return Beta factor for cost reduction
     */
    private static double calculateBetaFactor(int numCoPassengers) {
        return (double) numCoPassengers / (numCoPassengers + 1);
    }

    /**
     * Calculate alpha factor based on redirection time
     * @param redirectionTimeMinutes Redirection time in minutes
     * @return Alpha factor for cost reduction
     */
    private static double calculateAlphaFactor(double redirectionTimeMinutes) {
        // Handle cases below minimum or above maximum thresholds //TODO: Handle these cases in a more elegant way?
        if (redirectionTimeMinutes <= 5) {
            return 0.25;
        } else if (redirectionTimeMinutes >= 15) {
            return 0.75;
        } else {
            // Linear interpolation between 5 and 15 minutes
            return 0.25 + 0.05 * (redirectionTimeMinutes - 5);
        }
    }

    /**
     * Calculate shared UAM cost using mathematical model based on original UAM price, number of co-passengers,
     * and redirection time.
     *
     * @param originalPrice Original UAM price without sharing
     * @param numCoPassengers Number of co-passengers (excluding the passenger themselves)
     * @param redirectionTimeInSeconds Redirection time in seconds
     * @return Calculated shared UAM price
     */
    public static double calculateSharedUAMCost(double originalPrice,
                                                int numCoPassengers,
                                                double redirectionTimeInSeconds) {
        // Validate inputs
        if (numCoPassengers <= 0) {
            throw new IllegalArgumentException("Number of co-passengers must be larger than 0");
        }

        double redirectionTimeMinutes = redirectionTimeInSeconds / 60;
        if (redirectionTimeMinutes < 0) {
            throw new IllegalArgumentException("Redirection time cannot be negative");
        }

        double beta = calculateBetaFactor(numCoPassengers);
        double alpha = calculateAlphaFactor(redirectionTimeMinutes);

        return originalPrice * (1 - beta * alpha);
    }
    //TODO: Calculate total redirection time in minutes based on access and egress time changes.
    //TODO: UAM pricing scheme should also include fixed costs?
}
