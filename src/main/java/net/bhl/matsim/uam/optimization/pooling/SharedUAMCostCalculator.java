package net.bhl.matsim.uam.optimization.pooling;

import java.util.HashMap;
import java.util.Map;

/**
 * Utility class for calculating shared UAM costs based on number of co-passengers
 * and redirection time.
 */
public class SharedUAMCostCalculator {

    // Beta factor lookup table based on number of co-passengers
    private static final Map<Integer, Double> BETA_FACTORS = new HashMap<>() {{
        put(1, 0.5);   // 2 total passengers (1 co-passenger)
        put(2, 0.67);  // 3 total passengers (2 co-passengers)
        put(3, 0.75);  // 4 total passengers (3 co-passengers)
        put(4, 0.8);   // 5 total passengers (4 co-passengers)
    }};

    // Alpha factor lookup table based on redirection time in minutes
    private static final Map<Integer, Double> ALPHA_FACTORS = new HashMap<>() {{
        put(5, 0.25);   // 25% discount for 5 minutes redirection
        put(10, 0.50);  // 50% discount for 10 minutes redirection
        put(15, 0.75);  // 75% discount for 15 minutes redirection
    }};

    /**
     * Calculate shared UAM cost based on original UAM price, number of co-passengers,
     * and redirection time.
     *
     * @param originalPrice Original UAM price without sharing
     * @param numCoPassengers Number of co-passengers (excluding the passenger themselves)
     * @param redirectionTimeInSeconds Redirection time in seconds
     * @return Calculated shared UAM price
     */
    public static double calculateSharedUAMCost(double originalPrice, int numCoPassengers, double redirectionTimeInSeconds) {
        // Validate inputs
        if (numCoPassengers < 0 || numCoPassengers > 4) {
            throw new IllegalArgumentException("Number of co-passengers must be between 0 and 4");
        }

        double redirectionTimeMinutes = redirectionTimeInSeconds /  60;
        if (redirectionTimeMinutes < 0) {
            throw new IllegalArgumentException("Redirection time cannot be negative");
        }

        // Get beta factor based on number of co-passengers
        double betaFactor = BETA_FACTORS.get(numCoPassengers);

        // Get alpha factor based on redirection time
        double alphaFactor = getAlphaFactor(redirectionTimeMinutes);

        // Calculate total discount
        double discount = betaFactor * alphaFactor;

        // Calculate final price
        return originalPrice * (1 - discount);
    }

    /**
     * Get alpha factor based on redirection time.
     * Uses linear interpolation for times between defined values.
     *
     * @param redirectionTimeMinutes Redirection time in minutes
     * @return Interpolated alpha factor
     */
    private static double getAlphaFactor(double redirectionTimeMinutes) {
        // Handle cases below minimum or above maximum thresholds //TODO: Handle these cases in a more elegant way?
        if (redirectionTimeMinutes <= 5) {
            return ALPHA_FACTORS.get(5);
        }
        if (redirectionTimeMinutes >= 15) {
            return ALPHA_FACTORS.get(15);
        }

        // Find the surrounding defined time points for interpolation
        int lowerTime = 5 * (int)(redirectionTimeMinutes / 5);
        int upperTime = lowerTime + 5;

        // Get the corresponding alpha factors
        double lowerAlpha = ALPHA_FACTORS.get(lowerTime);
        double upperAlpha = ALPHA_FACTORS.get(upperTime);

        // Perform linear interpolation
        double ratio = (redirectionTimeMinutes - lowerTime) / (upperTime - lowerTime);
        return lowerAlpha + (upperAlpha - lowerAlpha) * ratio;
    }

    //TODO: Calculate total redirection time in minutes based on access and egress time changes.
    //TODO: UAM pricing scheme should also include fixed costs?
}
