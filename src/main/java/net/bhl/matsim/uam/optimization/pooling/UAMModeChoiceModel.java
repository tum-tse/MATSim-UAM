package net.bhl.matsim.uam.optimization.pooling;

import java.util.*;

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

    // Store the most recent simulation results
    private Map<String, Map<Integer, Boolean>> lastSimulationChoices; // Map<TripId, Map<ScenarioNumber, Choice>>
    private Map<String, Double> lastAcceptanceProbabilities; // Map<TripId, Probability>

    public UAMModeChoiceModel(Random random, int numberSimulations) {
        this.random = random;
        this.numberSimulations = numberSimulations;
        this.lastSimulationChoices = new HashMap<>();
        this.lastAcceptanceProbabilities = new HashMap<>();
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
     * Class to hold simulation results
     */
    public static class SimulationResult {
        private final Map<Integer, Boolean> choices; // Map of scenario number to choice
        private final double acceptanceProbability;
        private final String tripId;

        public SimulationResult(String tripId, Map<Integer, Boolean> choices, double acceptanceProbability) {
            this.tripId = tripId;
            this.choices = choices;
            this.acceptanceProbability = acceptanceProbability;
        }

        public Map<Integer, Boolean> getChoices() {
            return choices;
        }

        public double getAcceptanceProbability() {
            return acceptanceProbability;
        }

        public String getTripId() {
            return tripId;
        }
    }

    /**
     * Performs Monte Carlo simulation to determine acceptance probability
     * @return SimulationResult containing individual choices and overall probability of accepting shared UAM
     */
    public SimulationResult simulateChoices(
            String tripId,
            double sharedCost, double sharedInVehicleTime, double sharedRedirectionTime,
            double sharedWaitingTime, int coPassengers,
            double nonSharedCost, double nonSharedInVehicleTime, double nonSharedRedirectionTime) {

        int acceptanceCount = 0;
        Map<Integer, Boolean> choices = new HashMap<>();

        for (int i = 0; i < numberSimulations; i++) {
            // Add Gumbel-distributed error terms
            double epsilon1 = drawGumbel();
            double epsilon2 = drawGumbel();

            double sharedUtility = calculateSharedUAMUtility(
                    sharedCost, sharedInVehicleTime, sharedRedirectionTime+sharedWaitingTime, coPassengers) + epsilon1;

            double nonSharedUtility = calculateNonSharedUAMUtility(
                    nonSharedCost, nonSharedInVehicleTime, nonSharedRedirectionTime) + epsilon2;

            boolean accepted = sharedUtility > nonSharedUtility;
            choices.put(i, accepted);
            if (accepted) {
                acceptanceCount++;
            }
        }

        double probability = (double) acceptanceCount / numberSimulations;
        this.lastSimulationChoices.put(tripId, choices);
        this.lastAcceptanceProbabilities.put(tripId, probability);

        return new SimulationResult(tripId, choices, probability);
    }

    /**
     * Legacy method for backward compatibility
     * @return probability of accepting shared UAM
     */
    public double calculateAcceptanceProbability(
            String tripId,
            double sharedCost, double sharedInVehicleTime, double sharedRedirectionTime,
            double sharedWaitingTime, int coPassengers,
            double nonSharedCost, double nonSharedInVehicleTime, double nonSharedRedirectionTime) {

        SimulationResult result = simulateChoices(
                tripId,
                sharedCost, sharedInVehicleTime, sharedRedirectionTime,
                sharedWaitingTime, coPassengers,
                nonSharedCost, nonSharedInVehicleTime, nonSharedRedirectionTime);

        return result.getAcceptanceProbability();
    }

    /**
     * Get choices for a specific trip
     * @param tripId The ID of the trip to get choices for
     * @return Map of scenario numbers to choices for the specified trip
     */
    public Map<Integer, Boolean> getChoicesForTrip(String tripId) {
        return new HashMap<>(lastSimulationChoices.getOrDefault(tripId, null));
    }

    /**
     * Get acceptance probability for a specific trip
     * @param tripId The ID of the trip to get probability for
     * @return Acceptance probability for the specified trip
     */
    public double getAcceptanceProbabilityForTrip(String tripId) {
        return lastAcceptanceProbabilities.getOrDefault(tripId, null);
    }

    /**
     * Get all simulation choices for all trips
     * @return Map of trip IDs to their scenario choices
     */
    public Map<String, Map<Integer, Boolean>> getAllSimulationChoices() {
        return new HashMap<>(lastSimulationChoices);
    }

    /**
     * Get all acceptance probabilities
     * @return Map of trip IDs to their acceptance probabilities
     */
    public Map<String, Double> getAllAcceptanceProbabilities() {
        return new HashMap<>(lastAcceptanceProbabilities);
    }

    /**
     * Draws from standard Gumbel distribution using inverse transform method
     */
    private double drawGumbel() {
        double u = random.nextDouble();
        return -Math.log(-Math.log(u));
    }
}