package net.bhl.matsim.uam.optimization.pooling;

import net.bhl.matsim.uam.optimization.pooling.sn.EVTOLBatteryManager;

/**
 * Configuration class for sensitivity analysis parameters in UAM optimization
 */
public class SensitivityConfig {
    // Monte Carlo simulation parameters
    private int numSimulations;
    
    // Charging parameters
    private double chargingRateKwhPerSecond;
    
    // Fixed pooling parameters (as requested)
    private final double poolingTimeWindow = 3.0 ; // 3 minutes
    private final double originSearchRadius = 6000.0; // meters
    private final double destinationSearchRadius = 6000.0; // meters
    
    // Default constructor with default values
    public SensitivityConfig() {
        this.numSimulations = 1000; // Default value
        this.chargingRateKwhPerSecond = EVTOLBatteryManager.DEFAULT_CHARGING_RATE_KWH_PER_SECOND;
    }
    
    // Constructor with specific parameters
    public SensitivityConfig(int numSimulations, double chargingRateKwhPerSecond) {
        this.numSimulations = numSimulations;
        this.chargingRateKwhPerSecond = chargingRateKwhPerSecond;
    }
    
    // Static factory methods for creating specific configurations
    public static SensitivityConfig forMonteCarloAnalysis(int numSimulations) {
        SensitivityConfig config = new SensitivityConfig();
        config.setNumSimulations(numSimulations);
        return config;
    }
    
    public static SensitivityConfig forChargingRateAnalysis(double chargingRateKwhPerSecond) {
        SensitivityConfig config = new SensitivityConfig();
        config.setChargingRateKwhPerSecond(chargingRateKwhPerSecond);
        return config;
    }
    
    // Getters and setters
    public int getNumSimulations() {
        return numSimulations;
    }
    
    public void setNumSimulations(int numSimulations) {
        this.numSimulations = numSimulations;
    }
    
    public double getChargingRateKwhPerSecond() {
        return chargingRateKwhPerSecond;
    }
    
    public void setChargingRateKwhPerSecond(double chargingRateKwhPerSecond) {
        this.chargingRateKwhPerSecond = chargingRateKwhPerSecond;
    }
    
    public double getPoolingTimeWindow() {
        return poolingTimeWindow;
    }
    
    public double getOriginSearchRadius() {
        return originSearchRadius;
    }
    
    public double getDestinationSearchRadius() {
        return destinationSearchRadius;
    }
    
    // Utility method to create a descriptive string for output directories
    public String getParameterString() {
        return String.format("mc_%d_charging_%.4f", numSimulations, chargingRateKwhPerSecond);
    }
    
    @Override
    public String toString() {
        return String.format("SensitivityConfig{numSimulations=%d, chargingRateKwhPerSecond=%.6f, poolingTimeWindow=%.1f, originSearchRadius=%.1f, destinationSearchRadius=%.1f}", 
                numSimulations, chargingRateKwhPerSecond, poolingTimeWindow, originSearchRadius, destinationSearchRadius);
    }
}