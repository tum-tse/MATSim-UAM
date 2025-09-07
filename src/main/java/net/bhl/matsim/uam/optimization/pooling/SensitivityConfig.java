package net.bhl.matsim.uam.optimization.pooling;

import net.bhl.matsim.uam.optimization.pooling.sn.EVTOLBatteryManager;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Properties;

/**
 * Configuration class for sensitivity analysis parameters in UAM optimization
 * Reads configuration from a properties file.
 */
public class SensitivityConfig {
    // Monte Carlo simulation parameters
    private int numSimulations;
    
    // Charging parameters
    private double chargingRateKwhPerSecond;
    private boolean vehicleReuseStrategy;
    
    // Fixed pooling parameters (as requested)
    private double bufferStartTime;
    private double poolingTimeWindow; // 3 minutes
    private double originSearchRadius; // meters
    private double destinationSearchRadius; // meters
    
    // Constructor that reads from config file path
    public SensitivityConfig(String configFilePath) throws IOException {
        loadFromFile(configFilePath);
    }
    
    // Private method to load configuration from file
    private void loadFromFile(String configFilePath) throws IOException {
        // Check if config file exists
        if (!Files.exists(Paths.get(configFilePath))) {
            throw new IOException("Configuration file not found: " + configFilePath);
        }
        
        // Try to load from file and override defaults if successful
        Properties props = new Properties();
        
        try (FileInputStream fis = new FileInputStream(configFilePath)) {
            props.load(fis);
            
            // Override defaults with values from file (using Properties.getProperty with defaults)
            this.bufferStartTime = Double.parseDouble(
                    props.getProperty("bufferStartTime", "25200.0")
            );
            this.poolingTimeWindow = Double.parseDouble(
                    props.getProperty("poolingTimeWindow", "3.0")
            );
            this.originSearchRadius = Double.parseDouble(
                    props.getProperty("originSearchRadius", "6000.0")
            );
            this.destinationSearchRadius = Double.parseDouble(
                    props.getProperty("destinationSearchRadius", "6000.0")
            );

            this.numSimulations = Integer.parseInt(
                props.getProperty("numSimulations", "1000")
            );
            
            this.chargingRateKwhPerSecond = Double.parseDouble(
                props.getProperty("chargingRateKwhPerSecond", 
                    String.valueOf(EVTOLBatteryManager.DEFAULT_CHARGING_RATE_KWH_PER_SECOND))
            );

            this.vehicleReuseStrategy = Boolean.parseBoolean(
                    props.getProperty("vehicleReuseStrategy", "false")
            );

        }
    }
    
    // Static factory methods for creating specific configurations from file
    public static SensitivityConfig fromFile(String configFilePath) throws IOException {
        return new SensitivityConfig(configFilePath);
    }
    
    // Getters and setters
    public int getNumSimulations() {
        return numSimulations;
    }
    
    public double getChargingRateKwhPerSecond() {
        return chargingRateKwhPerSecond;
    }

    public boolean getVehicleReuseStrategy() {return vehicleReuseStrategy;}

    public double getBufferStartTime() {return bufferStartTime;}
    
    public double getPoolingTimeWindow() {
        return poolingTimeWindow;
    }
    
    public double getOriginSearchRadius() {
        return originSearchRadius;
    }
    
    public double getDestinationSearchRadius() {
        return destinationSearchRadius;
    }
    
    @Override
    public String toString() {
        return String.format("SensitivityConfig{numSimulations=%d, chargingRateKwhPerSecond=%.6f, vehicleReuseStrategy=%b, bufferStartTime=%.1f, poolingTimeWindow=%.1f, originSearchRadius=%.1f, destinationSearchRadius=%.1f}",
                numSimulations, chargingRateKwhPerSecond, vehicleReuseStrategy, bufferStartTime, poolingTimeWindow, originSearchRadius, destinationSearchRadius);
    }
}