package net.bhl.matsim.uam.optimization.pooling.sn;

import org.matsim.api.core.v01.Coord;

/**
 * Battery management for eVTOL vehicles
 * Handles charging, discharging, and energy consumption calculations
 */
public class EVTOLBatteryManager {

    // Constants from requirements (converted to proper units)
    public static final double BATTERY_CAPACITY_KWH = 125.0; // kWh
    public static final double DEFAULT_CHARGING_RATE_KWH_PER_SECOND = 10.42 / 60.0; // 10.42 kWh/min = 0.1737 kWh/s (default)
    public static final double VERTICAL_ENERGY_CONSUMPTION_KWH_PER_PASSENGER_PER_METER = 0.485 / 1000.0; // 0.485 kWh/pkm
    public static final double HORIZONTAL_ENERGY_CONSUMPTION_KWH_PER_PASSENGER_PER_METER = 6.961 / 1000.0; // 6.961 kWh/pkm
    public static final double EVTOL_ALTITUDE_METERS = 600.0; // From literature

    private double currentBatteryLevel; // in kWh
    private final double chargingRateKwhPerSecond; // Instance-specific charging rate

    public EVTOLBatteryManager() {
        this.currentBatteryLevel = BATTERY_CAPACITY_KWH; // Start fully charged
        this.chargingRateKwhPerSecond = DEFAULT_CHARGING_RATE_KWH_PER_SECOND;
    }

    public EVTOLBatteryManager(double initialBatteryLevel, double chargingRateKwhPerSecond) {
        this.currentBatteryLevel = Math.min(initialBatteryLevel, BATTERY_CAPACITY_KWH);
        this.chargingRateKwhPerSecond = chargingRateKwhPerSecond;
    }

    /**
     * Calculate energy consumption for a trip
     * @param origin Starting coordinate
     * @param destination Ending coordinate
     * @param numPassengers Number of passengers
     * @return Energy consumption in kWh
     */
    public double calculateEnergyConsumption(Coord origin, Coord destination, int numPassengers) {
        // Calculate horizontal distance
        double horizontalDistance = calculateHorizontalDistance(origin, destination);

        // Vertical distance is 2 * altitude (up and down)
        double verticalDistance = 2 * EVTOL_ALTITUDE_METERS;

        // Calculate energy consumption
        double horizontalEnergy = horizontalDistance * HORIZONTAL_ENERGY_CONSUMPTION_KWH_PER_PASSENGER_PER_METER;
        double verticalEnergy = 0;
        if (horizontalDistance!=0){
        verticalEnergy = verticalDistance * VERTICAL_ENERGY_CONSUMPTION_KWH_PER_PASSENGER_PER_METER;
        }

        return horizontalEnergy + verticalEnergy;
    }

    /**
     * Calculate energy consumption for deadheading (repositioning without passengers)
     * @param origin Starting coordinate
     * @param destination Ending coordinate
     * @return Energy consumption in kWh
     */
    public double calculateDeadheadingEnergyConsumption(Coord origin, Coord destination) {
        return calculateEnergyConsumption(origin, destination, 0); // No passengers, but aircraft still consumes energy
    }

    /**
     * Check if there's enough battery for a trip
     * @param origin Starting coordinate
     * @param destination Ending coordinate
     * @param numPassengers Number of passengers
     * @return true if battery is sufficient
     */
    public boolean hasSufficientBattery(Coord origin, Coord destination, int numPassengers) {
        double requiredEnergy = calculateEnergyConsumption(origin, destination, numPassengers);
        return currentBatteryLevel >= requiredEnergy;
    }

    /**
     * Check if there's enough battery for a trip including deadheading to next trip
     * @param currentTrip Current trip to execute
     * @param nextTripOrigin Origin of next trip (for deadheading calculation)
     * @return true if battery is sufficient for both current trip and deadheading
     */
    public boolean hasSufficientBatteryWithDeadheading(VehicleTrip currentTrip, Coord nextTripOrigin) {
        double currentTripEnergy = calculateEnergyConsumption(
                currentTrip.getOrigin(),
                currentTrip.getDestination(),
                currentTrip.getTotalPassengers()
        );

        double deadheadingEnergy = calculateDeadheadingEnergyConsumption(
                currentTrip.getDestination(),
                nextTripOrigin
        );

        return currentBatteryLevel >= (currentTripEnergy + deadheadingEnergy);
    }

    /**
     * Execute a trip (consume battery)
     * @param origin Starting coordinate
     * @param destination Ending coordinate
     * @param numPassengers Number of passengers
     * @return true if trip was executed successfully
     */
    public boolean executeTrip(Coord origin, Coord destination, int numPassengers) {
        double requiredEnergy = calculateEnergyConsumption(origin, destination, numPassengers);

        if (currentBatteryLevel >= requiredEnergy) {
            currentBatteryLevel -= requiredEnergy;
            return true;
        }
        return false;
    }

    /**
     * Execute deadheading (repositioning flight)
     * @param origin Starting coordinate
     * @param destination Ending coordinate
     * @return true if deadheading was executed successfully
     */
    public boolean executeDeadheading(Coord origin, Coord destination) {
        double requiredEnergy = calculateDeadheadingEnergyConsumption(origin, destination);

        if (currentBatteryLevel >= requiredEnergy) {
            currentBatteryLevel -= requiredEnergy;
            return true;
        }
        return false;
    }

    /**
     * Charge battery during connection time
     * @param connectionTimeSeconds Time available for charging in seconds
     * @return Amount of energy charged in kWh
     */
    public double chargeBattery(long connectionTimeSeconds) {
        double maxPossibleCharge = connectionTimeSeconds * chargingRateKwhPerSecond;
        double actualCharge = Math.min(maxPossibleCharge, BATTERY_CAPACITY_KWH - currentBatteryLevel);

        currentBatteryLevel += actualCharge;
        return actualCharge;
    }

    /**
     * Get time needed to charge to full capacity
     * @return Time in seconds needed for full charge
     */
    public long getTimeToFullCharge() {
        double energyNeeded = BATTERY_CAPACITY_KWH - currentBatteryLevel;
        return (long) (energyNeeded / chargingRateKwhPerSecond);
    }

    /**
     * Get minimum charging time needed for a specific trip
     * @param origin Starting coordinate of next trip
     * @param destination Ending coordinate of next trip
     * @param numPassengers Number of passengers for next trip
     * @return Minimum charging time in seconds, or -1 if trip is impossible even with full battery
     */
    public long getMinChargingTimeForTrip(Coord origin, Coord destination, int numPassengers) {
        double requiredEnergy = calculateEnergyConsumption(origin, destination, numPassengers);

        if (requiredEnergy > BATTERY_CAPACITY_KWH) {
            return -1; // Trip impossible even with full battery
        }

        if (currentBatteryLevel >= requiredEnergy) {
            return 0; // No charging needed
        }

        double energyToCharge = requiredEnergy - currentBatteryLevel;
        return (long) (energyToCharge / chargingRateKwhPerSecond);
    }

    private double calculateHorizontalDistance(Coord c1, Coord c2) {
        double dx = c1.getX() - c2.getX();
        double dy = c1.getY() - c2.getY();
        return Math.sqrt(dx * dx + dy * dy);
    }

    // Getters and setters
    public double getCurrentBatteryLevel() {
        return currentBatteryLevel;
    }

    public void setCurrentBatteryLevel(double batteryLevel) {
        this.currentBatteryLevel = Math.min(batteryLevel, BATTERY_CAPACITY_KWH);
    }

    public double getBatteryPercentage() {
        return (currentBatteryLevel / BATTERY_CAPACITY_KWH) * 100.0;
    }

    public boolean isBatteryFull() {
        return Math.abs(currentBatteryLevel - BATTERY_CAPACITY_KWH) < 0.01; // Allow small floating point errors
    }

    public EVTOLBatteryManager copy() {
        return new EVTOLBatteryManager(this.currentBatteryLevel, this.chargingRateKwhPerSecond);
    }
    
    /**
     * Get the charging rate for this instance
     * @return Instance charging rate in kWh per second
     */
    public double getChargingRateKwhPerSecond() {
        return chargingRateKwhPerSecond;
    }
}