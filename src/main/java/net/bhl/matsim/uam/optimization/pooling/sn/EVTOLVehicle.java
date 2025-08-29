package net.bhl.matsim.uam.optimization.pooling.sn;

import org.matsim.api.core.v01.Coord;

import java.util.*;

/**
 * Represents an eVTOL vehicle with comprehensive battery management capabilities
 * This class handles vehicle state, battery management, and route execution validation
 */
public class EVTOLVehicle {

    private String vehicleId;
    private EVTOLBatteryManager batteryManager;
    private List<VehicleTrip> assignedTrips;
    private Coord currentLocation;
    private long currentTime; // Current simulation time in seconds
    private boolean isAvailable;

    public EVTOLVehicle(String vehicleId) {
        this.vehicleId = vehicleId;
        this.batteryManager = new EVTOLBatteryManager(); // Start with full battery
        this.assignedTrips = new ArrayList<>();
        this.currentLocation = null;
        this.currentTime = 0;
        this.isAvailable = true;
    }

    public EVTOLVehicle(String vehicleId, Coord initialLocation, long initialTime) {
        this.vehicleId = vehicleId;
        this.batteryManager = new EVTOLBatteryManager();
        this.assignedTrips = new ArrayList<>();
        this.currentLocation = initialLocation;
        this.currentTime = initialTime;
        this.isAvailable = true;
    }

    /**
     * Check if vehicle can execute a trip considering current battery state
     * @param trip Trip to evaluate
     * @return true if vehicle can execute the trip
     */
    public boolean canExecuteTrip(VehicleTrip trip) {
        // Check if vehicle needs to reposition first (deadheading)
        if (currentLocation != null && !currentLocation.equals(trip.getOrigin())) {
            return batteryManager.hasSufficientBattery(currentLocation, trip.getOrigin(), 0) && // Deadheading
                    batteryManager.hasSufficientBattery(trip.getOrigin(), trip.getDestination(), trip.getTotalPassengers()); // Actual trip
        }

        return batteryManager.hasSufficientBattery(trip.getOrigin(), trip.getDestination(), trip.getTotalPassengers());
    }

    /**
     * Check if vehicle can execute a trip after charging for available connection time
     * @param trip Trip to evaluate
     * @param connectionTimeSeconds Available time for charging
     * @return true if vehicle can execute trip after charging
     */
    public boolean canExecuteTripAfterCharging(VehicleTrip trip, long connectionTimeSeconds) {
        // Create a copy of battery manager to simulate charging
        EVTOLBatteryManager tempBattery = batteryManager.copy();
        tempBattery.chargeBattery(connectionTimeSeconds);

        // Check if can execute with charged battery
        if (currentLocation != null && !currentLocation.equals(trip.getOrigin())) {
            return tempBattery.hasSufficientBattery(currentLocation, trip.getOrigin(), 0) && // Deadheading
                    tempBattery.hasSufficientBattery(trip.getOrigin(), trip.getDestination(), trip.getTotalPassengers()); // Actual trip
        }

        return tempBattery.hasSufficientBattery(trip.getOrigin(), trip.getDestination(), trip.getTotalPassengers());
    }

    /**
     * Check if vehicle can execute two consecutive trips
     * @param currentTrip First trip
     * @param nextTrip Second trip
     * @param connectionTimeSeconds Time between trips for charging
     * @return true if both trips can be executed
     */
    public boolean canExecuteConsecutiveTrips(VehicleTrip currentTrip, VehicleTrip nextTrip, long connectionTimeSeconds) {
        // Create a copy to simulate the operations
        EVTOLBatteryManager tempBattery = batteryManager.copy();

        // First, check deadheading to current trip if needed
        if (currentLocation != null && !currentLocation.equals(currentTrip.getOrigin())) {
            if (!tempBattery.hasSufficientBattery(currentLocation, currentTrip.getOrigin(), 0)) {
                return false;
            }
            tempBattery.executeDeadheading(currentLocation, currentTrip.getOrigin());
        }

        // Execute current trip
        if (!tempBattery.executeTrip(currentTrip.getOrigin(), currentTrip.getDestination(), currentTrip.getTotalPassengers())) {
            return false;
        }

        // Charge during connection time
        tempBattery.chargeBattery(connectionTimeSeconds);

        // Check deadheading to next trip
        if (!tempBattery.executeDeadheading(currentTrip.getDestination(), nextTrip.getOrigin())) {
            return false;
        }

        // Check if can execute next trip
        return tempBattery.hasSufficientBattery(nextTrip.getOrigin(), nextTrip.getDestination(), nextTrip.getTotalPassengers());
    }

    /**
     * Execute a trip and update vehicle state
     * @param trip Trip to execute
     * @return true if trip was executed successfully
     */
    public boolean executeTrip(VehicleTrip trip) {
        if (!isAvailable) {
            return false;
        }

        // Handle deadheading if needed
        if (currentLocation != null && !currentLocation.equals(trip.getOrigin())) {
            if (!batteryManager.executeDeadheading(currentLocation, trip.getOrigin())) {
                return false;
            }
        }

        // Execute the actual trip
        if (batteryManager.executeTrip(trip.getOrigin(), trip.getDestination(), trip.getTotalPassengers())) {
            assignedTrips.add(trip);
            currentLocation = trip.getDestination();
            currentTime = trip.getArrivalTime();
            return true;
        }

        return false;
    }

    /**
     * Charge battery during connection/waiting time
     * @param chargingTimeSeconds Time available for charging
     * @return Amount of energy charged
     */
    public double chargeBattery(long chargingTimeSeconds) {
        return batteryManager.chargeBattery(chargingTimeSeconds);
    }

    /**
     * Get the minimum connection time needed between two trips for battery feasibility
     * @param currentTrip First trip
     * @param nextTrip Second trip
     * @return Minimum connection time in seconds, or -1 if impossible
     */
    public long getMinimumConnectionTime(VehicleTrip currentTrip, VehicleTrip nextTrip) {
        // Create simulation copy
        EVTOLBatteryManager tempBattery = batteryManager.copy();

        // Execute current trip
        if (!tempBattery.executeTrip(currentTrip.getOrigin(), currentTrip.getDestination(), currentTrip.getTotalPassengers())) {
            return -1; // Can't even execute current trip
        }

        // Calculate energy needed for deadheading to next trip
        double deadheadingEnergy = tempBattery.calculateDeadheadingEnergyConsumption(
                currentTrip.getDestination(), nextTrip.getOrigin()
        );

        // Calculate energy needed for next trip
        double nextTripEnergy = tempBattery.calculateEnergyConsumption(
                nextTrip.getOrigin(), nextTrip.getDestination(), nextTrip.getTotalPassengers()
        );

        double totalEnergyNeeded = deadheadingEnergy + nextTripEnergy;

        if (tempBattery.getCurrentBatteryLevel() >= totalEnergyNeeded) {
            return 0; // No charging needed
        }

        double energyToCharge = totalEnergyNeeded - tempBattery.getCurrentBatteryLevel();
        return (long) (energyToCharge / EVTOLBatteryManager.CHARGING_RATE_KWH_PER_SECOND);
    }

    /**
     * Reset vehicle to initial state (for optimization iterations)
     */
    public void reset() {
        this.batteryManager = new EVTOLBatteryManager(); // Full battery
        this.assignedTrips.clear();
        this.isAvailable = true;
        this.currentTime = 0;
        // Keep currentLocation as it might be set initially
    }

    /**
     * Reset vehicle to specific location and time
     */
    public void reset(Coord location, long time) {
        reset();
        this.currentLocation = location;
        this.currentTime = time;
    }

    /**
     * Get total energy consumption for all assigned trips including deadheading
     */
    public double getTotalEnergyConsumption() {
        double total = 0;
        Coord previousLocation = currentLocation;

        for (VehicleTrip trip : assignedTrips) {
            // Add deadheading energy if needed
            if (previousLocation != null && !previousLocation.equals(trip.getOrigin())) {
                total += batteryManager.calculateDeadheadingEnergyConsumption(previousLocation, trip.getOrigin());
            }

            // Add trip energy
            total += batteryManager.calculateEnergyConsumption(trip.getOrigin(), trip.getDestination(), trip.getTotalPassengers());

            previousLocation = trip.getDestination();
        }

        return total;
    }

    /**
     * Get total flight distance for all assigned trips including deadheading
     */
    public double getTotalFlightDistance() {
        double totalDistance = 0;
        Coord previousLocation = currentLocation;

        for (VehicleTrip trip : assignedTrips) {
            // Add deadheading distance if needed
            if (previousLocation != null && !previousLocation.equals(trip.getOrigin())) {
                totalDistance += calculateDistance(previousLocation, trip.getOrigin());
            }

            // Add trip distance
            totalDistance += calculateDistance(trip.getOrigin(), trip.getDestination());

            previousLocation = trip.getDestination();
        }

        return totalDistance;
    }

    /**
     * Get detailed battery status for debugging and monitoring
     */
    public String getBatteryStatus() {
        return String.format("Vehicle %s: Battery %.2f kWh (%.1f%%), Location: %s, Available: %s, Trips: %d",
                vehicleId,
                batteryManager.getCurrentBatteryLevel(),
                batteryManager.getBatteryPercentage(),
                currentLocation != null ? String.format("(%.1f,%.1f)", currentLocation.getX(), currentLocation.getY()) : "null",
                isAvailable,
                assignedTrips.size()
        );
    }

    /**
     * Print detailed route analysis for this vehicle
     */
    public void printRouteAnalysis() {
        System.out.println("=== Route Analysis for Vehicle " + vehicleId + " ===");
        System.out.println(getBatteryStatus());

        if (assignedTrips.isEmpty()) {
            System.out.println("No trips assigned.");
            return;
        }

        double totalEnergy = 0;
        Coord prevLocation = currentLocation;

        for (int i = 0; i < assignedTrips.size(); i++) {
            VehicleTrip trip = assignedTrips.get(i);

            System.out.printf("\nTrip %d: %s\n", i + 1, trip.getId());

            // Deadheading if needed
            if (prevLocation != null && !prevLocation.equals(trip.getOrigin())) {
                double deadheadDistance = calculateDistance(prevLocation, trip.getOrigin());
                double deadheadEnergy = batteryManager.calculateDeadheadingEnergyConsumption(prevLocation, trip.getOrigin());
                totalEnergy += deadheadEnergy;
                System.out.printf("  Deadheading: %.1f m, %.2f kWh\n", deadheadDistance, deadheadEnergy);
            }

            // Trip itself
            double tripDistance = calculateDistance(trip.getOrigin(), trip.getDestination());
            double tripEnergy = batteryManager.calculateEnergyConsumption(trip.getOrigin(), trip.getDestination(), trip.getTotalPassengers());
            totalEnergy += tripEnergy;

            System.out.printf("  Trip: %.1f m, %d passengers, %.2f kWh\n",
                    tripDistance, trip.getTotalPassengers(), tripEnergy);

            // Charging time between trips
            if (i < assignedTrips.size() - 1) {
                VehicleTrip nextTrip = assignedTrips.get(i + 1);
                long chargingTime = nextTrip.getDepartureTime() - trip.getArrivalTime();
                double chargedEnergy = chargingTime * EVTOLBatteryManager.CHARGING_RATE_KWH_PER_SECOND;
                System.out.printf("  Charging: %d seconds, %.2f kWh charged\n", chargingTime, chargedEnergy);
            }

            prevLocation = trip.getDestination();
        }

        System.out.printf("\nTotal energy consumption: %.2f kWh (%.1f%% of capacity)\n",
                totalEnergy, (totalEnergy / EVTOLBatteryManager.BATTERY_CAPACITY_KWH) * 100);
        System.out.printf("Final battery level: %.2f kWh (%.1f%%)\n",
                batteryManager.getCurrentBatteryLevel(), batteryManager.getBatteryPercentage());
    }

    private double calculateDistance(Coord c1, Coord c2) {
        double dx = c1.getX() - c2.getX();
        double dy = c1.getY() - c2.getY();
        return Math.sqrt(dx * dx + dy * dy);
    }

    // Getters and setters
    public String getVehicleId() {
        return vehicleId;
    }

    public EVTOLBatteryManager getBatteryManager() {
        return batteryManager;
    }

    public List<VehicleTrip> getAssignedTrips() {
        return Collections.unmodifiableList(assignedTrips);
    }

    public Coord getCurrentLocation() {
        return currentLocation;
    }

    public void setCurrentLocation(Coord currentLocation) {
        this.currentLocation = currentLocation;
    }

    public long getCurrentTime() {
        return currentTime;
    }

    public void setCurrentTime(long currentTime) {
        this.currentTime = currentTime;
    }

    public boolean isAvailable() {
        return isAvailable;
    }

    public void setAvailable(boolean available) {
        this.isAvailable = available;
    }

    public double getCurrentBatteryPercentage() {
        return batteryManager.getBatteryPercentage();
    }

    public double getCurrentBatteryLevel() {
        return batteryManager.getCurrentBatteryLevel();
    }
}