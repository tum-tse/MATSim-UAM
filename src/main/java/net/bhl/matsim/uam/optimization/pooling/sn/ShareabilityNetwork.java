package net.bhl.matsim.uam.optimization.pooling.sn;

import org.matsim.api.core.v01.Coord;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Enhanced ShareabilityNetwork with eVTOL battery management
 * Core classes for UAM shareability network implementation with battery constraints
 */
public class ShareabilityNetwork {
    private List<VehicleTrip> trips;
    private Map<String, Set<String>> adjacencyList;
    private Map<String, EVTOLVehicle> vehiclePool; // Available vehicles
    private int maxConnectionTimeMinutes;
    private double flightSpeedMetersPerSecond;
    private double chargingRateKwhPerSecond; // Charging rate for vehicles
    private int nextVehicleId;
    private BatteryStatistics lastBatteryStatistics; // Store statistics from last optimization

    public ShareabilityNetwork(List<VehicleTrip> trips, int maxConnectionTimeMinutes,
                               double flightSpeedMetersPerSecond, double chargingRateKwhPerSecond) {
        this.trips = new ArrayList<>(trips);
        this.adjacencyList = new HashMap<>();
        this.vehiclePool = new HashMap<>();
        this.maxConnectionTimeMinutes = maxConnectionTimeMinutes;
        this.flightSpeedMetersPerSecond = flightSpeedMetersPerSecond;
        this.chargingRateKwhPerSecond = chargingRateKwhPerSecond;
        this.nextVehicleId = 1;
        this.lastBatteryStatistics = null;

        buildBatteryAwareNetwork();
    }

    /**
     * Build network considering both time and battery constraints
     */
    private void buildBatteryAwareNetwork() {
        // Sort trips by departure time
        trips.sort(Comparator.comparingLong(VehicleTrip::getDepartureTime));

        // Build adjacency list - connect trips that can be served sequentially with battery constraints
        for (int i = 0; i < trips.size(); i++) {
            VehicleTrip t1 = trips.get(i);
            adjacencyList.put(t1.getId(), new HashSet<>());

            for (int j = i + 1; j < trips.size(); j++) {
                VehicleTrip t2 = trips.get(j);

                // Check if t2 can be served after t1 considering battery constraints
                if (canServeSequentiallyWithBattery(t1, t2)) {
                    adjacencyList.get(t1.getId()).add(t2.getId());
                }
            }
        }
    }

    /**
     * Enhanced method to check if two trips can be served sequentially considering battery constraints
     */
    private boolean canServeSequentiallyWithBattery(VehicleTrip t1, VehicleTrip t2) {
        // First check time constraints (original logic)
        double distanceMeters = calculateDistance(t1.getDestination(), t2.getOrigin());
        long flightTimeSeconds = (long) (distanceMeters / flightSpeedMetersPerSecond);
        long earliestArrival = t1.getArrivalTime() + flightTimeSeconds;
        long connectionTime = (t2.getDepartureTime() - t1.getArrivalTime()) / 60; // Convert to minutes

        if (connectionTime > maxConnectionTimeMinutes || earliestArrival > t2.getDepartureTime()) {
            return false;
        }

        // Check battery constraints with a temporary vehicle
        EVTOLVehicle tempVehicle = new EVTOLVehicle("TEMP", chargingRateKwhPerSecond);
        long connectionTimeSeconds = t2.getDepartureTime() - t1.getArrivalTime();

        return tempVehicle.canExecuteConsecutiveTrips(t1, t2, connectionTimeSeconds);
    }

    /**
     * Enhanced optimal vehicle assignment considering battery constraints
     */
    public List<List<VehicleTrip>> findOptimalVehicleAssignments() {
        List<List<VehicleTrip>> vehicleRoutes = new ArrayList<>();
        Set<String> unassignedTrips = new HashSet<>(adjacencyList.keySet());

        // Reset vehicle pool for new optimization
        vehiclePool.clear();
        nextVehicleId = 1;

        while (!unassignedTrips.isEmpty()) {
            // Find the best route considering battery constraints
            BatteryAwareRouteResult routeResult = findBestBatteryAwareRoute(unassignedTrips);

            if (routeResult.route.isEmpty()) {
                // If no valid route found, assign single trip to new vehicle
                String tripId = unassignedTrips.iterator().next();
                VehicleTrip singleTrip = getTripById(tripId);
                if (singleTrip != null) {
                    vehicleRoutes.add(Arrays.asList(singleTrip));
                    unassignedTrips.remove(tripId);
                }
                continue;
            }

            // Add the found route
            vehicleRoutes.add(routeResult.route);

            // Remove assigned trips from unassigned set
            for (VehicleTrip trip : routeResult.route) {
                unassignedTrips.remove(trip.getId());
            }
        }

        // Automatically calculate and store battery statistics
        this.lastBatteryStatistics = getBatteryStatistics(vehicleRoutes);

        return vehicleRoutes;
    }

    /**
     * Find the best route considering battery constraints
     */
    private BatteryAwareRouteResult findBestBatteryAwareRoute(Set<String> availableTrips) {
        BatteryAwareRouteResult bestResult = new BatteryAwareRouteResult();

        // Try starting from each available trip
        for (String startTripId : availableTrips) {
            VehicleTrip startTrip = getTripById(startTripId);
            if (startTrip == null) continue;

            // Try to find the longest valid path from this starting trip
            BatteryAwareRouteResult currentResult = findLongestBatteryAwarePath(startTrip, availableTrips);

            // Choose the route with maximum trips (or other optimization criteria)
            if (currentResult.route.size() > bestResult.route.size()) {
                bestResult = currentResult;
            }
        }

        return bestResult;
    }

    /**
     * Find the longest path starting from a specific trip, considering battery constraints
     */
    private BatteryAwareRouteResult findLongestBatteryAwarePath(VehicleTrip startTrip, Set<String> availableTrips) {
        EVTOLVehicle vehicle = new EVTOLVehicle("V" + nextVehicleId++, chargingRateKwhPerSecond);
        List<VehicleTrip> route = new ArrayList<>();
        Set<String> remainingTrips = new HashSet<>(availableTrips);

        VehicleTrip currentTrip = startTrip;
        route.add(currentTrip);
        remainingTrips.remove(currentTrip.getId());

        // Simulate vehicle executing the first trip
        vehicle.setCurrentLocation(currentTrip.getOrigin());
        vehicle.setCurrentTime(currentTrip.getDepartureTime());
        vehicle.executeTrip(currentTrip);

        boolean foundNext = true;
        while (foundNext && !remainingTrips.isEmpty()) {
            foundNext = false;
            VehicleTrip bestNextTrip = null;
            long bestConnectionTime = Long.MAX_VALUE;

            // Find the best next trip that can be executed
            for (String nextTripId : remainingTrips) {
                if (!adjacencyList.get(currentTrip.getId()).contains(nextTripId)) {
                    continue; // Not connected in adjacency list
                }

                VehicleTrip nextTrip = getTripById(nextTripId);
                if (nextTrip == null) continue;

                long connectionTime = nextTrip.getDepartureTime() - currentTrip.getArrivalTime();

                // Check if vehicle can execute this trip after charging
                if (vehicle.canExecuteTripAfterCharging(nextTrip, connectionTime)) {
                    // Prefer trips with shorter connection time (less idle time)
                    if (connectionTime < bestConnectionTime) {
                        bestNextTrip = nextTrip;
                        bestConnectionTime = connectionTime;
                    }
                }
            }

            if (bestNextTrip != null) {
                // Charge during connection time
                vehicle.chargeBattery(bestConnectionTime);

                // Execute the trip
                if (vehicle.executeTrip(bestNextTrip)) {
                    route.add(bestNextTrip);
                    remainingTrips.remove(bestNextTrip.getId());
                    currentTrip = bestNextTrip;
                    foundNext = true;
                } else {
                    // This should not happen if canExecuteTripAfterCharging returned true
                    System.out.println("Warning: Could not execute trip after battery check passed");
                    break;
                }
            }
        }

        return new BatteryAwareRouteResult(route, vehicle);
    }

    /**
     * Alternative assignment method that tries to assign trips to existing vehicles first
     */
    public List<List<VehicleTrip>> findOptimalVehicleAssignmentsWithReuse() { // TODO: Need to be optimized!
        List<EVTOLVehicle> vehicles = new ArrayList<>();
        List<VehicleTrip> unassignedTrips = new ArrayList<>(trips);

        // Sort trips by departure time
        unassignedTrips.sort(Comparator.comparingLong(VehicleTrip::getDepartureTime));

        for (VehicleTrip trip : unassignedTrips) {
            boolean assigned = false;

            // Try to assign to existing vehicle
            for (EVTOLVehicle vehicle : vehicles) {
                if (canAssignTripToVehicle(vehicle, trip)) {
                    assignTripToVehicle(vehicle, trip);
                    assigned = true;
                    break;
                }
            }

            // If not assigned, create new vehicle
            if (!assigned) {
                EVTOLVehicle newVehicle = new EVTOLVehicle("V" + nextVehicleId++, chargingRateKwhPerSecond);
                newVehicle.setCurrentLocation(trip.getOrigin());
                newVehicle.setCurrentTime(trip.getDepartureTime());
                newVehicle.executeTrip(trip);
                vehicles.add(newVehicle);
            }
        }

        // Convert vehicles to routes
        List<List<VehicleTrip>> vehicleRoutes = vehicles.stream()
                .map(v -> v.getAssignedTrips())
                .collect(Collectors.toList());

        // Automatically calculate and store battery statistics
        this.lastBatteryStatistics = getBatteryStatistics(vehicleRoutes);

        return vehicleRoutes;
    }

    /**
     * Check if a trip can be assigned to a vehicle considering battery and time constraints
     */
    private boolean canAssignTripToVehicle(EVTOLVehicle vehicle, VehicleTrip trip) {
        if (vehicle.getAssignedTrips().isEmpty()) {
            return vehicle.canExecuteTrip(trip);
        }

        // Get last assigned trip
        List<VehicleTrip> assignedTrips = vehicle.getAssignedTrips();
        VehicleTrip lastTrip = assignedTrips.get(assignedTrips.size() - 1);

        // Check time constraints
        if (trip.getDepartureTime() <= lastTrip.getArrivalTime()) {
            return false; // Trip starts before vehicle is available
        }

        long connectionTime = trip.getDepartureTime() - lastTrip.getArrivalTime();
        if (connectionTime > maxConnectionTimeMinutes * 60) {
            return false; // Connection time too long
        }

        // Check battery constraints
        return vehicle.canExecuteTripAfterCharging(trip, connectionTime);
    }

    /**
     * Assign a trip to a vehicle and update vehicle state
     */
    private void assignTripToVehicle(EVTOLVehicle vehicle, VehicleTrip trip) {
        // Calculate charging time
        List<VehicleTrip> assignedTrips = vehicle.getAssignedTrips();
        if (!assignedTrips.isEmpty()) {
            VehicleTrip lastTrip = assignedTrips.get(assignedTrips.size() - 1);
            long chargingTime = trip.getDepartureTime() - lastTrip.getArrivalTime();
            vehicle.chargeBattery(chargingTime);
        }

        // Execute the trip
        vehicle.executeTrip(trip);
    }

    /**
     * Get statistics about battery usage and charging
     */
    public BatteryStatistics getBatteryStatistics(List<List<VehicleTrip>> routes) {
        BatteryStatistics stats = new BatteryStatistics();

        for (List<VehicleTrip> route : routes) {
            if (route.isEmpty()) continue;

            // Create a temporary vehicle to simulate the route
            EVTOLVehicle tempVehicle = new EVTOLVehicle("STATS_VEHICLE", chargingRateKwhPerSecond);
            double totalEnergyUsed = 0;
            double totalChargingTime = 0;
            long totalChargingTimeSeconds = 0;

            for (int i = 0; i < route.size(); i++) {
                VehicleTrip trip = route.get(i);

                // Calculate deadheading energy if needed (repositioning to trip start)
                if (tempVehicle.getCurrentLocation() != null &&
                        !tempVehicle.getCurrentLocation().equals(trip.getOrigin())) {
                    double deadheadEnergy = tempVehicle.getBatteryManager()
                            .calculateDeadheadingEnergyConsumption(
                                    tempVehicle.getCurrentLocation(), trip.getOrigin());
                    totalEnergyUsed += deadheadEnergy;

                    // Execute deadheading
                    tempVehicle.getBatteryManager().executeDeadheading(
                            tempVehicle.getCurrentLocation(), trip.getOrigin());
                }

                // Calculate trip energy consumption
                double tripEnergy = tempVehicle.getBatteryManager()
                        .calculateEnergyConsumption(trip.getOrigin(), trip.getDestination(),
                                trip.getTotalPassengers());
                totalEnergyUsed += tripEnergy;

                // Execute the trip
                tempVehicle.executeTrip(trip);

                // Calculate charging time between trips (if not the last trip)
                if (i < route.size() - 1) {
                    VehicleTrip nextTrip = route.get(i + 1);
                    long connectionTimeSeconds = nextTrip.getDepartureTime() - trip.getArrivalTime();
                    totalChargingTimeSeconds += connectionTimeSeconds;

                    // Simulate charging
                    tempVehicle.chargeBattery(connectionTimeSeconds);
                }
            }

            // Convert charging time to minutes
            totalChargingTime = totalChargingTimeSeconds / 60.0;

            stats.addVehicleStats(totalEnergyUsed, totalChargingTime);
        }

        return stats;
    }

    private double calculateDistance(Coord l1, Coord l2) {
        double dx = l1.getX() - l2.getX();
        double dy = l1.getY() - l2.getY();
        return Math.sqrt(dx * dx + dy * dy);
    }

    private VehicleTrip getTripById(String id) {
        return trips.stream()
                .filter(t -> t.getId().equals(id))
                .findFirst()
                .orElse(null);
    }

    // Getter methods
    public List<VehicleTrip> getTrips() {
        return Collections.unmodifiableList(trips);
    }

    public Map<String, Set<String>> getAdjacencyList() {
        return Collections.unmodifiableMap(adjacencyList);
    }

    /**
     * Get the battery statistics from the last optimization
     * @return BatteryStatistics object, or null if no optimization has been run
     */
    public BatteryStatistics getLastBatteryStatistics() {
        return lastBatteryStatistics;
    }

    /**
     * Helper class for route finding results
     */
    private static class BatteryAwareRouteResult {
        List<VehicleTrip> route;
        EVTOLVehicle vehicle;

        public BatteryAwareRouteResult() {
            this.route = new ArrayList<>();
            this.vehicle = null;
        }

        public BatteryAwareRouteResult(List<VehicleTrip> route, EVTOLVehicle vehicle) {
            this.route = new ArrayList<>(route);
            this.vehicle = vehicle;
        }
    }

    /**
     * Helper class for battery statistics
     */
    public static class BatteryStatistics {
        private int totalVehicles = 0;
        private double totalEnergyUsed = 0;
        private double totalChargingTimeMinutes = 0;
        private double maxEnergyUsedByVehicle = 0;
        private double minEnergyUsedByVehicle = Double.MAX_VALUE;

        public void addVehicleStats(double energyUsed, double chargingTimeMinutes) {
            totalVehicles++;
            totalEnergyUsed += energyUsed;
            totalChargingTimeMinutes += chargingTimeMinutes;
            maxEnergyUsedByVehicle = Math.max(maxEnergyUsedByVehicle, energyUsed);
            if (minEnergyUsedByVehicle == Double.MAX_VALUE) {
                minEnergyUsedByVehicle = energyUsed;
            } else {
                minEnergyUsedByVehicle = Math.min(minEnergyUsedByVehicle, energyUsed);
            }
        }

        public double getAverageEnergyPerVehicle() {
            return totalVehicles > 0 ? totalEnergyUsed / totalVehicles : 0;
        }

        public double getAverageChargingTimePerVehicle() {
            return totalVehicles > 0 ? totalChargingTimeMinutes / totalVehicles : 0;
        }

        // Getters
        public int getTotalVehicles() { return totalVehicles; }
        public double getTotalEnergyUsed() { return totalEnergyUsed; }
        public double getTotalChargingTimeMinutes() { return totalChargingTimeMinutes; }
        public double getMaxEnergyUsedByVehicle() { return maxEnergyUsedByVehicle; }
        public double getMinEnergyUsedByVehicle() {
            return minEnergyUsedByVehicle == Double.MAX_VALUE ? 0.0 : minEnergyUsedByVehicle;
        }

        public void printStatistics() {
            System.out.println("=== eVTOL Battery Statistics ===");
            System.out.println("Total vehicles: " + totalVehicles);
            System.out.printf("Total energy used: %.2f kWh\n", totalEnergyUsed);
            System.out.printf("Average energy per vehicle: %.2f kWh (%.1f%% of capacity)\n",
                    getAverageEnergyPerVehicle(),
                    (getAverageEnergyPerVehicle() / EVTOLBatteryManager.BATTERY_CAPACITY_KWH) * 100);
            System.out.printf("Max energy used by single vehicle: %.2f kWh (%.1f%% of capacity)\n",
                    maxEnergyUsedByVehicle,
                    (maxEnergyUsedByVehicle / EVTOLBatteryManager.BATTERY_CAPACITY_KWH) * 100);

            if (totalVehicles > 1) {
                System.out.printf("Min energy used by single vehicle: %.2f kWh (%.1f%% of capacity)\n",
                        getMinEnergyUsedByVehicle(),
                        (getMinEnergyUsedByVehicle() / EVTOLBatteryManager.BATTERY_CAPACITY_KWH) * 100);
            }

            System.out.printf("Total charging time: %.2f minutes (%.2f hours)\n",
                    totalChargingTimeMinutes, totalChargingTimeMinutes / 60.0);
            System.out.printf("Average charging time per vehicle: %.2f minutes\n", getAverageChargingTimePerVehicle());

            // Battery efficiency insights
            double avgBatteryUtilization = (getAverageEnergyPerVehicle() / EVTOLBatteryManager.BATTERY_CAPACITY_KWH) * 100;
            if (avgBatteryUtilization > 90) {
                System.out.println("⚠️  High battery utilization - consider more vehicles or longer charging times");
            } else if (avgBatteryUtilization < 30) {
                System.out.println("ℹ️  Low battery utilization - opportunity for route consolidation");
            } else {
                System.out.println("✅ Good battery utilization");
            }
        }
    }
}