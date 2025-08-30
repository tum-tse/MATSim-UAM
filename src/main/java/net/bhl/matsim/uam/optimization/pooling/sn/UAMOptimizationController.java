package net.bhl.matsim.uam.optimization.pooling.sn;

import net.bhl.matsim.uam.optimization.Vertiport;
import net.bhl.matsim.uam.optimization.pooling.SensitivityConfig;
import net.bhl.matsim.uam.optimization.utils.TripItemForOptimization;
import org.matsim.api.core.v01.Coord;
import org.matsim.api.core.v01.Id;
import org.matsim.contrib.dvrp.fleet.DvrpVehicle;

import java.util.*;

/**
 * Enhanced UAMOptimizationController with battery-aware optimization for eVTOL operations
 * This class maintains backward compatibility while adding comprehensive battery management
 */
public class UAMOptimizationController {

    private List<VehicleTrip> vehicleTrips;
    private double maxDetourRatio;
    private int maxPassengersPerVehicle;
    private int maxConnectionTimeMinutes;
    private double flightSpeedMetersPerSecond;
    private double chargingRateKwhPerSecond; // Charging rate for vehicles
    private Map<Id<DvrpVehicle>, Vertiport> vehicleOriginStationMap;
    private Map<Id<DvrpVehicle>, Vertiport> vehicleDestinationStationMap;

    // Constructor for Map<Integer, List<TripItemForOptimization>> input (original signature)
    public UAMOptimizationController(Map<Integer, List<TripItemForOptimization>> vehicleAssignments,
                                     double maxDetourRatio,
                                     int maxPassengersPerVehicle,
                                     int maxConnectionTimeMinutes,
                                     double flightSpeedMetersPerSecond,
                                     Map<Id<DvrpVehicle>, Vertiport> vehicleOriginStationMap,
                                     Map<Id<DvrpVehicle>, Vertiport> vehicleDestinationStationMap,
                                     SensitivityConfig sensitivityConfig) {

        this.maxDetourRatio = maxDetourRatio;
        this.maxPassengersPerVehicle = maxPassengersPerVehicle;
        this.maxConnectionTimeMinutes = maxConnectionTimeMinutes;
        this.flightSpeedMetersPerSecond = flightSpeedMetersPerSecond;
        this.chargingRateKwhPerSecond = sensitivityConfig.getChargingRateKwhPerSecond();

        // Convert vehicle assignments directly to vehicle trips
        this.vehicleTrips = convertAssignmentsToTrips(vehicleAssignments);

        this.vehicleOriginStationMap = vehicleOriginStationMap;
        this.vehicleDestinationStationMap = vehicleDestinationStationMap;
    }

    // Constructor for direct VehicleTrip list input (new functionality)
    public UAMOptimizationController(List<VehicleTrip> vehicleTrips,
                                     double maxDetourRatio,
                                     int maxPassengersPerVehicle,
                                     int maxConnectionTimeMinutes,
                                     double flightSpeedMetersPerSecond) {

        this.vehicleTrips = new ArrayList<>(vehicleTrips);
        this.maxDetourRatio = maxDetourRatio;
        this.maxPassengersPerVehicle = maxPassengersPerVehicle;
        this.maxConnectionTimeMinutes = maxConnectionTimeMinutes;
        this.flightSpeedMetersPerSecond = flightSpeedMetersPerSecond;
        this.chargingRateKwhPerSecond = EVTOLBatteryManager.DEFAULT_CHARGING_RATE_KWH_PER_SECOND;
        this.vehicleOriginStationMap = null;
        this.vehicleDestinationStationMap = null;
    }

    /**
     * Main optimization method with battery-aware routing
     * This is the enhanced version that automatically includes battery statistics
     * @return OptimizationResult with battery statistics
     */
    public OptimizationResult optimize() {
        // Build shareability network with battery-aware optimization
        ShareabilityNetwork network = new ShareabilityNetwork(
                vehicleTrips,
                maxConnectionTimeMinutes,
                flightSpeedMetersPerSecond,
                chargingRateKwhPerSecond
        );

        // Find optimal vehicle assignments (automatically calculates battery statistics)
        List<List<VehicleTrip>> vehicleRoutes = network.findOptimalVehicleAssignments();

        // Get battery statistics from the network (calculated automatically during optimization)
        ShareabilityNetwork.BatteryStatistics batteryStats = network.getLastBatteryStatistics();

        // Calculate other metrics
        int fleetSize = vehicleRoutes.size();
        int vtolOperations = calculateVtolOperations(vehicleRoutes);

        // Return enhanced result with battery statistics
        return new OptimizationResult(vehicleRoutes, fleetSize, vtolOperations, batteryStats, chargingRateKwhPerSecond);
    }

    /**
     * Alternative optimization method using vehicle reuse strategy
     * @return OptimizationResult with battery statistics
     */
    public OptimizationResult optimizeWithVehicleReuse() {
        // Build shareability network with battery-aware optimization
        ShareabilityNetwork network = new ShareabilityNetwork(
                vehicleTrips,
                maxConnectionTimeMinutes,
                flightSpeedMetersPerSecond,
                chargingRateKwhPerSecond
        );

        // Use alternative assignment strategy
        List<List<VehicleTrip>> vehicleRoutes = network.findOptimalVehicleAssignmentsWithReuse();

        // Get battery statistics from the network
        ShareabilityNetwork.BatteryStatistics batteryStats = network.getLastBatteryStatistics();

        // Calculate other metrics
        int fleetSize = vehicleRoutes.size();
        int vtolOperations = calculateVtolOperations(vehicleRoutes);

        // Return enhanced result with battery statistics
        return new OptimizationResult(vehicleRoutes, fleetSize, vtolOperations, batteryStats, chargingRateKwhPerSecond);
    }

    /**
     * Legacy optimization method for backward compatibility
     * This method maintains the exact same signature as the original but now includes battery statistics
     * @return OptimizationResult with battery statistics included
     */
    public OptimizationResult optimizeWithBatteryAwareness() {
        return optimize();
    }

    /**
     * Method to calculate total VTOL operations (original implementation maintained)
     */
    private int calculateVtolOperations(List<List<VehicleTrip>> vehicleRoutes) {
        int totalOperations = 0;

        // Each route represents one vehicle
        for (List<VehicleTrip> route : vehicleRoutes) {
            if (!route.isEmpty()) {
                // Each route has at least one takeoff and one landing
                int operationsForRoute = 2;  // Initial takeoff and final landing

                // Add intermediate takeoff+landing for each connection between consecutive trips
                operationsForRoute += (route.size() - 1) * 2;

                totalOperations += operationsForRoute;
            }
        }

        return totalOperations;
    }

    /**
     * Convert vehicle assignments to VehicleTrip objects (original implementation maintained)
     */
    private List<VehicleTrip> convertAssignmentsToTrips(Map<Integer, List<TripItemForOptimization>> vehicleAssignments) {
        List<VehicleTrip> trips = new ArrayList<>();

        for (Map.Entry<Integer, List<TripItemForOptimization>> entry : vehicleAssignments.entrySet()) {
            if (entry.getValue().isEmpty()) continue;

            // Get all trips for this vehicle
            List<TripItemForOptimization> assignedTrips = entry.getValue();

            // Find latest departure time
            long latestDepartureTime = assignedTrips.stream()
                    .mapToLong(trip -> (long) trip.departureTime)
                    .max()
                    .orElse(0);

            Coord origin;
            Coord destination;
            if(vehicleOriginStationMap == null && vehicleDestinationStationMap == null) {
                // Use first trip's origin/destination as reference
                TripItemForOptimization firstTrip = assignedTrips.get(0);
                origin = new Coord(firstTrip.accessVertiport.coord.getX(),
                        firstTrip.accessVertiport.coord.getY());
                destination = new Coord(firstTrip.egressVertiport.coord.getX(),
                        firstTrip.egressVertiport.coord.getY());
            }else{
                origin = new Coord(vehicleOriginStationMap.get(Id.create(entry.getKey(), DvrpVehicle.class)).coord.getX(),
                        vehicleOriginStationMap.get(Id.create(entry.getKey(), DvrpVehicle.class)).coord.getY());
                destination = new Coord(vehicleDestinationStationMap.get(Id.create(entry.getKey(), DvrpVehicle.class)).coord.getX(),
                        vehicleDestinationStationMap.get(Id.create(entry.getKey(), DvrpVehicle.class)).coord.getY());
            }

            // Calculate arrival time
            double distance = calculateDistance(origin, destination);
            long flightTime = (long)(distance / flightSpeedMetersPerSecond);
            long arrivalTime = latestDepartureTime + flightTime;

            // Create pooled vehicle trip
            VehicleTrip vehicleTrip = new VehicleTrip(
                    "V" + entry.getKey(),
                    origin,
                    destination,
                    latestDepartureTime,
                    arrivalTime,
                    assignedTrips.size()
            );

            // Add individual trips as pooled trips
            for (TripItemForOptimization trip : assignedTrips) {
                VehicleTrip individualTrip = new VehicleTrip(
                        "T" + trip.tripID,
                        origin,
                        destination,
                        (long) trip.departureTime,
                        (long) (trip.departureTime + flightTime),
                        1
                );
                vehicleTrip.addPooledTrip(individualTrip);
            }

            trips.add(vehicleTrip);
        }

        return trips;
    }

    /**
     * Get battery statistics for a given set of routes without running full optimization
     * @param routes Vehicle routes to analyze
     * @return BatteryStatistics object
     */
    public ShareabilityNetwork.BatteryStatistics getBatteryStatisticsForRoutes(List<List<VehicleTrip>> routes) {
        // Create a temporary network just to use its battery calculation methods
        ShareabilityNetwork tempNetwork = new ShareabilityNetwork(
                vehicleTrips, maxConnectionTimeMinutes, flightSpeedMetersPerSecond, chargingRateKwhPerSecond
        );
        return tempNetwork.getBatteryStatistics(routes);
    }

    /**
     * Print battery usage analysis for current vehicle trips
     */
    public void printBatteryAnalysis() {
        System.out.println("=== eVTOL Battery Analysis for Current Trips ===");

        EVTOLBatteryManager tempBattery = new EVTOLBatteryManager(EVTOLBatteryManager.BATTERY_CAPACITY_KWH, chargingRateKwhPerSecond);
        double totalEnergyNeeded = 0;

        System.out.println("Individual trip energy requirements:");
        for (VehicleTrip trip : vehicleTrips) {
            double tripEnergy = tempBattery.calculateEnergyConsumption(
                    trip.getOrigin(), trip.getDestination(), trip.getTotalPassengers()
            );
            totalEnergyNeeded += tripEnergy;

            System.out.printf("Trip %s: %.2f kWh (Passengers: %d, Distance: %.1f m)\n",
                    trip.getId(), tripEnergy, trip.getTotalPassengers(),
                    calculateDistance(trip.getOrigin(), trip.getDestination()));
        }

        System.out.printf("\nBattery Analysis Summary:\n");
        System.out.printf("Total energy needed for all trips: %.2f kWh\n", totalEnergyNeeded);
        System.out.printf("Minimum vehicles needed (energy constraint): %d\n",
                (int) Math.ceil(totalEnergyNeeded / EVTOLBatteryManager.BATTERY_CAPACITY_KWH));
        System.out.printf("Battery capacity per vehicle: %.2f kWh\n", EVTOLBatteryManager.BATTERY_CAPACITY_KWH);
        System.out.printf("Charging rate: %.3f kWh/second (%.2f kWh/minute)\n",
                tempBattery.getChargingRateKwhPerSecond(),
                tempBattery.getChargingRateKwhPerSecond() * 60);
        System.out.printf("Vertical energy consumption: %.6f kWh/passenger/meter\n",
                EVTOLBatteryManager.VERTICAL_ENERGY_CONSUMPTION_KWH_PER_PASSENGER_PER_METER);
        System.out.printf("Horizontal energy consumption: %.6f kWh/passenger/meter\n",
                EVTOLBatteryManager.HORIZONTAL_ENERGY_CONSUMPTION_KWH_PER_PASSENGER_PER_METER);
    }

    /**
     * Print configuration summary
     */
    public void printConfiguration() {
        System.out.println("=== UAM Optimization Controller Configuration ===");
        System.out.println("Number of trips: " + vehicleTrips.size());
        System.out.println("Max detour ratio: " + maxDetourRatio);
        System.out.println("Max passengers per vehicle: " + maxPassengersPerVehicle);
        System.out.println("Max connection time: " + maxConnectionTimeMinutes + " minutes");
        System.out.println("Flight speed: " + flightSpeedMetersPerSecond + " m/s (" +
                (flightSpeedMetersPerSecond * 3.6) + " km/h)");
        System.out.println("Has origin/destination station maps: " +
                (vehicleOriginStationMap != null && vehicleDestinationStationMap != null));
        System.out.println("Battery-aware optimization: ENABLED");
    }

    /**
     * Run a complete analysis including configuration, battery analysis, and optimization
     */
    public OptimizationResult runCompleteAnalysis() {
        System.out.println("=== Complete eVTOL Optimization Analysis ===\n");

        printConfiguration();
        System.out.println();

        printBatteryAnalysis();
        System.out.println();

        System.out.println("Running battery-aware optimization...\n");
        OptimizationResult result = optimize();

        result.printSummary();

        return result;
    }

    private double calculateDistance(Coord l1, Coord l2) {
        double dx = l1.getX() - l2.getX();
        double dy = l1.getY() - l2.getY();
        return Math.sqrt(dx * dx + dy * dy);
    }

    // Getters (original interface maintained)
    public List<VehicleTrip> getVehicleTrips() {
        return Collections.unmodifiableList(vehicleTrips);
    }

    public double getMaxDetourRatio() {
        return maxDetourRatio;
    }

    public int getMaxPassengersPerVehicle() {
        return maxPassengersPerVehicle;
    }

    public int getMaxConnectionTimeMinutes() {
        return maxConnectionTimeMinutes;
    }

    public double getFlightSpeedMetersPerSecond() {
        return flightSpeedMetersPerSecond;
    }
}