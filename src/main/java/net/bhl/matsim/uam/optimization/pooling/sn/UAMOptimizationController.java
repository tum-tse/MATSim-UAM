package net.bhl.matsim.uam.optimization.pooling.sn;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

public class UAMOptimizationController {
    private final List<PassengerTrip> originalTrips;
    private final double maxDetourRatio;
    private final int maxPassengersPerVehicle;
    private final int maxConnectionTimeMinutes;
    private final double flightSpeedMetersPerSecond;

    public UAMOptimizationController(List<PassengerTrip> trips, double maxDetourRatio,
                                     int maxPassengersPerVehicle, int maxConnectionTimeMinutes,
                                     double flightSpeedMetersPerSecond) {
        this.originalTrips = new ArrayList<>(trips);
        this.maxDetourRatio = maxDetourRatio;
        this.maxPassengersPerVehicle = maxPassengersPerVehicle;
        this.maxConnectionTimeMinutes = maxConnectionTimeMinutes;
        this.flightSpeedMetersPerSecond = flightSpeedMetersPerSecond;
    }

    public OptimizationResult optimize() {
        // Step 1: Create trip pools
        List<PassengerTrip> pooledTrips = createTripPools();

        // Step 2: Build shareability network with pooled trips
        ShareabilityNetwork network = new ShareabilityNetwork(
                pooledTrips,
                maxConnectionTimeMinutes,
                flightSpeedMetersPerSecond
        );

        // Step 3: Find optimal vehicle assignments
        List<List<PassengerTrip>> vehicleRoutes = network.findOptimalVehicleAssignments();

        // Step 4: Calculate metrics
        //double totalFlightDistance = network.calculateTotalFlightDistance(vehicleRoutes);
        double totalDeadheadingFlightDistance = network.calculateTotalDeadheadingFlightDistance(vehicleRoutes);
        int fleetSize = vehicleRoutes.size();

        return new OptimizationResult(vehicleRoutes,
                //totalFlightDistance,
                totalDeadheadingFlightDistance,
                fleetSize);
    }

    private List<PassengerTrip> createTripPools() {
        List<PassengerTrip> result = new ArrayList<>();
        List<PassengerTrip> unassignedTrips = new ArrayList<>(originalTrips);

        // Sort trips by departure time
        unassignedTrips.sort(Comparator.comparingLong(PassengerTrip::getDepartureTime));

        while (!unassignedTrips.isEmpty()) {
            PassengerTrip currentTrip = unassignedTrips.remove(0);
            TripPool currentPool = new TripPool(maxPassengersPerVehicle, maxDetourRatio);
            currentPool.addTrip(currentTrip);

            // Try to find compatible trips to add to the pool
            Iterator<PassengerTrip> iterator = unassignedTrips.iterator();
            while (iterator.hasNext()) {
                PassengerTrip candidate = iterator.next();
                if (currentPool.canAddTrip(candidate)) {
                    currentPool.addTrip(candidate);
                    iterator.remove();
                }
            }

            // Add pooled trip to result
            result.add(currentPool.toPooledTrip());
        }

        return result;
    }
}
