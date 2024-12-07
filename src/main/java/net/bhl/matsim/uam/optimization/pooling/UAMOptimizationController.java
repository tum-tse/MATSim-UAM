package net.bhl.matsim.uam.optimization.pooling;

import org.matsim.api.core.v01.Coord;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

public class UAMOptimizationController {
    private final List<UAMTrip> originalTrips;
    private final double maxDetourRatio;
    private final int maxPassengersPerVehicle;
    private final int maxConnectionTimeMinutes;
    private final double flightSpeedMetersPerSecond;

    public UAMOptimizationController(List<UAMTrip> trips, double maxDetourRatio,
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
        List<UAMTrip> pooledTrips = createTripPools();

        // Step 2: Build shareability network with pooled trips
        ShareabilityNetwork network = new ShareabilityNetwork(
                pooledTrips,
                maxConnectionTimeMinutes,
                flightSpeedMetersPerSecond
        );

        // Step 3: Find optimal vehicle assignments
        List<List<UAMTrip>> vehicleRoutes = network.findOptimalVehicleAssignments();

        // Step 4: Calculate metrics
        double totalFlightDistance = network.calculateTotalFlightDistance(vehicleRoutes);
        int fleetSize = vehicleRoutes.size();

        return new OptimizationResult(vehicleRoutes, totalFlightDistance, fleetSize);
    }

    private List<UAMTrip> createTripPools() {
        List<UAMTrip> result = new ArrayList<>();
        List<UAMTrip> unassignedTrips = new ArrayList<>(originalTrips);

        // Sort trips by departure time
        unassignedTrips.sort(Comparator.comparingLong(UAMTrip::getDepartureTime));

        while (!unassignedTrips.isEmpty()) {
            UAMTrip currentTrip = unassignedTrips.remove(0);
            TripPool currentPool = new TripPool(maxPassengersPerVehicle, maxDetourRatio);
            currentPool.addTrip(currentTrip);

            // Try to find compatible trips to add to the pool
            Iterator<UAMTrip> iterator = unassignedTrips.iterator();
            while (iterator.hasNext()) {
                UAMTrip candidate = iterator.next();
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
