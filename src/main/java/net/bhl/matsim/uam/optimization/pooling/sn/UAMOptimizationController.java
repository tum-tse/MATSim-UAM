package net.bhl.matsim.uam.optimization.pooling.sn;

import java.util.*;

import net.bhl.matsim.uam.optimization.utils.TripItemForOptimization;
import org.matsim.api.core.v01.Coord;

public class UAMOptimizationController {
    private final List<VehicleTrip> vehicleTrips;
    private final double maxDetourRatio;
    private final int maxPassengersPerVehicle;
    private final int maxConnectionTimeMinutes;
    private final double flightSpeedMetersPerSecond;

    // Constructor for non-pooled trips
    public UAMOptimizationController(List<VehicleTrip> trips, double maxDetourRatio,
                                     int maxPassengersPerVehicle, int maxConnectionTimeMinutes,
                                     double flightSpeedMetersPerSecond) {
        this.maxDetourRatio = maxDetourRatio;
        this.maxPassengersPerVehicle = maxPassengersPerVehicle;
        this.maxConnectionTimeMinutes = maxConnectionTimeMinutes;
        this.flightSpeedMetersPerSecond = flightSpeedMetersPerSecond;

        // Create trip pools from individual trips
        this.vehicleTrips = createTripPools(trips);
    }

    // Constructor for already pooled vehicle assignments
    public UAMOptimizationController(Map<Integer, List<TripItemForOptimization>> vehicleAssignments,
                                     double maxDetourRatio,
                                     int maxPassengersPerVehicle,
                                     int maxConnectionTimeMinutes,
                                     double flightSpeedMetersPerSecond) {
        this.maxDetourRatio = maxDetourRatio;
        this.maxPassengersPerVehicle = maxPassengersPerVehicle;
        this.maxConnectionTimeMinutes = maxConnectionTimeMinutes;
        this.flightSpeedMetersPerSecond = flightSpeedMetersPerSecond;

        // Convert vehicle assignments directly to vehicle trips
        this.vehicleTrips = convertAssignmentsToTrips(vehicleAssignments);
    }

    public OptimizationResult optimize() {
        // Build shareability network with existing vehicle trips
        ShareabilityNetwork network = new ShareabilityNetwork(
                vehicleTrips,
                maxConnectionTimeMinutes,
                flightSpeedMetersPerSecond
        );

        // Find optimal vehicle assignments
        List<List<VehicleTrip>> vehicleRoutes = network.findOptimalVehicleAssignments();

        // Calculate metrics
        //double totalFlightDistance = network.calculateTotalFlightDistance(vehicleRoutes);
        double totalDeadheadingFlightDistance = network.calculateTotalDeadheadingFlightDistance(vehicleRoutes);
        int fleetSize = vehicleRoutes.size();

        return new OptimizationResult(vehicleRoutes,
                //totalFlightDistance,
                totalDeadheadingFlightDistance,
                fleetSize);
    }

    private List<VehicleTrip> createTripPools(List<VehicleTrip> trips) {
        List<VehicleTrip> result = new ArrayList<>();
        List<VehicleTrip> unassignedTrips = new ArrayList<>(trips);

        // Sort trips by departure time
        unassignedTrips.sort(Comparator.comparingLong(VehicleTrip::getDepartureTime));

        while (!unassignedTrips.isEmpty()) {
            VehicleTrip currentTrip = unassignedTrips.remove(0);
            TripPool currentPool = new TripPool(maxPassengersPerVehicle, maxDetourRatio);
            currentPool.addTrip(currentTrip);

            // Try to find compatible trips to add to the pool
            Iterator<VehicleTrip> iterator = unassignedTrips.iterator();
            while (iterator.hasNext()) {
                VehicleTrip candidate = iterator.next();
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

            // Use first trip's origin/destination as reference
            TripItemForOptimization firstTrip = assignedTrips.get(0);
            Coord origin = new Coord(firstTrip.accessVertiport.coord.getX(),
                    firstTrip.accessVertiport.coord.getY());
            Coord destination = new Coord(firstTrip.egressVertiport.coord.getX(),
                    firstTrip.egressVertiport.coord.getY());

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

    private double calculateDistance(Coord l1, Coord l2) {
        double dx = l1.getX() - l2.getX();
        double dy = l1.getY() - l2.getY();
        return Math.sqrt(dx * dx + dy * dy);
    }
}
