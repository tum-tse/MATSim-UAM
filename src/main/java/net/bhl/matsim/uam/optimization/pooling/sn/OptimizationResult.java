package net.bhl.matsim.uam.optimization.pooling.sn;

import java.util.List;

public class OptimizationResult {
    private final List<List<VehicleTrip>> vehicleRoutes;
    //private final double totalFlightDistance;
    private final double totalDeadheadingFlightDistance;
    private final int fleetSize;

    public OptimizationResult(List<List<VehicleTrip>> vehicleRoutes,
                              //double totalFlightDistance,
                              double totalDeadheadingFlightDistance,
                              int fleetSize) {
        this.vehicleRoutes = vehicleRoutes;
        //this.totalFlightDistance = totalFlightDistance;
        this.totalDeadheadingFlightDistance = totalDeadheadingFlightDistance;
        this.fleetSize = fleetSize;
    }

    // Getters

    public void printSummary() {
        System.out.println("Optimization Results:");
        System.out.println("Fleet Size: " + fleetSize);
        //System.out.println("Total Flight Distance: " + totalFlightDistance);
        System.out.println("Total Deadheading Flight Distance: " + totalDeadheadingFlightDistance);
        System.out.println("\nVehicle Routes:");

        for (int i = 0; i < vehicleRoutes.size(); i++) {
            System.out.println("\nVehicle " + (i + 1) + ":");
            List<VehicleTrip> route = vehicleRoutes.get(i);

            for (VehicleTrip trip : route) {
                if (trip.isPooledTrip()) {
                    System.out.println("  Pooled Trip " + trip.getId() + ":");
                    for (VehicleTrip pooledTrip : trip.getPooledTrips()) {
                        System.out.println("    - Trip " + pooledTrip.getId());
                    }
                } else {
                    System.out.println("  Trip " + trip.getId());
                }
            }
        }
    }
}
