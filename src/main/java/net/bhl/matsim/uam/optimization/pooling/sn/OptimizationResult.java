package net.bhl.matsim.uam.optimization.pooling.sn;

import java.util.List;

public class OptimizationResult {
    private final List<List<UAMTrip>> vehicleRoutes;
    private final double totalFlightDistance;
    private final int fleetSize;

    public OptimizationResult(List<List<UAMTrip>> vehicleRoutes,
                              double totalFlightDistance,
                              int fleetSize) {
        this.vehicleRoutes = vehicleRoutes;
        this.totalFlightDistance = totalFlightDistance;
        this.fleetSize = fleetSize;
    }

    // Getters

    public void printSummary() {
        System.out.println("Optimization Results:");
        System.out.println("Fleet Size: " + fleetSize);
        System.out.println("Total Flight Distance: " + totalFlightDistance);
        System.out.println("\nVehicle Routes:");

        for (int i = 0; i < vehicleRoutes.size(); i++) {
            System.out.println("\nVehicle " + (i + 1) + ":");
            List<UAMTrip> route = vehicleRoutes.get(i);

            for (UAMTrip trip : route) {
                if (trip.isPooledTrip()) {
                    System.out.println("  Pooled Trip " + trip.getId() + ":");
                    for (UAMTrip pooledTrip : trip.getPooledTrips()) {
                        System.out.println("    - Trip " + pooledTrip.getId());
                    }
                } else {
                    System.out.println("  Trip " + trip.getId());
                }
            }
        }
    }
}
