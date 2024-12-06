package net.bhl.matsim.uam.optimization.pooling;

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

// Example usage
public class Main {
    public static void main(String[] args) {
        // Create sample trips
        List<UAMTrip> trips = new ArrayList<>();
        trips.add(new UAMTrip("T1",
                new Location(0, 0), new Location(10, 10),
                0, 600, 2));
        trips.add(new UAMTrip("T2",
                new Location(1, 1), new Location(11, 11),
                100, 700, 1));
        trips.add(new UAMTrip("T3",
                new Location(20, 20), new Location(30, 30),
                800, 1400, 2));
        // Add more trips...

        // Create and run optimizer
        UAMOptimizationController optimizer = new UAMOptimizationController(
                trips,
                0.3, // maxDetourRatio
                4,   // maxPassengersPerVehicle
                30,  // maxConnectionTimeMinutes
                50.0 // flightSpeedMetersPerSecond
        );

        OptimizationResult result = optimizer.optimize();
        result.printSummary();
    }
}
