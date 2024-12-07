package net.bhl.matsim.uam.optimization.pooling;

import org.matsim.api.core.v01.Coord;

import java.util.ArrayList;
import java.util.List;

// Example usage
public class RunShareabilityNetwork {
    public static void main(String[] args) {
        // Create sample trips
        List<UAMTrip> trips = new ArrayList<>();
        trips.add(new UAMTrip("T1",
                new Coord(0, 0), new Coord(10, 10),
                0, 600, 2));
        trips.add(new UAMTrip("T2",
                new Coord(1, 1), new Coord(11, 11),
                100, 700, 1));
        trips.add(new UAMTrip("T3",
                new Coord(20, 20), new Coord(30, 30),
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
