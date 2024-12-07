package net.bhl.matsim.uam.optimization.pooling.sn;

import org.matsim.api.core.v01.Coord;

import java.util.ArrayList;
import java.util.List;

// Example usage
public class RunShareabilityNetwork {
    public static void main(String[] args) {
        // Create sample trips with overlapping time windows
        List<UAMTrip> trips = new ArrayList<>();

        // Group 1: Early morning trips (7:00-7:30 AM)
        trips.add(new UAMTrip("T1",
                new Coord(0, 0), new Coord(10, 10),
                25200, 25800, 2));  // 7:00-7:10
        trips.add(new UAMTrip("T2",
                new Coord(1, 1), new Coord(11, 11),
                25300, 25900, 1));  // 7:01-7:11
        trips.add(new UAMTrip("T3",
                new Coord(2, 2), new Coord(12, 12),
                25400, 26000, 2));  // 7:02-7:12

        // Group 2: Concurrent trips requiring multiple vehicles (7:15-7:45)
        trips.add(new UAMTrip("T4",
                new Coord(20, 20), new Coord(30, 30),
                26100, 26700, 2));  // 7:15-7:25
        trips.add(new UAMTrip("T5",
                new Coord(21, 21), new Coord(31, 31),
                26100, 26700, 2));  // 7:15-7:25
        trips.add(new UAMTrip("T6",
                new Coord(22, 22), new Coord(32, 32),
                26100, 26700, 1));  // 7:15-7:25

        // Group 3: Potential pooling candidates (7:30-8:00)
        trips.add(new UAMTrip("T7",
                new Coord(5, 5), new Coord(15, 15),
                27000, 27600, 1));  // 7:30-7:40
        trips.add(new UAMTrip("T8",
                new Coord(6, 6), new Coord(16, 16),
                27060, 27660, 2));  // 7:31-7:41
        trips.add(new UAMTrip("T9",
                new Coord(7, 7), new Coord(17, 17),
                27120, 27720, 1));  // 7:32-7:42

        // Group 4: Mid-morning trips (8:00-8:30)
        trips.add(new UAMTrip("T10",
                new Coord(40, 40), new Coord(50, 50),
                28800, 29400, 2));  // 8:00-8:10
        trips.add(new UAMTrip("T11",
                new Coord(41, 41), new Coord(51, 51),
                28860, 29460, 2));  // 8:01-8:11
        trips.add(new UAMTrip("T12",
                new Coord(42, 42), new Coord(52, 52),
                28920, 29520, 1));  // 8:02-8:12

        // Group 5: Overlapping peaks (8:15-8:45)
        trips.add(new UAMTrip("T13",
                new Coord(25, 25), new Coord(35, 35),
                29700, 30300, 2));  // 8:15-8:25
        trips.add(new UAMTrip("T14",
                new Coord(26, 26), new Coord(36, 36),
                29700, 30300, 1));  // 8:15-8:25
        trips.add(new UAMTrip("T15",
                new Coord(27, 27), new Coord(37, 37),
                29700, 30300, 2));  // 8:15-8:25

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
