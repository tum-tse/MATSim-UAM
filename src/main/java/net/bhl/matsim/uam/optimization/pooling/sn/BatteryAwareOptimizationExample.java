package net.bhl.matsim.uam.optimization.pooling.sn;

import net.bhl.matsim.uam.optimization.pooling.MultiObjectiveNSGAII;
import org.matsim.api.core.v01.Coord;

import java.util.ArrayList;
import java.util.List;

import static net.bhl.matsim.uam.optimization.pooling.sn.EVTOLBatteryManager.DEFAULT_CHARGING_RATE_KWH_PER_SECOND;

/**
 * Complete example demonstrating the battery-aware eVTOL optimization system
 * Shows how battery statistics are automatically calculated and printed with printSummary()
 */
public class BatteryAwareOptimizationExample {

    private static final int VEHICLE_CAPACITY = MultiObjectiveNSGAII.VEHICLE_CAPACITY;
    private static final double VEHICLE_CRUISE_SPEED = MultiObjectiveNSGAII.VEHICLE_CRUISE_SPEED;
    private static final double MAX_DETOUR_RATIO = MultiObjectiveNSGAII.MAX_DETOUR_RATIO;
    private static final int MAX_CONNECTION_TIME_MINUTES = MultiObjectiveNSGAII.MAX_CONNECTION_TIME_MINUTES;

    public static void main(String[] args) {
        System.out.println("=".repeat(80));
        System.out.println("COMPLETE eVTOL BATTERY-AWARE OPTIMIZATION EXAMPLE");
        System.out.println("=".repeat(80));
        System.out.println();

        // Create example trips that demonstrate battery constraints
        List<VehicleTrip> trips = createExampleTrips();

        // Print trip information
        System.out.println("=== INPUT TRIPS ===");
        for (VehicleTrip trip : trips) {
            System.out.printf("Trip %s: Origin(%.0f,%.0f) -> Destination(%.0f,%.0f), " +
                            "Departure: %d s, Passengers: %d\n",
                    trip.getId(),
                    trip.getOrigin().getX(), trip.getOrigin().getY(),
                    trip.getDestination().getX(), trip.getDestination().getY(),
                    trip.getDepartureTime(), trip.getTotalPassengers());
        }
        System.out.println();

        // Create and configure the battery-aware optimizer
        UAMOptimizationController optimizer = new UAMOptimizationController(
                trips,
                MAX_DETOUR_RATIO,
                VEHICLE_CAPACITY,
                MAX_CONNECTION_TIME_MINUTES,
                VEHICLE_CRUISE_SPEED,
                DEFAULT_CHARGING_RATE_KWH_PER_SECOND,
                false
        );

        // Run complete analysis
        OptimizationResult result = optimizer.runCompleteAnalysis();

        System.out.println("\n" + "=".repeat(60));
        System.out.println("TESTING DIFFERENT OPTIMIZATION STRATEGIES");
        System.out.println("=".repeat(60));

        // Compare different optimization strategies
        System.out.println("\n1. Standard Battery-Aware Optimization:");
        OptimizationResult standardResult = optimizer.optimize();
        System.out.println(standardResult.getSummaryString());

        // Create and configure the battery-aware optimizer
        UAMOptimizationController optimizerWithvehicleReuseStrategy = new UAMOptimizationController(
                trips,
                MAX_DETOUR_RATIO,
                VEHICLE_CAPACITY,
                MAX_CONNECTION_TIME_MINUTES,
                VEHICLE_CRUISE_SPEED,
                DEFAULT_CHARGING_RATE_KWH_PER_SECOND,
                true
        );
        System.out.println("\n2. Vehicle Reuse Strategy:");
        OptimizationResult reuseResult = optimizerWithvehicleReuseStrategy.optimize();
        System.out.println(reuseResult.getSummaryString());

        // Show detailed comparison
        System.out.println("\n=== STRATEGY COMPARISON ===");
        if (standardResult.hasBatteryStatistics() && reuseResult.hasBatteryStatistics()) {
            var standardStats = standardResult.getBatteryStatistics();
            var reuseStats = reuseResult.getBatteryStatistics();

            System.out.printf("Standard Strategy: %d vehicles, %.2f kWh total energy\n",
                    standardStats.getTotalVehicles(), standardStats.getTotalEnergyUsed());
            System.out.printf("Reuse Strategy: %d vehicles, %.2f kWh total energy\n",
                    reuseStats.getTotalVehicles(), reuseStats.getTotalEnergyUsed());

            double energySavings = standardStats.getTotalEnergyUsed() - reuseStats.getTotalEnergyUsed();
            int vehicleSavings = standardStats.getTotalVehicles() - reuseStats.getTotalVehicles();

            if (vehicleSavings > 0) {
                System.out.printf("✅ Vehicle reuse saves %d vehicles and %.2f kWh energy\n",
                        vehicleSavings, energySavings);
            } else if (vehicleSavings < 0) {
                System.out.printf("⚠️ Vehicle reuse requires %d more vehicles but saves %.2f kWh energy\n",
                        -vehicleSavings, energySavings);
            } else {
                System.out.println("⚖️ Both strategies use the same number of vehicles");
            }
        }

        // Demonstrate detailed route analysis
        System.out.println("\n=== DETAILED ROUTE ANALYSIS ===");
        standardResult.printDetailedVehicleRoutes();

        // Test battery statistics printing directly
        System.out.println("\n=== BATTERY STATISTICS ONLY ===");
        standardResult.printBatteryStatistics();

        // Demonstrate edge cases
        System.out.println("\n" + "=".repeat(50));
        System.out.println("TESTING EDGE CASES");
        System.out.println("=".repeat(50));

        testHighEnergyTrips();
        testLowEnergyTrips();

        System.out.println("\n" + "=".repeat(80));
        System.out.println("EXAMPLE COMPLETED SUCCESSFULLY!");
        System.out.println("All battery statistics are automatically calculated and printed.");
        System.out.println("=".repeat(80));
    }

    /**
     * Create example trips that demonstrate different battery scenarios
     */
    private static List<VehicleTrip> createExampleTrips() {
        List<VehicleTrip> trips = new ArrayList<>();

        // Morning peak trips (8:00-8:30) - Mixed energy requirements
        trips.add(new VehicleTrip("T1",
                new Coord(10, 10), new Coord(80, 80), // Long distance
                28800, 30000, 3));  // 8:00-8:20 (high energy due to distance + passengers)

        trips.add(new VehicleTrip("T2",
                new Coord(82, 82), new Coord(15, 15), // Return trip, can use same vehicle
                30300, 32100, 2));  // 8:25-8:55 (return journey)

        trips.add(new VehicleTrip("T3",
                new Coord(20, 20), new Coord(25, 25), // Short distance
                28900, 29200, 1));  // 8:01-8:06 (low energy, single passenger)

        // Mid-morning trips (8:30-9:00) - High energy scenarios
        trips.add(new VehicleTrip("T4",
                new Coord(5, 5), new Coord(95, 95), // Very long distance
                30600, 32400, 4));  // 8:30-8:40 (maximum passengers, maximum distance)

        trips.add(new VehicleTrip("T5",
                new Coord(90, 90), new Coord(10, 10), // Another long return trip
                32700, 34500, 3));  // 8:55-9:15 (can reuse previous vehicle after charging)

        // Overlapping trips requiring multiple vehicles
        trips.add(new VehicleTrip("T6",
                new Coord(40, 40), new Coord(60, 60), // Medium distance
                30600, 31200, 2));  // 8:30-8:40 (same time as T4, needs different vehicle)

        trips.add(new VehicleTrip("T7",
                new Coord(30, 30), new Coord(70, 70), // Medium-long distance
                31500, 32700, 3));  // 8:45-9:05 (moderate energy requirements)

        return trips;
    }

    /**
     * Test scenario with trips requiring high energy consumption
     */
    private static void testHighEnergyTrips() {
        System.out.println("\n--- Testing High Energy Consumption Scenario ---");

        List<VehicleTrip> highEnergyTrips = new ArrayList<>();
        // Create trips that will stress the battery system
        highEnergyTrips.add(new VehicleTrip("HE1",
                new Coord(0, 0), new Coord(100, 100), // Very long distance
                30000, 32000, 4)); // Maximum passengers

        highEnergyTrips.add(new VehicleTrip("HE2",
                new Coord(0, 100), new Coord(100, 0), // Cross pattern, long distance
                32500, 34500, 4)); // Maximum passengers

        UAMOptimizationController highEnergyOptimizer = new UAMOptimizationController(
                highEnergyTrips,
                MAX_DETOUR_RATIO,
                VEHICLE_CAPACITY,
                MAX_CONNECTION_TIME_MINUTES,
                VEHICLE_CRUISE_SPEED,
                DEFAULT_CHARGING_RATE_KWH_PER_SECOND,
                false
        );

        OptimizationResult highEnergyResult = highEnergyOptimizer.optimize();
        System.out.println("High energy scenario: " + highEnergyResult.getSummaryString());

        if (highEnergyResult.hasBatteryStatistics()) {
            var stats = highEnergyResult.getBatteryStatistics();
            double maxUtilization = (stats.getMaxEnergyUsedByVehicle() / EVTOLBatteryManager.BATTERY_CAPACITY_KWH) * 100;
            if (maxUtilization > 80) {
                System.out.printf("⚠️  High battery utilization detected: %.1f%%\n", maxUtilization);
            }
        }
    }

    /**
     * Test scenario with trips requiring low energy consumption
     */
    private static void testLowEnergyTrips() {
        System.out.println("\n--- Testing Low Energy Consumption Scenario ---");

        List<VehicleTrip> lowEnergyTrips = new ArrayList<>();
        // Create short trips that allow for good vehicle utilization
        for (int i = 1; i <= 8; i++) {
            lowEnergyTrips.add(new VehicleTrip("LE" + i,
                    new Coord(10 + i * 2, 10 + i * 2),
                    new Coord(15 + i * 2, 15 + i * 2), // Short distances
                    30000 + i * 600, // 10-minute intervals
                    30000 + i * 600 + 300, // 5-minute trips
                    1)); // Single passengers
        }

        UAMOptimizationController lowEnergyOptimizer = new UAMOptimizationController(
                lowEnergyTrips,
                MAX_DETOUR_RATIO,
                VEHICLE_CAPACITY,
                MAX_CONNECTION_TIME_MINUTES,
                VEHICLE_CRUISE_SPEED,
                DEFAULT_CHARGING_RATE_KWH_PER_SECOND,
                false
        );

        OptimizationResult lowEnergyResult = lowEnergyOptimizer.optimize();
        System.out.println("Low energy scenario: " + lowEnergyResult.getSummaryString());

        if (lowEnergyResult.hasBatteryStatistics()) {
            var stats = lowEnergyResult.getBatteryStatistics();
            double avgUtilization = (stats.getAverageEnergyPerVehicle() / EVTOLBatteryManager.BATTERY_CAPACITY_KWH) * 100;
            if (avgUtilization < 20) {
                System.out.printf("ℹ️  Low battery utilization: %.1f%% - good potential for route consolidation\n", avgUtilization);
            }
        }
    }
}