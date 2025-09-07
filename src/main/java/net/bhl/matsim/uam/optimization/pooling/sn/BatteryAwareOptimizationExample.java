package net.bhl.matsim.uam.optimization.pooling.sn;

import net.bhl.matsim.uam.optimization.pooling.MultiObjectiveNSGAII;
import org.matsim.api.core.v01.Coord;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static net.bhl.matsim.uam.optimization.pooling.sn.EVTOLBatteryManager.DEFAULT_CHARGING_RATE_KWH_PER_SECOND;

/**
 * Simplified battery-aware eVTOL optimization system for extreme stress testing only
 */
public class BatteryAwareOptimizationExample {

    // Essential constants for extreme stress test
    private static final int VEHICLE_CAPACITY = MultiObjectiveNSGAII.VEHICLE_CAPACITY;
    private static final double VEHICLE_CRUISE_SPEED = MultiObjectiveNSGAII.VEHICLE_CRUISE_SPEED;
    private static final double MAX_DETOUR_RATIO = MultiObjectiveNSGAII.MAX_DETOUR_RATIO;
    private static final int MAX_CONNECTION_TIME_MINUTES = MultiObjectiveNSGAII.MAX_CONNECTION_TIME_MINUTES;

    // Only the charging rates used in extreme stress test
    private static final double[] CRITICAL_CHARGING_RATES = {
            0 / 60.0,       // No charging (0C)
            4.17 / 60.0,    // Slow charging (2C)
            6.24 / 60.0,    // Default charging rate (3C)
            8.34 / 60.0,    // Fast charging (4C)
            16.68 / 60.0    // Maximum charging (8C)
    };

    public static void main(String[] args) {
        System.out.println("=".repeat(80));
        System.out.println("eVTOL EXTREME STRESS TEST");
        System.out.println("=".repeat(80));

        runExtremeStressTest();

        System.out.println("\n" + "=".repeat(80));
        System.out.println("EXTREME STRESS TEST COMPLETED!");
        System.out.println("=".repeat(80));
    }

    /**
     * Run extreme stress test with maximum trip density
     */
    private static void runExtremeStressTest() {
        System.out.println("=".repeat(80));
        System.out.println("EXTREME STRESS TEST: MAXIMUM SYSTEM LOAD");
        System.out.println("=".repeat(80));

        // Create extreme scenario: 100 trips in 60-minute window
        List<VehicleTrip> extremeTrips = createExtremeStressTestTrips();

        System.out.printf("Generated %d trips in 60-minute peak demand window\n", extremeTrips.size());
        printTripDistribution(extremeTrips);

        // Test with different charging strategies
        String[] criticalRateNames = {
                "No charging", "Slow charging", "Default Rate", "Fast charging", "Maximum charging"
        };

        System.out.println("\n--- Critical Charging Rate Comparison ---");
        for (int i = 0; i < CRITICAL_CHARGING_RATES.length; i++) {
            try {
                long startTime = System.currentTimeMillis();

                UAMOptimizationController optimizer = new UAMOptimizationController(
                        extremeTrips,
                        MAX_DETOUR_RATIO,
                        VEHICLE_CAPACITY,
                        MAX_CONNECTION_TIME_MINUTES,
                        VEHICLE_CRUISE_SPEED,
                        CRITICAL_CHARGING_RATES[i],
                        false
                );

                OptimizationResult result = optimizer.optimize();
                long duration = System.currentTimeMillis() - startTime;

                System.out.printf("%-15s: %s (%.1fs)\n",
                        criticalRateNames[i], result.getSummaryString(), duration / 1000.0);

                if (result.hasBatteryStatistics()) {
                    result.printBatteryStatistics();
                }

            } catch (Exception e) {
                System.out.printf("%-15s: FAILED - %s\n", criticalRateNames[i], e.getMessage());
            }
            System.out.println();
        }
    }

    /**
     * Create extreme stress test with very high trip density
     */
    private static List<VehicleTrip> createExtremeStressTestTrips() {
        List<VehicleTrip> trips = new ArrayList<>();
        Random random = new Random(12345); // Fixed seed for reproducibility

        int baseTime = 28800; // 8:00 AM
        int timeWindow = 3600;  // 60 minutes (extreme density)
        int tripCount = 100;

        for (int i = 1; i <= tripCount; i++) {
            // Create high-density urban scenario with shorter distances
            double centerX = 75000; // City center
            double centerY = 75000;
            double urbanRadius = 10000; // 30km radius urban area

            // Origins clustered around urban centers
            double originAngle = random.nextDouble() * 2 * Math.PI;
            double originRadius = random.nextGaussian() * urbanRadius * 0.3 + urbanRadius * 0.5;
            Coord origin = new Coord(
                    centerX + originRadius * Math.cos(originAngle),
                    centerY + originRadius * Math.sin(originAngle)
            );

            // Destinations with preference for city center and airports
            double destAngle = random.nextDouble() * 2 * Math.PI;
            double destRadius = random.nextGaussian() * urbanRadius * 0.4 + urbanRadius * 0.6;
            Coord destination = new Coord(
                    centerX + destRadius * Math.cos(destAngle),
                    centerY + destRadius * Math.sin(destAngle)
            );

            // Peak hour departure times (concentrated)
            int departureTime = baseTime + random.nextInt(timeWindow);

            // Calculate trip duration based on distance
            double distance = Math.sqrt(Math.pow(destination.getX() - origin.getX(), 2) +
                    Math.pow(destination.getY() - origin.getY(), 2));
            int tripDuration = (int) (distance / VEHICLE_CRUISE_SPEED) + 300 + random.nextInt(300);

            // Higher passenger loads in extreme scenario
            int passengers = random.nextDouble() < 0.7 ? VEHICLE_CAPACITY - 1 + random.nextInt(2) : 1;

            trips.add(new VehicleTrip("EX" + i,
                    origin,
                    destination,
                    departureTime,
                    departureTime + tripDuration,
                    passengers));
        }

        return trips;
    }

    /**
     * Print distribution of trips by time and distance
     */
    private static void printTripDistribution(List<VehicleTrip> trips) {
        int shortTrips = 0, mediumTrips = 0, longTrips = 0;
        int singlePassenger = 0, multiplePassengers = 0;

        for (VehicleTrip trip : trips) {
            double distance = Math.sqrt(
                    Math.pow(trip.getDestination().getX() - trip.getOrigin().getX(), 2) +
                            Math.pow(trip.getDestination().getY() - trip.getOrigin().getY(), 2)
            );

            if (distance < 20000) shortTrips++;
            else if (distance < 50000) mediumTrips++;
            else longTrips++;

            if (trip.getTotalPassengers() == 1) singlePassenger++;
            else multiplePassengers++;
        }

        System.out.printf("Trip Distribution: %d short (<20km), %d medium (20-50km), %d long (>50km)\n",
                shortTrips, mediumTrips, longTrips);
        System.out.printf("Passenger Load: %d single-passenger, %d multi-passenger\n",
                singlePassenger, multiplePassengers);
    }
}