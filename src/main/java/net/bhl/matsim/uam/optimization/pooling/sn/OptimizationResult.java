package net.bhl.matsim.uam.optimization.pooling.sn;

import java.util.List;

/**
 * Enhanced OptimizationResult that includes battery statistics for eVTOL operations
 * This class maintains backward compatibility with the original OptimizationResult
 */
public class OptimizationResult {
    private final List<List<VehicleTrip>> vehicleRoutes;
    private final int fleetSize;
    private final int vtolOperations;
    private final ShareabilityNetwork.BatteryStatistics batteryStatistics; // New field for battery stats
    private final double chargingRateKwhPerSecond; // Charging rate used for this optimization

    // Enhanced constructor with battery statistics
    public OptimizationResult(List<List<VehicleTrip>> vehicleRoutes,
                              int fleetSize,
                              int vtolOperations,
                              ShareabilityNetwork.BatteryStatistics batteryStatistics,
                              double chargingRateKwhPerSecond) {
        this.vehicleRoutes = vehicleRoutes;
        this.fleetSize = fleetSize;
        this.vtolOperations = vtolOperations;
        this.batteryStatistics = batteryStatistics;
        this.chargingRateKwhPerSecond = chargingRateKwhPerSecond;
    }

    // Backward compatibility constructor (original signature)
    public OptimizationResult(List<List<VehicleTrip>> vehicleRoutes,
                              int fleetSize,
                              int vtolOperations) {
        this.vehicleRoutes = vehicleRoutes;
        this.fleetSize = fleetSize;
        this.vtolOperations = vtolOperations;
        this.batteryStatistics = null; // No battery statistics provided
        this.chargingRateKwhPerSecond = Double.MIN_VALUE;
    }

    // Getters (maintaining original interface)
    public List<List<VehicleTrip>> getVehicleRoutes() {
        return vehicleRoutes;
    }

    public int getFleetSize() {
        return fleetSize;
    }

    public int getVtolOperations() {
        return vtolOperations;
    }

    public ShareabilityNetwork.BatteryStatistics getBatteryStatistics() {
        return batteryStatistics;
    }

    /**
     * Enhanced printSummary that maintains original format but adds battery statistics
     * This replaces the original printSummary method while maintaining backward compatibility
     */
    public void printSummary() {
        System.out.println("Optimization Results:");
        System.out.println("Fleet Size: " + fleetSize);
        System.out.println("Total VTOL Operations: " + vtolOperations);

        // Print battery statistics if available (new functionality)
        if (batteryStatistics != null) {
            System.out.println();
            batteryStatistics.printStatistics();
        }

        // Original vehicle routes printing (maintained from original)
        System.out.println("\nVehicle Routes:");

        for (int i = 0; i < vehicleRoutes.size(); i++) {
            System.out.println("\nVehicle " + (i + 1) + ":");
            List<VehicleTrip> route = vehicleRoutes.get(i);

            if (route.isEmpty()) {
                System.out.println("  No trips assigned");
                continue;
            }

            for (int j = 0; j < route.size(); j++) {
                VehicleTrip trip = route.get(j);

                if (trip.isPooledTrip()) {
                    System.out.println("  Pooled Trip " + trip.getId() + ":");
                    for (VehicleTrip pooledTrip : trip.getPooledTrips()) {
                        System.out.println("    - Trip " + pooledTrip.getId());
                    }
                } else {
                    System.out.println("  Trip " + trip.getId());
                }

                // Add connection time information for battery-aware operations
                if (j < route.size() - 1) {
                    VehicleTrip nextTrip = route.get(j + 1);
                    long connectionTime = nextTrip.getDepartureTime() - trip.getArrivalTime();
                    System.out.printf("    (Connection time to next trip: %d seconds - charging time)\n", connectionTime);
                }
            }
        }
    }

    /**
     * Print only battery statistics (for quick access)
     */
    public void printBatteryStatistics() {
        if (batteryStatistics != null) {
            batteryStatistics.printStatistics();
        } else {
            System.out.println("Battery statistics not available.");
        }
    }

    /**
     * Print detailed vehicle routes information with battery context
     */
    public void printDetailedVehicleRoutes() {
        System.out.println("=== Detailed Vehicle Routes with Battery Information ===");

        for (int i = 0; i < vehicleRoutes.size(); i++) {
            System.out.println("\nVehicle " + (i + 1) + ":");
            List<VehicleTrip> route = vehicleRoutes.get(i);

            if (route.isEmpty()) {
                System.out.println("  No trips assigned");
                continue;
            }

            // Calculate total energy for this route
            double totalRouteEnergy = 0;
            EVTOLBatteryManager tempBattery = new EVTOLBatteryManager(EVTOLBatteryManager.BATTERY_CAPACITY_KWH, chargingRateKwhPerSecond);

            for (int j = 0; j < route.size(); j++) {
                VehicleTrip trip = route.get(j);
                System.out.printf("  Trip %d: %s", j + 1, trip.getId());

                if (trip.isPooledTrip()) {
                    System.out.println(" (Pooled with " + trip.getPooledTrips().size() + " individual trips):");
                    for (VehicleTrip pooledTrip : trip.getPooledTrips()) {
                        System.out.println("    - " + pooledTrip.getId() +
                                " (Passengers: " + pooledTrip.getNumPassengers() + ")");
                    }
                } else {
                    System.out.println(" (Passengers: " + trip.getNumPassengers() + ")");
                }

                // Print trip details
                System.out.printf("    Origin: (%.1f, %.1f), Destination: (%.1f, %.1f)\n",
                        trip.getOrigin().getX(), trip.getOrigin().getY(),
                        trip.getDestination().getX(), trip.getDestination().getY());
                System.out.printf("    Departure: %d s, Arrival: %d s, Duration: %d s\n",
                        trip.getDepartureTime(), trip.getArrivalTime(),
                        trip.getArrivalTime() - trip.getDepartureTime());

                // Calculate and show energy consumption
                double tripEnergy = tempBattery.calculateEnergyConsumption(
                        trip.getOrigin(), trip.getDestination(), trip.getTotalPassengers());
                totalRouteEnergy += tripEnergy;
                System.out.printf("    Energy consumption: %.2f kWh\n", tripEnergy);

                // Show connection time to next trip
                if (j < route.size() - 1) {
                    VehicleTrip nextTrip = route.get(j + 1);
                    long connectionTime = nextTrip.getDepartureTime() - trip.getArrivalTime();
                    double chargedEnergy = connectionTime * tempBattery.getChargingRateKwhPerSecond();
                    System.out.printf("    Connection time to next trip: %d s (%.1f min)\n",
                            connectionTime, connectionTime / 60.0);
                    System.out.printf("    Energy charged during connection: %.2f kWh\n", chargedEnergy);
                }
            }

            System.out.printf("  Total route energy: %.2f kWh (%.1f%% of battery capacity)\n",
                    totalRouteEnergy,
                    (totalRouteEnergy / EVTOLBatteryManager.BATTERY_CAPACITY_KWH) * 100);
        }
    }

    /**
     * Get a summary string of key metrics
     */
    public String getSummaryString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Fleet: ").append(fleetSize).append(" vehicles, ");
        sb.append("VTOL Ops: ").append(vtolOperations);

        if (batteryStatistics != null) {
            sb.append(", Total Energy: ").append(String.format("%.2f kWh",
                    batteryStatistics.getTotalEnergyUsed()));
            sb.append(", Avg Energy/Vehicle: ").append(String.format("%.2f kWh",
                    batteryStatistics.getAverageEnergyPerVehicle()));
            sb.append(", Total Charging: ").append(String.format("%.1f min",
                    batteryStatistics.getTotalChargingTimeMinutes()));
        }

        return sb.toString();
    }

    /**
     * Check if battery statistics are available
     */
    public boolean hasBatteryStatistics() {
        return batteryStatistics != null;
    }

    /**
     * Get detailed performance metrics for analysis
     */
    public PerformanceMetrics getPerformanceMetrics() {
        return new PerformanceMetrics(
                fleetSize,
                vtolOperations,
                batteryStatistics != null ? batteryStatistics.getTotalEnergyUsed() : 0.0,
                batteryStatistics != null ? batteryStatistics.getAverageEnergyPerVehicle() : 0.0,
                batteryStatistics != null ? batteryStatistics.getTotalChargingTimeMinutes() : 0.0,
                batteryStatistics != null ? batteryStatistics.getAverageChargingTimePerVehicle() : 0.0
        );
    }

    /**
     * Helper class for performance metrics
     */
    public static class PerformanceMetrics {
        public final int fleetSize;
        public final int vtolOperations;
        public final double totalEnergyUsed;
        public final double averageEnergyPerVehicle;
        public final double totalChargingTime;
        public final double averageChargingTimePerVehicle;

        public PerformanceMetrics(int fleetSize, int vtolOperations,
                                  double totalEnergyUsed, double averageEnergyPerVehicle,
                                  double totalChargingTime, double averageChargingTimePerVehicle) {
            this.fleetSize = fleetSize;
            this.vtolOperations = vtolOperations;
            this.totalEnergyUsed = totalEnergyUsed;
            this.averageEnergyPerVehicle = averageEnergyPerVehicle;
            this.totalChargingTime = totalChargingTime;
            this.averageChargingTimePerVehicle = averageChargingTimePerVehicle;
        }

        @Override
        public String toString() {
            return String.format("PerformanceMetrics{fleetSize=%d, vtolOps=%d, totalEnergy=%.2f kWh, " +
                            "avgEnergy=%.2f kWh, totalCharging=%.1f min, avgCharging=%.1f min}",
                    fleetSize, vtolOperations, totalEnergyUsed, averageEnergyPerVehicle,
                    totalChargingTime, averageChargingTimePerVehicle);
        }
    }
}