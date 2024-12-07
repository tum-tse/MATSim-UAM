package net.bhl.matsim.uam.optimization.pooling.sn;

import org.matsim.api.core.v01.Coord;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Core classes for UAM shareability network implementation
 * Represents a shareability network for UAM trips.
 * The network is represented as an adjacency list, where each trip is a node and an edge between two nodes indicates that the corresponding trips can be served sequentially.
 * The network is built based on the trips' departure and arrival times, and the maximum connection time.
 * The network can be used to find optimal vehicle assignments for the trips.
 */
public class ShareabilityNetwork {
    private List<PassengerTrip> trips;
    private Map<String, Set<String>> adjacencyList;
    //private Map<String, UAMVehicle> vehicleAssignments;
    private int maxConnectionTimeMinutes;
    private double flightSpeedMetersPerSecond;

    public ShareabilityNetwork(List<PassengerTrip> trips, int maxConnectionTimeMinutes,
                               double flightSpeedMetersPerSecond) {
        this.trips = new ArrayList<>(trips);
        this.adjacencyList = new HashMap<>();
        //this.vehicleAssignments = new HashMap<>();
        this.maxConnectionTimeMinutes = maxConnectionTimeMinutes;
        this.flightSpeedMetersPerSecond = flightSpeedMetersPerSecond;

        buildNetwork();
    }

    private void buildNetwork() {
        // Sort trips by departure time
        trips.sort(Comparator.comparingLong(PassengerTrip::getDepartureTime));

        // Build adjacency list - connect trips that can be served sequentially
        for (int i = 0; i < trips.size(); i++) {
            PassengerTrip t1 = trips.get(i);
            adjacencyList.put(t1.getId(), new HashSet<>());

            for (int j = i + 1; j < trips.size(); j++) {
                PassengerTrip t2 = trips.get(j);

                // Check if t2 can be served after t1
                if (canServeSequentially(t1, t2)) {
                    adjacencyList.get(t1.getId()).add(t2.getId());
                }
            }
        }
    }

    private boolean canServeSequentially(PassengerTrip t1, PassengerTrip t2) {
        // Calculate flight time between t1's destination and t2's origin
        double distanceMeters = calculateDistance(t1.getDestination(), t2.getOrigin());
        long flightTimeSeconds = (long) (distanceMeters / flightSpeedMetersPerSecond);

        // Calculate earliest possible arrival time at t2's origin
        long earliestArrival = t1.getArrivalTime() + flightTimeSeconds;

        // Check if connection time is within limit and vehicle can arrive before t2's departure
        long connectionTime = (t2.getDepartureTime() - t1.getArrivalTime()) / 60; // Convert to minutes
        return connectionTime <= maxConnectionTimeMinutes &&
                earliestArrival <= t2.getDepartureTime();
    }

    private double calculateDistance(Coord l1, Coord l2) {
        // Euclidean distance calculation
        double dx = l1.getX() - l2.getX();
        double dy = l1.getY() - l2.getY();
        return Math.sqrt(dx * dx + dy * dy);
    }

    public List<List<PassengerTrip>> findOptimalVehicleAssignments() {
        List<List<PassengerTrip>> vehicleRoutes = new ArrayList<>();
        Set<String> unassignedTrips = new HashSet<>(adjacencyList.keySet());

        while (!unassignedTrips.isEmpty()) {
            // Find the longest possible path starting from the earliest unassigned trip
            List<String> path = findLongestPath(unassignedTrips);

            // Convert path to list of trips
            List<PassengerTrip> route = path.stream()
                    .map(id -> trips.stream()
                            .filter(t -> t.getId().equals(id))
                            .findFirst()
                            .orElseThrow())
                    .collect(Collectors.toList());

            vehicleRoutes.add(route);
            unassignedTrips.removeAll(path);
        }

        return vehicleRoutes;
    }

    private List<String> findLongestPath(Set<String> availableTrips) {
        // Find the earliest unassigned trip
        String start = availableTrips.stream()
                .min((id1, id2) -> {
                    PassengerTrip t1 = trips.stream().filter(t -> t.getId().equals(id1)).findFirst().orElseThrow();
                    PassengerTrip t2 = trips.stream().filter(t -> t.getId().equals(id2)).findFirst().orElseThrow();
                    return Long.compare(t1.getDepartureTime(), t2.getDepartureTime());
                })
                .orElseThrow();

        // Use DFS to find the longest possible path from start
        List<String> longestPath = new ArrayList<>();
        findLongestPathDFS(start, availableTrips, new ArrayList<>(), longestPath);
        return longestPath;
    }

    private void findLongestPathDFS(String current, Set<String> availableTrips,
                                    List<String> currentPath, List<String> longestPath) {
        currentPath.add(current);

        if (currentPath.size() > longestPath.size()) {
            longestPath.clear();
            longestPath.addAll(currentPath);
        }

        // Explore all possible next trips
        for (String next : adjacencyList.get(current)) {
            if (availableTrips.contains(next)) {
                findLongestPathDFS(next, availableTrips, currentPath, longestPath);
            }
        }

        currentPath.remove(currentPath.size() - 1);
    }

    public double calculateTotalFlightDistance(List<List<PassengerTrip>> vehicleRoutes) {
        double totalDistance = 0;

        for (List<PassengerTrip> route : vehicleRoutes) {
            // Add flight distances for trips
            for (PassengerTrip trip : route) {
                totalDistance += calculateDistance(trip.getOrigin(), trip.getDestination());
            }

            // Add deadheading distances between consecutive trips
            for (int i = 0; i < route.size() - 1; i++) {
                PassengerTrip t1 = route.get(i);
                PassengerTrip t2 = route.get(i + 1);
                totalDistance += calculateDistance(t1.getDestination(), t2.getOrigin());
            }
        }

        return totalDistance;
    }
    public double calculateTotalDeadheadingFlightDistance(List<List<PassengerTrip>> vehicleRoutes) {
        double totalDeadheadingDistance = 0;

        for (List<PassengerTrip> route : vehicleRoutes) {
            // Add deadheading distances between consecutive trips
            for (int i = 0; i < route.size() - 1; i++) {
                PassengerTrip t1 = route.get(i);
                PassengerTrip t2 = route.get(i + 1);
                totalDeadheadingDistance += calculateDistance(t1.getDestination(), t2.getOrigin());
            }
        }

        return totalDeadheadingDistance;
    }
}
