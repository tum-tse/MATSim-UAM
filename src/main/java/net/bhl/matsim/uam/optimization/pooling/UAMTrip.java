package net.bhl.matsim.uam.optimization.pooling;

// Core classes for UAM shareability network implementation

public class UAMTrip {
    private String id;
    private Location origin;
    private Location destination;
    private long departureTime;
    private long arrivalTime;
    private int numPassengers;
    private List<UAMTrip> pooledTrips; // If this is a dummy merged trip

    public UAMTrip(String id, Location origin, Location destination,
                   long departureTime, long arrivalTime, int numPassengers) {
        this.id = id;
        this.origin = origin;
        this.destination = destination;
        this.departureTime = departureTime;
        this.arrivalTime = arrivalTime;
        this.numPassengers = numPassengers;
        this.pooledTrips = new ArrayList<>();
    }

    // Getters and setters

    public void addPooledTrip(UAMTrip trip) {
        pooledTrips.add(trip);
    }

    public boolean isPooledTrip() {
        return !pooledTrips.isEmpty();
    }

    public List<UAMTrip> getPooledTrips() {
        return Collections.unmodifiableList(pooledTrips);
    }

    public int getTotalPassengers() {
        if (isPooledTrip()) {
            return pooledTrips.stream().mapToInt(UAMTrip::getNumPassengers).sum();
        }
        return numPassengers;
    }
}

public class ShareabilityNetwork {
    private List<UAMTrip> trips;
    private Map<String, Set<String>> adjacencyList;
    private Map<String, UAMVehicle> vehicleAssignments;
    private int maxConnectionTimeMinutes;
    private double flightSpeedMetersPerSecond;

    public ShareabilityNetwork(List<UAMTrip> trips, int maxConnectionTimeMinutes,
                               double flightSpeedMetersPerSecond) {
        this.trips = new ArrayList<>(trips);
        this.adjacencyList = new HashMap<>();
        this.vehicleAssignments = new HashMap<>();
        this.maxConnectionTimeMinutes = maxConnectionTimeMinutes;
        this.flightSpeedMetersPerSecond = flightSpeedMetersPerSecond;

        buildNetwork();
    }

    private void buildNetwork() {
        // Sort trips by departure time
        trips.sort(Comparator.comparingLong(UAMTrip::getDepartureTime));

        // Build adjacency list - connect trips that can be served sequentially
        for (int i = 0; i < trips.size(); i++) {
            UAMTrip t1 = trips.get(i);
            adjacencyList.put(t1.getId(), new HashSet<>());

            for (int j = i + 1; j < trips.size(); j++) {
                UAMTrip t2 = trips.get(j);

                // Check if t2 can be served after t1
                if (canServeSequentially(t1, t2)) {
                    adjacencyList.get(t1.getId()).add(t2.getId());
                }
            }
        }
    }

    private boolean canServeSequentially(UAMTrip t1, UAMTrip t2) {
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

    private double calculateDistance(Location l1, Location l2) {
        // Euclidean distance calculation
        double dx = l1.getX() - l2.getX();
        double dy = l1.getY() - l2.getY();
        return Math.sqrt(dx * dx + dy * dy);
    }

    public List<List<UAMTrip>> findOptimalVehicleAssignments() {
        List<List<UAMTrip>> vehicleRoutes = new ArrayList<>();
        Set<String> unassignedTrips = new HashSet<>(adjacencyList.keySet());

        while (!unassignedTrips.isEmpty()) {
            // Find longest possible path starting from earliest unassigned trip
            List<String> path = findLongestPath(unassignedTrips);

            // Convert path to list of trips
            List<UAMTrip> route = path.stream()
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
        // Find earliest unassigned trip
        String start = availableTrips.stream()
                .min((id1, id2) -> {
                    UAMTrip t1 = trips.stream().filter(t -> t.getId().equals(id1)).findFirst().orElseThrow();
                    UAMTrip t2 = trips.stream().filter(t -> t.getId().equals(id2)).findFirst().orElseThrow();
                    return Long.compare(t1.getDepartureTime(), t2.getDepartureTime());
                })
                .orElseThrow();

        // Use DFS to find longest possible path from start
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

    public double calculateTotalFlightDistance(List<List<UAMTrip>> vehicleRoutes) {
        double totalDistance = 0;

        for (List<UAMTrip> route : vehicleRoutes) {
            // Add flight distances for trips
            for (UAMTrip trip : route) {
                totalDistance += calculateDistance(trip.getOrigin(), trip.getDestination());
            }

            // Add deadheading distances between consecutive trips
            for (int i = 0; i < route.size() - 1; i++) {
                UAMTrip t1 = route.get(i);
                UAMTrip t2 = route.get(i + 1);
                totalDistance += calculateDistance(t1.getDestination(), t2.getOrigin());
            }
        }

        return totalDistance;
    }
}