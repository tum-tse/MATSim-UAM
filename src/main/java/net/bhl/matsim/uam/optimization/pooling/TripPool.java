package net.bhl.matsim.uam.optimization.pooling;

public class TripPool {
    private List<UAMTrip> trips;
    private final int maxPassengers;
    private final double maxDetourRatio;
    private Location pooledOrigin;
    private Location pooledDestination;
    private long pooledDepartureTime;
    private long pooledArrivalTime;

    public TripPool(int maxPassengers, double maxDetourRatio) {
        this.trips = new ArrayList<>();
        this.maxPassengers = maxPassengers;
        this.maxDetourRatio = maxDetourRatio;
    }

    public boolean canAddTrip(UAMTrip trip) {
        // Check total passengers
        int totalPassengers = getTotalPassengers() + trip.getNumPassengers();
        if (totalPassengers > maxPassengers) {
            return false;
        }

        if (trips.isEmpty()) {
            return true;
        }

        // Calculate new pooled locations and times
        Location newOrigin = calculatePooledLocation(
                Stream.concat(trips.stream(), Stream.of(trip))
                        .map(UAMTrip::getOrigin)
                        .collect(Collectors.toList())
        );

        Location newDestination = calculatePooledLocation(
                Stream.concat(trips.stream(), Stream.of(trip))
                        .map(UAMTrip::getDestination)
                        .collect(Collectors.toList())
        );

        // Check detour ratio for all trips including new one
        for (UAMTrip t : trips) {
            double originalDistance = calculateDistance(t.getOrigin(), t.getDestination());
            double newDistance = calculateDistance(newOrigin, t.getOrigin()) +
                    calculateDistance(newOrigin, newDestination) +
                    calculateDistance(newDestination, t.getDestination());

            if (newDistance / originalDistance > (1 + maxDetourRatio)) {
                return false;
            }
        }

        // Check detour ratio for new trip
        double originalDistance = calculateDistance(trip.getOrigin(), trip.getDestination());
        double newDistance = calculateDistance(newOrigin, trip.getOrigin()) +
                calculateDistance(newOrigin, newDestination) +
                calculateDistance(newDestination, trip.getDestination());

        return newDistance / originalDistance <= (1 + maxDetourRatio);
    }

    public void addTrip(UAMTrip trip) {
        if (!canAddTrip(trip)) {
            throw new IllegalArgumentException("Cannot add trip to pool");
        }

        trips.add(trip);
        updatePooledProperties();
    }

    private void updatePooledProperties() {
        // Update pooled locations
        pooledOrigin = calculatePooledLocation(
                trips.stream().map(UAMTrip::getOrigin).collect(Collectors.toList())
        );

        pooledDestination = calculatePooledLocation(
                trips.stream().map(UAMTrip::getDestination).collect(Collectors.toList())
        );

        // Update pooled times
        pooledDepartureTime = trips.stream()
                .mapToLong(UAMTrip::getDepartureTime)
                .min()
                .orElseThrow();

        pooledArrivalTime = trips.stream()
                .mapToLong(UAMTrip::getArrivalTime)
                .max()
                .orElseThrow();
    }

    private Location calculatePooledLocation(List<Location> locations) {
        // Calculate centroid
        double sumX = locations.stream().mapToDouble(Location::getX).sum();
        double sumY = locations.stream().mapToDouble(Location::getY).sum();
        return new Location(sumX / locations.size(), sumY / locations.size());
    }

    public UAMTrip toPooledTrip() {
        if (trips.isEmpty()) {
            throw new IllegalStateException("No trips in pool");
        }

        UAMTrip pooledTrip = new UAMTrip(
                "POOL_" + trips.stream().map(UAMTrip::getId).collect(Collectors.joining("_")),
                pooledOrigin,
                pooledDestination,
                pooledDepartureTime,
                pooledArrivalTime,
                getTotalPassengers()
        );

        // Add individual trips to pooled trip
        trips.forEach(pooledTrip::addPooledTrip);

        return pooledTrip;
    }

    private int getTotalPassengers() {
        return trips.stream().mapToInt(UAMTrip::getNumPassengers).sum();
    }

    private double calculateDistance(Location l1, Location l2) {
        double dx = l1.getX() - l2.getX();
        double dy = l1.getY() - l2.getY();
        return Math.sqrt(dx * dx + dy * dy);
    }
}

public class UAMVehicle {
    private String id;
    private Location currentLocation;
    private List<UAMTrip> assignedTrips;
    private long nextAvailableTime;
    private final int capacity;

    public UAMVehicle(String id, Location initialLocation, int capacity) {
        this.id = id;
        this.currentLocation = initialLocation;
        this.capacity = capacity;
        this.assignedTrips = new ArrayList<>();
        this.nextAvailableTime = 0;
    }

    // Getters and setters

    public boolean canAssignTrip(UAMTrip trip) {
        return trip.getTotalPassengers() <= capacity;
    }

    public void assignTrip(UAMTrip trip) {
        if (!canAssignTrip(trip)) {
            throw new IllegalArgumentException("Trip exceeds vehicle capacity");
        }

        assignedTrips.add(trip);
        // Update current location and next available time
        currentLocation = trip.getDestination();
        nextAvailableTime = trip.getArrivalTime();
    }
}
