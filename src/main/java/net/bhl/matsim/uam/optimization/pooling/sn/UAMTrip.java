package net.bhl.matsim.uam.optimization.pooling.sn;

// Core classes for UAM shareability network implementation

import org.matsim.api.core.v01.Coord;

import java.util.*;

public class UAMTrip {
    private String id;
    private Coord origin;
    private Coord destination;
    private long departureTime;
    private long arrivalTime;
    private int numPassengers;
    private List<UAMTrip> pooledTrips; // If this is a dummy merged trip

    public UAMTrip(String id, Coord origin, Coord destination,
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

    int getNumPassengers() {
        return numPassengers;
    }

    public long getDepartureTime() {
        return departureTime;
    }

    public long getArrivalTime() {
        return arrivalTime;
    }

    public String getId() {
        return id;
    }

    public Coord getDestination() {
        return destination;
    }

    public Coord getOrigin() {
        return origin;
    }
}
