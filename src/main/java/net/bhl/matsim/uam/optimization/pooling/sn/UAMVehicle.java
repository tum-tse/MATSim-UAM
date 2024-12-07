package net.bhl.matsim.uam.optimization.pooling.sn;

import org.matsim.api.core.v01.Coord;

import java.util.ArrayList;
import java.util.List;

public class UAMVehicle {
    private String id;
    private Coord currentLocation;
    private List<UAMTrip> assignedTrips;
    private long nextAvailableTime;
    private final int capacity;

    public UAMVehicle(String id, Coord initialLocation, int capacity) {
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
