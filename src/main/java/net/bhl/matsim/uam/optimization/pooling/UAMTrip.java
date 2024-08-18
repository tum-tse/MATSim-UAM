package net.bhl.matsim.uam.optimization.pooling;

import net.bhl.matsim.uam.infrastructure.UAMStation;

public static class UAMTrip {
    private final String tripId;
    private final double originX, originY, destX, destY, departureTime, flightDistance;
    private UAMStation origStation, destStation; // Changed to Integer to handle null values
    private final String purpose, income;
    private double accessTimeToPooledStation; // Time for access to the station

    public UAMTrip(String tripId, double originX, double originY, double destX, double destY, double departureTime, double flightDistance, UAMStation origStation, UAMStation destStation, String purpose, String income) {
        this.tripId = tripId;
        this.originX = originX;
        this.originY = originY;
        this.destX = destX;
        this.destY = destY;
        this.departureTime = departureTime;
        this.flightDistance = flightDistance;
        this.origStation = origStation;
        this.destStation = destStation;
        this.purpose = purpose;
        this.income = income;
    }

    // TODO: Use MATSim to calculate the routes and travel distances
    public double calculateAccessTeleportationDistance(UAMStation station) {
        return Math.sqrt(Math.pow(originX - station.getLocationLink().getCoord().getX(), 2) + Math.pow(originY - station.getLocationLink().getCoord().getY(), 2));
    }
    // TODO: Use MATSim to calculate the routes and travel times
    public double calculateAccessTeleportationTime(UAMStation station) {
        double distance = calculateAccessTeleportationDistance(station);
        accessTimeToPooledStation = distance / TELEPORTATION_SPEED;
        return accessTimeToPooledStation;
    }
    // TODO: Use MATSim to calculate the routes and travel distances
    public double calculateEgressTeleportationDistance(UAMStation station) {
        return Math.sqrt(Math.pow(destX - station.getLocationLink().getCoord().getX(), 2) + Math.pow(destY - station.getLocationLink().getCoord().getY(), 2));
    }
    // TODO: Use MATSim to calculate the routes and travel times
    public double calculateEgressTeleportationTime(UAMStation station) {
        double distance = calculateEgressTeleportationDistance(station);
        return distance / TELEPORTATION_SPEED;
    }
    // TODO: Use MATSim to calculate the routes and travel distances
    public double calculateFlightDistance(UAMStation originStation, UAMStation destStation) {
        return Math.sqrt(Math.pow(originStation.getLocationLink().getCoord().getX() - destStation.getLocationLink().getCoord().getX(), 2) + Math.pow(originStation.getLocationLink().getCoord().getY() - destStation.getLocationLink().getCoord().getY(), 2));
    }

    // setOriginStation
    public void setOriginStation(UAMStation station) {
        this.origStation = station;
    }
    // setDestinationStation
    public void setDestinationStation(UAMStation station) {
        this.destStation = station;
    }

    // getDepartureTime
    public double getDepartureTime() {
        return departureTime;
    }
    //getOriginX
    public double getOriginX() {
        return originX;
    }
    //getOriginY
    public double getOriginY() {
        return originY;
    }
    //getOriginStation
    public UAMStation getOriginStation() {
        return origStation;
    }
    //getDestinationStation
    public UAMStation getDestinationStation() {
        return destStation;
    }
    //getTripId
    public String getTripId() {
        return tripId;
    }
    //getFlightDistance
    public double getFlightDistance() {
        return flightDistance;
    }
}
