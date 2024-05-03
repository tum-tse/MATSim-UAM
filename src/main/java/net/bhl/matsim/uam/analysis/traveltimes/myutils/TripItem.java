package net.bhl.matsim.uam.analysis.traveltimes.myutils;

import org.matsim.api.core.v01.Coord;


public class TripItem {
    public String tripId;
    public Coord origin;
    public Coord destination;
    public double departureTime;
    public double travelTime;
    public double ptWaitingTime;
    public double ptInVehicleTime;
    public int ptTransfers;
    public String description;

    public double distance;
    public double flightDistance;
    public double accessTime;
    public double accessDistance;
    public double flightTime;
    public double egressTime;
    public double egressDistance;
    public double processTime;
    public String accessMode;
    public String egressMode;
    public String originStation;
    public String destinationStation;
}