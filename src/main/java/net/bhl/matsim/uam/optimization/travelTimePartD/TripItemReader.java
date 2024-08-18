package net.bhl.matsim.uam.optimization.travelTimePartD;

import org.matsim.api.core.v01.Coord;
import org.matsim.contrib.util.CSVReaders;
import org.matsim.core.utils.misc.Time;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class TripItemReader {
    public static List<TripItem> getTripItems(String tripsInput) throws IOException {
        List<TripItem> trips = new ArrayList<>();
        List<String[]> rows = CSVReaders.readCSV(tripsInput);
        for (String[] row : rows.subList(1, rows.size())) {
            int j = 0;
            TripItem trip = new TripItem();
            trip.tripID = row[0];
            trip.origin = new Coord(Double.parseDouble(row[2]), Double.parseDouble(row[3]));
            trip.destination = new Coord(Double.parseDouble(row[5]), Double.parseDouble(row[6]));
            trip.departureTime = Time.parseTime(row[15])*60; // Convert to seconds
            trips.add(trip);
        }
        return trips;
    }
}
