package net.bhl.matsim.uam.analysis.traveltimes.utils;

import org.matsim.api.core.v01.Coord;
import org.matsim.contrib.util.CSVReaders;
import org.matsim.core.utils.misc.Time;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class TripItemReader {

	public static List<TripItem> getTripItemsSimple (String tripsInput) throws IOException{
		List<TripItem> trips = new ArrayList<>();
		List<String[]> rows = CSVReaders.readCSV(tripsInput);
		for (String[] row : rows.subList(1, rows.size())) {
			int j = 0;
			TripItem trip = new TripItem();
			trip.tripID = Integer.parseInt(row[j++]);
			trip.personID = Integer.parseInt(row[j++]);
			trip.origin = new Coord(Double.parseDouble(row[j++]), Double.parseDouble(row[j++]));
			trip.destination = new Coord(Double.parseDouble(row[j++]), Double.parseDouble(row[j++]));
			trip.departureTime = Time.parseTime(row[j]);
			trips.add(trip);
		}
		return trips;
	}
	public static List<TripItem> getTripItems(String tripsInput) throws IOException {
		List<TripItem> trips = new ArrayList<>();
		List<String[]> rows = CSVReaders.readCSV(tripsInput);
		for (String[] row : rows.subList(1, rows.size())) {
			int j = 0;
			TripItem trip = new TripItem();
			trip.tripID = Integer.parseInt(row[j++]);
			trip.personID = Integer.parseInt(row[j++]);
			trip.origin = new Coord(Double.parseDouble(row[j++]), Double.parseDouble(row[j++]));
			trip.destination = new Coord(Double.parseDouble(row[j++]), Double.parseDouble(row[j++]));
			trip.departureTime = Time.parseTime(row[j++]);
			trip.ptTravelTime = Double.parseDouble(row[j++]);
			trip.ptTripLength = Double.parseDouble(row[j++]);
			trip.ptInvehicleTime = Double.parseDouble(row[j++]);
			trip.ptWaitingTime = Double.parseDouble(row[j++]);
			trip.carTravelTime = Double.parseDouble(row[j++]);
			trip.carTripLength = Double.parseDouble(row[j++]);
			trip.tripPurpose = row[j++];
			trip.carTravelCost = Double.parseDouble(row[j++]);
			trip.ptTravelCost = Double.parseDouble(row[j++]);
			trip.carUtility = Double.parseDouble(row[j++]);
			trip.ptUtility = Double.parseDouble(row[j++]);
			trip.UAMUtilityFix = Double.parseDouble(row[j++]);
			trip.carGeneralizedCost=Double.parseDouble(row[j++]);
			trip.ptGeneralizedCost=Double.parseDouble(row[j++]);
			trip.VOT=Double.parseDouble(row[j])/2080/3600;
			trips.add(trip);
		}


		return trips;
	}
}



