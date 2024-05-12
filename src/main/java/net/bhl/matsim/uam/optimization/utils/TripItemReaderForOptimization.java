package net.bhl.matsim.uam.optimization.utils;

import net.bhl.matsim.uam.optimization.ModeDecider;
import net.bhl.matsim.uam.optimization.utils.TripItemForOptimization;
import org.matsim.api.core.v01.Coord;
import org.matsim.contrib.util.CSVReaders;
import org.matsim.core.utils.misc.Time;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class TripItemReaderForOptimization {
	public static final double CAR_EMISSION_FACTOR = 0.4;
	public static final double PT_EMISSION_FACTOR = 0.1;
	public static final double CARBON_EQUIVALENCE= 0.1;

	public static List<TripItemForOptimization> getTripItems (String tripsInput) throws IOException{
		List<TripItemForOptimization> trips = new ArrayList<>();
		List<String[]> rows = CSVReaders.readCSV(tripsInput);
		for (String[] row : rows.subList(1, rows.size())) {
			int j = 0;
			TripItemForOptimization trip = new TripItemForOptimization();
			trip.tripID = row[j++];
			trip.personID = Integer.parseInt(row[j++]);
			trip.origin = new Coord(Double.parseDouble(row[j++]), Double.parseDouble(row[j++]));
			trip.destination = new Coord(Double.parseDouble(row[j++]), Double.parseDouble(row[j++]));
			trip.departureTime = Time.parseTime(row[j]);
			trips.add(trip);
		}
		return trips;
	}
	public List<TripItemForOptimization> getTripItemsForOptimization(String tripsInput) throws IOException {
		List<TripItemForOptimization> trips = new ArrayList<>();
		List<String[]> rows = CSVReaders.readCSV(tripsInput);
		for (String[] row : rows.subList(1, rows.size())) {
			int j = 0;
			TripItemForOptimization trip = new TripItemForOptimization();
			trip.tripID = row[j++];
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
			trip.VOT=Double.parseDouble(row[j++])/2080/3600;
			if (row[j].contains("auto"))
				trip.currentMode=0;
			else
				trip.currentMode=1;
			trip.carEmission=CAR_EMISSION_FACTOR*trip.carTripLength/1000;
			trip.ptEmission=CAR_EMISSION_FACTOR*trip.ptTripLength/1000;
		    Double [] probabilities = ModeDecider.calculateModeProbability(-9999, trip.carUtility, trip.ptUtility).toArray(new Double[3]);
			trip.currentGeneralizedCost=trip.carGeneralizedCost*probabilities[1]+trip.ptGeneralizedCost*probabilities[2];
			trips.add(trip);
		}


		return trips;
	}
}



