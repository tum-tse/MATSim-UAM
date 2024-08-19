package net.bhl.matsim.uam.optimization.utils;

import net.bhl.matsim.uam.optimization.ModeDecider;
import org.matsim.api.core.v01.Coord;
import org.matsim.contrib.util.CSVReaders;
import org.matsim.core.utils.misc.Time;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class TripItemReaderForOptimization {
	public double CAR_EMISSION_FACTOR ;
	public double PT_EMISSION_FACTOR ;
	public double CARBON_EQUIVALENCE;
	public double car_utility_mean;
	public double car_utility_sigma;
	public double pt_utility_mean;
	public double pt_utility_sigma;
	public long randomSeed;
    public int sampleSize;
	public TripItemReaderForOptimization(ScenarioSpecific scenarioSpecific) {
		this.scenarioSpecific = scenarioSpecific;
		CAR_EMISSION_FACTOR = scenarioSpecific.car_emission_factor;
		PT_EMISSION_FACTOR = scenarioSpecific.pt_emission_factor;
		CARBON_EQUIVALENCE = scenarioSpecific.carbon_equivalent_cost;
		car_utility_mean = scenarioSpecific.car_utility_mean;
		car_utility_sigma = scenarioSpecific.car_utility_sigma;
		pt_utility_mean = scenarioSpecific.pt_utility_mean;
		pt_utility_sigma = scenarioSpecific.pt_utility_sigma;
		randomSeed = scenarioSpecific.random_seed_RT;
		sampleSize = scenarioSpecific.MonteCarlosampleSize;
	}

	public TripItemReaderForOptimization() {
		CARBON_EQUIVALENCE = 0;
		CAR_EMISSION_FACTOR = 0;
		PT_EMISSION_FACTOR = 0;
	}

	public static ScenarioSpecific scenarioSpecific;

	public List<TripItemForOptimization> getTripItems (String tripsInput) throws IOException{
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
			trip.income=Double.parseDouble(row[j++]);
			trip.VOT=trip.income/2080/3600;
			trip.VOT_Less_Than_50km=Double.parseDouble(row[j++]);
			trip.VOT_More_Than_50km=Double.parseDouble(row[j++]);
			if (row[j++].contains("auto"))
				trip.currentMode=0;
			else
				trip.currentMode=1;
			trip.age=Integer.parseInt(row[j++]);
			trip.carOwnership=Integer.parseInt(row[j++]);
			trip.hh_id=row[j++];
			trip.HH_income=Double.parseDouble(row[j++]);
			trip.gender= Integer.parseInt(row[j++]);
			trip.relationship= row[j++];
			trip.occupation=Integer.parseInt(row[j++]);
			trip.driverLicense= Boolean.parseBoolean(row[j++]);
			trip.education=Integer.parseInt(row[j++]);
			trip.autos=Integer.parseInt(row[j]);
			trip.carEmission=CAR_EMISSION_FACTOR*trip.carTripLength/1000;
			trip.ptEmission=CAR_EMISSION_FACTOR*trip.ptTripLength/1000;
			ModeDecider modeDecider = new ModeDecider(-9999, trip.carUtility, trip.ptUtility, car_utility_mean, car_utility_sigma, pt_utility_mean, pt_utility_sigma, new Random(randomSeed));
		    Double [] probabilities = modeDecider.sample(sampleSize);
			trip.carProbabilityBefore=probabilities[1];
			trip.ptProbabilityBefore=probabilities[2];
			trip.currentGeneralizedCost=trip.carGeneralizedCost*probabilities[1]+trip.ptGeneralizedCost*probabilities[2];
			trip.currentEmission=trip.carEmission*probabilities[1]+trip.ptEmission*probabilities[2];
			trip.currentTravelTime=trip.carTravelTime*probabilities[1]+trip.ptTravelTime*probabilities[2];

			trips.add(trip);
		}


		return trips;
	}
}



