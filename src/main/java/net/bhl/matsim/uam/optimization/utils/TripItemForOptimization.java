package net.bhl.matsim.uam.optimization.utils;

import net.bhl.matsim.uam.optimization.Vertiport;
import org.matsim.api.core.v01.Coord;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;


public class TripItemForOptimization implements java.io.Serializable{

	public String tripID;
	public int personID;
	public Coord origin;
	public Coord destination;
	public double departureTime;
	public double travelTime;
	public double UAMUtilityFix;
	public double UAMUtilityVar;
	public double uamTravelTime;
	public double VOT;

	public double savedGeneralizedCost;
	public List<Double> tempSavedGeneralizedCosts=new ArrayList<>();
	public HashMap<List<Integer>,Double> tempSavedGeneralizedCostsMap=new HashMap<>();
	public int currentMode; // 0: car, 1: pt
	/*

	public boolean carAvailable;
	public double accessDistance;
	public double egressDistance;
	public double flightDistance;
	public double uamWaitingTime;

	public HashMap<Vertiport, HashMap<String,Double>> accessTimeAndDistanceToAllVertiportCandidates;
	public HashMap<Vertiport, HashMap<String,Double>> egressTimeAndDistanceToAllVertiportCandidates;
	*/
	public String tripPurpose;

    public double carEmission;
	public double ptEmission;
	public double UAMCost;
	public double UAMGeneralizedCost;
	public double carGeneralizedCost;
	public double ptGeneralizedCost;
	public double currentGeneralizedCost;
	public double currentEmission;
	public double carTravelTime;
	public double ptTravelTime;
	public double carTravelCost;
	public double ptTravelCost;
	public double ptWaitingTime;
	public double ptInvehicleTime;
	public double ptTripLength;
	public double carTripLength;
	public boolean isPTAvailable;
	public boolean isUAMAvailable;
	public String description;

	public double distance;
	public double accessTime;
	public double flightTime;
	public double egressTime;
	public double processTime;
	public String accessMode;
	public String egressMode;
	public String originStation;
	public String destinationStation;
	public String mode;
	public double carUtility;
	public double ptUtility;
	public double uamUtility;
	public double carProbability;
	public double ptProbability;
	public double uamProbability;
	public List<Vertiport> originNeighborVertiportCandidates=new ArrayList<>();
	public List<Vertiport> destinationNeighborVertiportCandidates= new ArrayList<>();
	public HashMap<Vertiport,HashMap<String,Double>> originNeighborVertiportCandidatesTimeAndDistance=new HashMap<>();
	public HashMap<Vertiport,HashMap<String,Double>> destinationNeighborVertiportCandidatesTimeAndDistance=new HashMap<>();
	public List<Vertiport> originNeighborVertiports=new ArrayList<>();
	public List<Vertiport> destinationNeighborVertiports= new ArrayList<>();
	public HashMap<Vertiport,HashMap<String,Double>> originNeighborVertiportsTimeAndDistance=new HashMap<>();
	public HashMap<Vertiport,HashMap<String,Double>> destinationNeighborVertiportsTimeAndDistance=new HashMap<>();
	public Vertiport accessVertiport;
	public Vertiport egressVertiport;
}
