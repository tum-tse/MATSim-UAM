package net.bhl.matsim.uam.optimization.utils;

import net.bhl.matsim.uam.optimization.Vertiport;
import org.matsim.api.core.v01.Coord;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;


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
	public double VOT_Less_Than_50km;
	public double VOT_More_Than_50km;

	public double savedGeneralizedCost;
	public List<Double> tempSavedGeneralizedCosts=new ArrayList<>();
	public ConcurrentHashMap<List<Integer>,Double> tempSavedGeneralizedCostsMap=new ConcurrentHashMap<>();
	public double savedEmission;
	public List<Double> tempSavedEmission=new ArrayList<>();
	public ConcurrentHashMap<List<Integer>,Double> tempSavedEmissionMap=new ConcurrentHashMap<>();
	public double savedTravelTime;
	public List<Double> tempSavedTravelTime=new ArrayList<>();
	public ConcurrentHashMap<List<Integer>,Double> tempSavedTravelTimeMap=new ConcurrentHashMap<>();
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
	public double uamEmission;
	public double UAMCost;
	public double UAMGeneralizedCost;
	public double carGeneralizedCost;
	public double ptGeneralizedCost;
	public double currentGeneralizedCost;
	public double currentEmission;
	public double currentTravelTime;
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
	public double tempAccessTime;
	public double flightTime;
	public double tempFlightTime;
	public double egressTime;
	public double tempEgressTime;
	public double processTime;
	public String accessMode;
	public String tempAccessMode;
	public String egressMode;
	public String tempEgressMode;
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
	public Vertiport tempAccessVertiport;
	public Vertiport tempEgressVertiport;
	public double HH_income;
	public int age;
	public double carProbabilityBefore;
	public double ptProbabilityBefore;
	public double tempCarProbability;
	public double tempPTProbability;
	public double tempUAMProbability;
	public boolean tempIsUAMAvailable;
	public double tempUamTravelTime;
	public double tempUAMCost;
	public double tempUamEmission;
	public double tempUAMGeneralizedCost;
	public double tempUAMUtilityVar;
	public double tempUamUtility;
    public int carOwnership;
	public String hh_id;
	public int gender;
	public String relationship;
	public int occupation;
	public int education;
	public boolean driverLicense;
	public double income;
	public int autos;

}
