package net.bhl.matsim.uam.optimization;

import org.matsim.api.core.v01.Coord;

import java.util.HashMap;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;

public class Vertiport implements java.io.Serializable{
    public int ID;
    public Coord coord;
    public double constructionCost;
    public int totalCapacity;
    // double array of saturation rates for each hour of the day with the length of 24
    public HashMap<Integer,Double> saturationRates=new HashMap<>();
    public ConcurrentHashMap<Integer,Double> concurrentSaturationRates=new ConcurrentHashMap<>();
    public ConcurrentHashMap<Integer,Double> tempConcurrentSaturationRates=new ConcurrentHashMap<>();
    public ConcurrentHashMap<Integer,Double> concurrentAccessDemand=new ConcurrentHashMap<>();
    public ConcurrentHashMap<Integer,Double> tempConcurrentAccessDemand=new ConcurrentHashMap<>();
    public ConcurrentHashMap<Integer,Double> concurrentEgressDemand=new ConcurrentHashMap<>();
    public ConcurrentHashMap<Integer,Double> tempConcurrentEgressDemand=new ConcurrentHashMap<>();
    public double tempTotalDemand;
    public double totalDemand;
    public double tempWaitingAreaCapacity;
    public double waitingAreaCapacity;
    public double tempAvailableCapacity;
    public double availableCapacity;
    public double tempMax15MinDemand;
    public double max15MinDemand;
    public double tempMaxSaturationRate;
    public double maxSaturationRate;
    public List<Vertiport> neighbors;
}
