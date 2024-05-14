package net.bhl.matsim.uam.optimization;

import org.matsim.api.core.v01.Coord;

import java.util.HashMap;
import java.util.List;

public class Vertiport implements java.io.Serializable{
    public int ID;
    public Coord coord;
    public double constructionCost;
    public double capacity;
    // double array of saturation rates for each hour of the day with the length of 24
    public HashMap<Integer,Double> saturationRates=new HashMap<>();
    public double maxSaturationRate;
    public List<Vertiport> neighbors;
}
