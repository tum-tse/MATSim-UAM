package net.bhl.matsim.uam.analysis.traveltimes.traveldisutility;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.matsim.api.core.v01.Id;
import org.matsim.api.core.v01.network.Link;
import org.matsim.api.core.v01.population.Person;
import org.matsim.core.gbl.Gbl;
import org.matsim.core.router.util.LinkToLinkTravelTime;
import org.matsim.core.router.util.TravelDisutility;
import org.matsim.core.router.util.TravelTime;
import org.matsim.vehicles.Vehicle;

import org.matsim.contrib.util.CSVReaders;

import com.google.common.collect.Iterables;

/**<p>
 * CostCalculator and TravelTimeCalculator for Links based on freespeed on links and
 * distance costs if set.  It sets the <em> function </em> that is to be used with calls
 * <tt>getLinkTravelTime( link, time)</tt> and <tt>getLinkTravelCost( link, time )</tt>.
 * </p><p>
 * The unit of "cost" is defined by the input: if the marginal utilities are given in "utils per second", then
 * cost is in "utils"; if the marginal utilities are given in "euros per second", then cost is in "euros".
 * When the CharyparNagelScoringFunction is used, the values come from the config file, where one is also free to
 * interpret the units.
 * </p>
 *
 * @author haowuintub (Hao Wu) based on FreespeedTravelTimeAndDisutility
 * Note: the travel time used for calculating the time point of querying the disutility is based on the free speed travel time!
 */
public class NoiseBasedTravelDisutility implements TravelDisutility, TravelTime, LinkToLinkTravelTime {
    private static final Logger log = LogManager.getLogger(NoiseBasedTravelDisutility.class);

    private final double travelCostFactor;
    private final double marginalUtlOfDistance;
    //private static int wrnCnt = 0 ;
    private final String noiseEmissionResults;
    private final double timeBin; // unit in seconds. e.g. 3600.0
    private final int simulationEndHour; // simulation hour counts start from 1. The simulation end hour could be for example 34
    private Map<Double, Map<Id<Link>, Double>> noiseEmissionsMap;
    private Map<Double, Double> maxNoiseEmissionValue;
    /**
     *
     * @param scaledMarginalUtilityOfTraveling Must be scaled, i.e. per second.  Usually negative.
     * @param scaledMarginalUtilityOfPerforming Must be scaled, i.e. per second.  Usually positive.
     * @param scaledMarginalUtilityOfDistance Must be scaled, i.e. per meter.  Usually negative.
     * @param noiseEmissionResults Please use MATSim noise contrib to generate the noise emission results.
     */
    public NoiseBasedTravelDisutility(double scaledMarginalUtilityOfTraveling, double scaledMarginalUtilityOfPerforming,
                                            double scaledMarginalUtilityOfDistance, String noiseEmissionResults, double timeBin, int simulationEndHour){
        // usually, the travel-utility should be negative (it's a disutility)
        // but for the cost, the cost should be positive.
        this.travelCostFactor = -scaledMarginalUtilityOfTraveling + scaledMarginalUtilityOfPerforming;

        //if ( wrnCnt < 1 ) {
            //wrnCnt++ ;
            if (this.travelCostFactor <= 0) {
                log.warn("The travel cost in " + this.getClass().getName() + " under normal circumstances should be > 0. " +
                        "Currently, it is " + this.travelCostFactor + "." +
                        "That is the sum of the costs for traveling and the opportunity costs." +
                        " Please adjust the parameters" +
                        "'traveling' and 'performing' in the module 'planCalcScore' in your config file to be" +
                        " lower or equal than 0 when added.");
                log.warn(Gbl.ONLYONCE) ;
            }
        //}

        this.marginalUtlOfDistance = scaledMarginalUtilityOfDistance;
        this.noiseEmissionResults = noiseEmissionResults;
        this.timeBin = timeBin;
        this.simulationEndHour = simulationEndHour;
        this.readNoiseEmissionResults();
    }

/*    public NoiseBasedTravelDisutility(PlanCalcScoreConfigGroup cnScoringGroup){
        this(cnScoringGroup.getModes().get(TransportMode.car).getMarginalUtilityOfTraveling() / 3600.0, cnScoringGroup.getPerforming_utils_hr() / 3600.0,
//				cnScoringGroup.getMarginalUtlOfDistanceCar());
                cnScoringGroup.getModes().get(TransportMode.car).getMonetaryDistanceRate() *cnScoringGroup.getMarginalUtilityOfMoney());
    }*/

    // read the necessary input data
    public void readNoiseEmissionResults() {

        for (int i=1;i<=simulationEndHour;i++) {
            String noiseEmissionResultsFile = noiseEmissionResults + "emission_" + this.timeBin*i + ".csv";
            List<String[]> nodes = CSVReaders.readCSV(noiseEmissionResultsFile);

            Map<Id<Link>, Double> newEntry = new HashMap<>();
            // read noise emission values of the links by the time of the day
            for (String[] line : Iterables.skip(nodes, 1)) { // skip CSV header
                if (newEntry.containsKey(Id.createLinkId(line[0]))){
                    throw new RuntimeException("The link " + Id.createLinkId(line[0]).toString() + "has multiple noise emission values!");
                } else {
                    newEntry.put(Id.createLinkId(line[0]), Double.parseDouble(line[5]));
                }
            }
            this.noiseEmissionsMap.put(timeBin*i, newEntry);

            // calculate the maximal value
            Map.Entry<Id<Link>, Double> maxEntry = null;
            for (Map.Entry<Id<Link>, Double> entry : newEntry.entrySet()) {
                if (maxEntry == null || entry.getValue().compareTo(maxEntry.getValue()) > 0) {
                    maxEntry = entry;
                }
            }
            if(maxEntry != null) {
                maxNoiseEmissionValue.put(timeBin * i, maxEntry.getValue());
            } else {
                throw new RuntimeException("For the timeBin " + timeBin*i + ": The maximal value is null!");
            }

        }
    }

    @Override
    public double getLinkTravelDisutility(Link link, double time, Person person, Vehicle vehicle) {
        double defaultCost = 0.;
/*        if (this.marginalUtlOfDistance == 0.0) {
            defaultCost = (link.getLength() / link.getFreespeed(time)) * this.travelCostFactor;
        } else {
            defaultCost = (link.getLength() / link.getFreespeed(time)) * this.travelCostFactor - this.marginalUtlOfDistance * link.getLength();
        }*/

        // add noise based cost
        double totalCost = defaultCost;
        double timeWindow = Math.floor(time/timeBin)*timeBin;
        double maximalNoiseEmissionValueForThisTimePoint = maxNoiseEmissionValue.get(timeWindow);
        totalCost = maximalNoiseEmissionValueForThisTimePoint - this.noiseEmissionsMap.get(timeWindow).get(link.getId())
        + totalCost;
        return totalCost;
    }

    // ToDo: This method is currently unused!
    @Override
    public double getLinkMinimumTravelDisutility(Link link) {
        return 0;
    }

    // The following codes are only for the travel times!
    @Override
    public double getLinkTravelTime(Link link, double time, Person person, Vehicle vehicle) {
        return link.getLength() / link.getFreespeed(time);
    }

    /**
     * If travelling freespeed the turning move travel time is not relevant
     */
    @Override
    public double getLinkToLinkTravelTime(Link fromLink, Link toLink, double time, Person person, Vehicle vehicle) {
        return this.getLinkTravelTime(fromLink, time, null, null);
    }
}
