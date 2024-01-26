package net.bhl.matsim.uam.optimization;

import net.bhl.matsim.uam.analysis.traveltimes.utils.TripItem;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static net.bhl.matsim.uam.optimization.VertiportOptimizerGreedyForwards.*;

public class testGreedy {
    public static final double flightSpeed= 350/3.6; // m/s
    public static final double UAM_PROCESS_TIME= 10*60; // s
    public static final double takeOffLandingTime= 60; // s
    private static String vertiportCandidateFile;
    public static final double CURRENT_TOTAL_GENERALIZED_COST= 7703039.047913355;
    private static String fileName;

    public static void main(String[] args) throws IOException {
        List<Integer> chosenVertiportID=List.of(49, 158, 84, 32, 142, 170, 155, 18, 21, 24, 48, 133, 119, 143, 56, 167, 169, 183, 157, 15, 112, 101, 182, 20, 194, 131, 120, 47, 33, 83, 96, 98, 185, 62, 14, 177, 11, 136, 60, 37, 25, 91, 72, 68, 58, 173, 132, 191, 35, 127, 63, 34, 90, 121, 176, 4, 125, 19, 171, 3, 36, 102, 134, 103, 80, 40, 28, 181, 123, 95, 43, 187, 31, 9);

        // Provide the file via program arguments
        if (args.length > 0) {
            fileName = args[0];
            vertiportCandidateFile = args[1];
        }
        // Get the object TripItem from the serialized file

        System.out.println("Loading the vertiport candidates...");
        VertiportReader vertiportReader = new VertiportReader();
        List<Vertiport> vertiportsCandidates = VertiportReader.getVertiports(vertiportCandidateFile);
        // Test
        // Only the first 50 vertiports are used for testing
        //      vertiportsCandidates = vertiportsCandidates.subList(0, 50);
        System.out.println("Finished loading the vertiport candidates.");
        System.out.println("Loading the trips...");
        List<TripItem> deserializedTripItems = deserializeTripItems(fileName);
        System.out.println("Finished loading the trips.");
        System.out.println("Calculating the generalized cost of each trip...");
        System.out.println(calculateSelectionScore(vertiportsCandidates,chosenVertiportID,deserializedTripItems));
    }
        public static double calculateSelectionScore( List<Vertiport> vertiportCandidate,List<Integer> chosenVertiportID, List< TripItem > deserializedTripItems) {
            // 实现适应度函数的具体逻辑
            double sumGeneralizedCost=0.0;
            double sumVertiportConstructionCost=0.0;
            double savedGeneralizedCost=0.0;
            double score=0.0;
            List<TripItem> uamAvailableTrips = new ArrayList<>();
            // calculate the sum of vertiport construction cost
            for (Integer vertiportID:chosenVertiportID) {
                sumVertiportConstructionCost=sumVertiportConstructionCost+vertiportCandidate.get(vertiportID).constructionCost;
            }

            for (TripItem tripItem : deserializedTripItems)
            {
                tripItem.isUAMAvailable = false;
                tripItem.uamUtility=-9999;
                List<Vertiport> originNeighbourVertiports = findAvailableNeighbourVertiports(chosenVertiportID, tripItem.originNeighborVertiportCandidates);
                List<Vertiport> destinationNeighbourVertiports = findAvailableNeighbourVertiports(chosenVertiportID, tripItem.destinationNeighborVertiportCandidates);
                if (originNeighbourVertiports.size() > 0 && destinationNeighbourVertiports.size() > 0) {
                    tripItem.isUAMAvailable = true;
                    tripItem.originNeighborVertiports = originNeighbourVertiports;
                    tripItem.destinationNeighborVertiports = destinationNeighbourVertiports;
                    uamAvailableTrips.add(tripItem);}
            }

            for (TripItem tripItem : uamAvailableTrips) {

                    for (Vertiport vertiport : tripItem.originNeighborVertiports) {
                        tripItem.originNeighborVertiportsTimeAndDistance.put(vertiport,tripItem.originNeighborVertiportCandidatesTimeAndDistance.get(vertiport));
                    }
                    for (Vertiport vertiport : tripItem.destinationNeighborVertiports) {
                        tripItem.destinationNeighborVertiportsTimeAndDistance.put(vertiport,tripItem.destinationNeighborVertiportCandidatesTimeAndDistance.get(vertiport));
                    }

                    // Initialize the Vertiport Allocation
                    // Find the vertiport pair for a trip with the lowest uam generalized cost, the access and ergress vertiport could not be the same
                    double lowestUAMGeneralizedCost = Double.MAX_VALUE;
                    Vertiport originVertiport = null;
                    Vertiport destinationVertiport = null;
                    for (Vertiport origin : tripItem.originNeighborVertiports) {
                        for (Vertiport destination :  tripItem.destinationNeighborVertiports) {
                            if (origin.ID != destination.ID) {
                                double accessTime = tripItem.originNeighborVertiportsTimeAndDistance.get(origin).get("travelTime");
                                double egressTime = tripItem.destinationNeighborVertiportsTimeAndDistance.get(destination).get("travelTime");
                                double accessDistance = tripItem.originNeighborVertiportsTimeAndDistance.get(origin).get("distance");
                                double egressDistance = tripItem.destinationNeighborVertiportsTimeAndDistance.get(destination).get("distance");
                                double accessCost =0;
                                double egressCost =0;
                                double flightDistance= calculateEuciDistance(origin.coord,destination.coord);
                                double flightTime=flightDistance/flightSpeed+takeOffLandingTime;
                                double flightCost=6.1+ calculateEuciDistance(origin.coord,destination.coord)/1000*0.6;
                                double uamTravelTime=accessTime+egressTime+flightTime+UAM_PROCESS_TIME;
                                if (tripItem.accessMode.equals("car") ){
                                    accessCost=accessDistance/1000*0.42;
                                }
                                if (tripItem.egressMode.equals("car") ){
                                    egressCost=egressDistance/1000*0.42;
                                }
                                double UAMCost=accessCost+egressCost+flightCost;
                                double UAMGeneralizedCost=UAMCost+uamTravelTime*tripItem.VOT;
                                if (UAMGeneralizedCost < lowestUAMGeneralizedCost) {
                                    tripItem.uamTravelTime=uamTravelTime;
                                    tripItem.UAMCost=UAMCost;
                                    tripItem.UAMUtilityVar=-2.48*UAMCost/100-4.28*flightTime/6000-6.79*(uamTravelTime-flightTime)/6000;
                                    tripItem.uamUtility=tripItem.UAMUtilityFix+tripItem.UAMUtilityVar;
                                    tripItem.UAMGeneralizedCost=UAMGeneralizedCost;
                                    tripItem.accessVertiport = origin;
                                    tripItem.egressVertiport = destination;
                                    lowestUAMGeneralizedCost = UAMGeneralizedCost;
                                }
                            }
                        }
                    }
                    // determine the probability of mode choice of each trip
                    tripItem.uamProbability=calculateModeProbability(tripItem.uamUtility,tripItem.carUtility,tripItem.ptUtility).get(0);
                    tripItem.carProbability=calculateModeProbability(tripItem.uamUtility,tripItem.carUtility,tripItem.ptUtility).get(1);
                    tripItem.ptProbability=calculateModeProbability(tripItem.uamUtility,tripItem.carUtility,tripItem.ptUtility).get(2);
                    double generalizedCostOneTripBefore=tripItem.carGeneralizedCost*calculateModeProbability(-9999,tripItem.carUtility,tripItem.ptUtility).get(1)+tripItem.ptGeneralizedCost*calculateModeProbability(-9999,tripItem.carUtility,tripItem.ptUtility).get(2);
                    double generalizedCostOneTripAfter=tripItem.UAMGeneralizedCost*tripItem.uamProbability+tripItem.carGeneralizedCost*tripItem.carProbability+tripItem.ptGeneralizedCost*tripItem.ptProbability;
                    double savedGeneralizedCostOneTrip=generalizedCostOneTripBefore-generalizedCostOneTripAfter;
                    if (savedGeneralizedCostOneTrip<0){
                        savedGeneralizedCostOneTrip=0;
                    }
                    savedGeneralizedCost=savedGeneralizedCost+savedGeneralizedCostOneTrip;


                }



            /*
            int carCount=0;
            int ptCount=0;
            int uamCount=0;
            // generate 100 scenarios of each trip
            int j=0;
            for (int modeChoiceIterator=0;modeChoiceIterator<100;modeChoiceIterator++){
                j=modeChoiceIterator;
                ModeDecider modeDecider=new ModeDecider(tripItem.carProbability,tripItem.ptProbability,tripItem.uamProbability);
                String mode=modeDecider.decideMode();
                if (mode.equals("car")){
                    generalizedCostOneTrip=generalizedCostOneTrip+tripItem.carGeneralizedCost;
                    carCount++;
                }
                if (mode.equals("pt")){
                    generalizedCostOneTrip=generalizedCostOneTrip+tripItem.ptGeneralizedCost;
                    ptCount++;
                }
                if (mode.equals("uam")){
                    generalizedCostOneTrip=generalizedCostOneTrip+tripItem.UAMGeneralizedCost;
                    uamCount++;
                }
                if (tripItem.carProbability>0.9999 || tripItem.ptProbability>0.9999 || tripItem.uamProbability>0.9999){
                    carCount=carCount*100;
                    ptCount=ptCount*100;
                    uamCount=uamCount*100;
                    break;
                }
            }
            // calculate the average generalized cost of each trip
            generalizedCostOneTrip=generalizedCostOneTrip/(j+1);

             */




            score=savedGeneralizedCost*2;
            return score;
        }
    }

