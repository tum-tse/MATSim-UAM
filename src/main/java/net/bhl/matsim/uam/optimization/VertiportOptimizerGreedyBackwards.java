package net.bhl.matsim.uam.optimization;

import net.bhl.matsim.uam.analysis.traveltimes.utils.TripItem;
import org.apache.log4j.Logger;
import org.matsim.api.core.v01.Coord;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.*;

public class VertiportOptimizerGreedyBackwards {
    public static final double flightSpeed= 350/3.6; // m/s
    public static final double UAM_PROCESS_TIME= 10*60; // s
    public static final double takeOffLandingTime= 60; // s
    private static String vertiportCandidateFile;
    public static final double CURRENT_TOTAL_GENERALIZED_COST= 7703039.047913355;
    private static String fileName;
    public static final Logger log = Logger.getLogger(VertiportOptimizerGreedyBackwards.class);

    public static double calculateEuciDistance(Coord coord1, Coord coord2) {
        double euciDistance = Math.sqrt(Math.pow(coord1.getX() - coord2.getX(), 2) + Math.pow(coord1.getY() - coord2.getY(), 2));
        return euciDistance;
    }
    public static void main(String[] args) throws IOException {
        // Provide the file via program arguments
        if (args.length > 0) {
            fileName = args[0];
            vertiportCandidateFile = args[1];
        }
        // Get the object TripItem from the serialized file

        log.info("Loading the vertiport candidates...");
        VertiportReader vertiportReader = new VertiportReader();
        List<Vertiport> vertiportsCandidates = VertiportReader.getVertiports(vertiportCandidateFile);
        // Test
        // Only the first 50 vertiports are used for testing
        //      vertiportsCandidates = vertiportsCandidates.subList(0, 50);
        log.info("Finished loading the vertiport candidates.");
        log.info("Loading the trips...");
        List<TripItem> deserializedTripItems = deserializeTripItems(fileName);
        log.info("Finished loading the trips.");
// Test
        // Only the first 5000 trips are used for testing
//        deserializedTripItems = deserializedTripItems.subList(0, 5000);
        HashMap<Integer, Double> vertiportsScore = new HashMap<>();
        // Initialize the vertiportsScore with each vertiport candidate has a score of 0
        for (Vertiport vertiport : vertiportsCandidates) {
            vertiportsScore.put(vertiport.ID, 0.0);
        }
        // Assume that all vertiport candidates are selected
        List<Integer> chosenVertiportID = new ArrayList<>();
        for (Vertiport vertiport : vertiportsCandidates) {
            chosenVertiportID.add(vertiport.ID);
        }
        log.info("Start Iterate through the trips");
        int count = 0;
        List<TripItem> uamAvailableTrips = new ArrayList<>();
        for (TripItem tripItem : deserializedTripItems) {
            tripItem.isUAMAvailable = false;
            tripItem.uamUtility = -9999;
            List<Vertiport> originNeighbourVertiports = findAvailableNeighbourVertiports(chosenVertiportID, tripItem.originNeighborVertiportCandidates);
            List<Vertiport> destinationNeighbourVertiports = findAvailableNeighbourVertiports(chosenVertiportID, tripItem.destinationNeighborVertiportCandidates);
            if (originNeighbourVertiports.size() > 0 && destinationNeighbourVertiports.size() > 0) {
                if (originNeighbourVertiports.size() > 1 || destinationNeighbourVertiports.size() > 1) {
                    tripItem.isUAMAvailable = true;
                    tripItem.originNeighborVertiports = originNeighbourVertiports;
                    tripItem.destinationNeighborVertiports = destinationNeighbourVertiports;
                    uamAvailableTrips.add(tripItem);
                } else if (! (originNeighbourVertiports.get(0).equals(destinationNeighbourVertiports.get(0)))) {
                    tripItem.isUAMAvailable = true;
                    tripItem.originNeighborVertiports = originNeighbourVertiports;
                    tripItem.destinationNeighborVertiports = destinationNeighbourVertiports;
                    uamAvailableTrips.add(tripItem);
                }

            }

        }
        log.info("Finished iterating through the trips. Number of trips available for UAM: " + uamAvailableTrips.size());
        log.info("Start calculating the score of each vertiport.");
        for (TripItem tripItem : uamAvailableTrips) {

            for (Vertiport vertiport : tripItem.originNeighborVertiports) {
                tripItem.originNeighborVertiportsTimeAndDistance.put(vertiport, tripItem.originNeighborVertiportCandidatesTimeAndDistance.get(vertiport));
            }
            for (Vertiport vertiport : tripItem.destinationNeighborVertiports) {
                tripItem.destinationNeighborVertiportsTimeAndDistance.put(vertiport, tripItem.destinationNeighborVertiportCandidatesTimeAndDistance.get(vertiport));
            }

            // Initialize the Vertiport Allocation
            // Find the vertiport pair for a trip with the lowest uam generalized cost, the access and ergress vertiport could not be the same
            double lowestUAMGeneralizedCost = Double.MAX_VALUE;
            for (Vertiport origin : tripItem.originNeighborVertiports) {
                for (Vertiport destination : tripItem.destinationNeighborVertiports) {
                    if (origin.ID != destination.ID) {
                        double accessTime = tripItem.originNeighborVertiportsTimeAndDistance.get(origin).get("travelTime");
                        double egressTime = tripItem.destinationNeighborVertiportsTimeAndDistance.get(destination).get("travelTime");
                        double accessDistance = tripItem.originNeighborVertiportsTimeAndDistance.get(origin).get("distance");
                        double egressDistance = tripItem.destinationNeighborVertiportsTimeAndDistance.get(destination).get("distance");
                        double accessCost = 0;
                        double egressCost = 0;
                        double flightDistance = calculateEuciDistance(origin.coord, destination.coord);
                        double flightTime = flightDistance / flightSpeed + takeOffLandingTime;
                        double flightCost = 6.1 + calculateEuciDistance(origin.coord, destination.coord) / 1000 * 0.6;
                        double uamTravelTime = accessTime + egressTime + flightTime + UAM_PROCESS_TIME;
                        if (tripItem.accessMode.equals("car")) {
                            accessCost = accessDistance / 1000 * 0.42;
                        }
                        if (tripItem.egressMode.equals("car")) {
                            egressCost = egressDistance / 1000 * 0.42;
                        }
                        double UAMCost = accessCost + egressCost + flightCost;
                        double UAMGeneralizedCost = UAMCost + uamTravelTime * tripItem.VOT;
                        if (UAMGeneralizedCost < lowestUAMGeneralizedCost) {
                            tripItem.uamTravelTime = uamTravelTime;
                            tripItem.UAMCost = UAMCost;
                            tripItem.UAMUtilityVar = -2.48 * UAMCost / 100 - 4.28 * flightTime / 6000 - 6.79 * (uamTravelTime - flightTime) / 6000;
                            tripItem.uamUtility = tripItem.UAMUtilityFix + tripItem.UAMUtilityVar;
                            tripItem.UAMGeneralizedCost = UAMGeneralizedCost;
                            tripItem.accessVertiport = origin;
                            tripItem.egressVertiport = destination;
                            lowestUAMGeneralizedCost = UAMGeneralizedCost;
                        }
                    }
                }
            }
            // determine the probability of mode choice of each trip
            tripItem.uamProbability = calculateModeProbability(tripItem.uamUtility, tripItem.carUtility, tripItem.ptUtility).get(0);
            tripItem.carProbability = calculateModeProbability(tripItem.uamUtility, tripItem.carUtility, tripItem.ptUtility).get(1);
            tripItem.ptProbability = calculateModeProbability(tripItem.uamUtility, tripItem.carUtility, tripItem.ptUtility).get(2);
            double generalizedCostOneTripBefore = tripItem.carGeneralizedCost * calculateModeProbability(-9999, tripItem.carUtility, tripItem.ptUtility).get(1) + tripItem.ptGeneralizedCost * calculateModeProbability(-9999, tripItem.carUtility, tripItem.ptUtility).get(2);
            double generalizedCostOneTripAfter = tripItem.UAMGeneralizedCost * tripItem.uamProbability + tripItem.carGeneralizedCost * tripItem.carProbability + tripItem.ptGeneralizedCost * tripItem.ptProbability;
            double savedGeneralizedCostOneTrip = generalizedCostOneTripBefore - generalizedCostOneTripAfter;
            if (savedGeneralizedCostOneTrip < 0) {
                savedGeneralizedCostOneTrip = 0;
            }
            if (tripItem.tripPurpose.startsWith("H")){
                savedGeneralizedCostOneTrip=savedGeneralizedCostOneTrip*2;
            }
            double newAccessVertiportScore=vertiportsScore.get(tripItem.accessVertiport.ID)+savedGeneralizedCostOneTrip/2;
            double newEgressVertiportScore=vertiportsScore.get(tripItem.egressVertiport.ID)+savedGeneralizedCostOneTrip/2;
            vertiportsScore.put(tripItem.accessVertiport.ID,newAccessVertiportScore);
            vertiportsScore.put(tripItem.egressVertiport.ID,newEgressVertiportScore);

            count++;
            if (count % 10000 == 0) {
                log.info("Finished calculating the score of " + count + "/"+ uamAvailableTrips.size()+ " trips.");
            }
         //   if (count >34200 ){
           //     System.out.println("Finished calculating the score of " + count + "/"+ uamAvailableTrips.size()+ " trips.");
            //}
        }
        log.info("Finished calculating the score of each vertiport.");

        log.info("Start sorting the vertiports based on the score.");
        List<Map.Entry<Integer, Double>> list = new ArrayList<Map.Entry<Integer, Double>>(vertiportsScore.entrySet());
        Collections.sort(list, new Comparator<Map.Entry<Integer, Double>>() {
            // descending order
            public int compare(Map.Entry<Integer, Double> o1, Map.Entry<Integer, Double> o2) {
                return o2.getValue().compareTo(o1.getValue());
            }

        });
        log.info("Finished sorting the vertiports based on the score.");
        log.info("Start selecting the vertiports.");
        List<Integer> selectedVertiportID = new ArrayList<>();
        // select the first 74 vertiports
        for (int i = 0; i < 74; i++) {
            selectedVertiportID.add(list.get(i).getKey());
        }

        log.info("Finished selecting the vertiports.");
        log.info("The selected vertiports are: " + selectedVertiportID);
        log.info("The score of selected vertiports are: "+ calculateSelectionScore(vertiportsCandidates,selectedVertiportID,deserializedTripItems));

    }





    public static List<TripItem> deserializeTripItems(String fileName) {
        List<TripItem> tripItems = new ArrayList<>();

        try  {
            FileInputStream fileIn = new FileInputStream(fileName);
            ObjectInputStream objectIn = new ObjectInputStream(fileIn)  ;
            tripItems = (List<TripItem>) objectIn.readObject();


        } catch (Exception e) {
            e.printStackTrace();
        }

        return tripItems;
    }
    public static List<Vertiport> findAvailableNeighbourVertiports(List<Integer> intList, List<Vertiport> vertiportList) {
        List<Vertiport> duplicates = new ArrayList<>();

        for (Vertiport vertiport : vertiportList) {
            if (intList.contains(vertiport.ID)) {
                duplicates.add(vertiport);
            } else {
                continue;
            }
        }

        return duplicates;
    }
    public static List<Integer> getRandomNumbers(int min, int max, int n) {
        List<Integer> list = new ArrayList<Integer>();
        for (int i = 0; i < n; i++) {
            int num = (int) (Math.random() * (max - min)) + min;
            if (!list.contains(num)) {
                list.add(num);
            } else {
                i--;
            }
        }
        return list;
    }

    public static List<Double> calculateModeProbability(double UAM_Utlility, double carUtility, double ptUtility){


        double sumUtilityExponential = Math.exp(carUtility) + Math.exp(ptUtility) + Math.exp(UAM_Utlility);
        double UAMProbability=Math.exp(UAM_Utlility)/ sumUtilityExponential;
        double carProbability=Math.exp(carUtility)/ sumUtilityExponential;
        double ptProbability=Math.exp(ptUtility)/ sumUtilityExponential;
        List<Double> modeProbability=new ArrayList<>();
        modeProbability.add(UAMProbability);
        modeProbability.add(carProbability);
        modeProbability.add(ptProbability);
        return modeProbability;
    }

    public static double calculateSelectionScore( List<Vertiport> vertiportCandidate,List<Integer> chosenVertiportID, List<TripItem> deserializedTripItems) {
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
            if (tripItem.tripPurpose.startsWith("H")){
                savedGeneralizedCostOneTrip=savedGeneralizedCostOneTrip*2;
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




        score=savedGeneralizedCost;
        return score;
    }






}