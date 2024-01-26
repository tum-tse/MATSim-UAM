package net.bhl.matsim.uam.optimization;

import net.bhl.matsim.uam.analysis.traveltimes.utils.TripItem;
import org.apache.log4j.Logger;
import org.matsim.api.core.v01.Coord;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;



public class SimulatedAnnealing {
    // This class provides the simulated annealing algorithm to solve the vertiport siting problem.

    public static final Logger log = Logger.getLogger(SimulatedAnnealing.class);
    public static final int NUM_OF_CANDIDATE_VERTIPORTS=200;
    public static final int NUM_OF_SELECTED_VERTIPORTS=74;
    public static final double INITIAL_TEMPERATURE=5000;
    public static final double FINAL_TEMPERATURE=0.1;
    public static final double ANNEALING_RATE=0.999;
    public static String serializedTripItemFile;
    public static String vertiportCandidateFile;
    public static final double flightSpeed= 350/3.6; // m/s
    public static final double UAM_PROCESS_TIME= 10*60; // s
    public static final double takeOffLandingTime= 60; // s


    public static void main(String[] args) throws IOException {
        // Provide the file via program arguments
        if (args.length > 0) {
            serializedTripItemFile = args[0];
            vertiportCandidateFile = args[1];
        }
        // Get the object TripItem from the serialized file

        log.info("Loading the vertiport candidates...");
        VertiportReader vertiportReader = new VertiportReader();
        List<Vertiport> vertiportsCandidates = VertiportReader.getVertiports(vertiportCandidateFile);

        log.info("Finished loading the vertiport candidates.");
        log.info("Loading the trips...");
        List<TripItem> deserializedTripItems = deserializeTripItems(serializedTripItemFile);
        log.info("Finished loading the trips.");

        List<Integer> currentSolution = generateRandomSolution();
        double currentEnergey = calculateFitness(currentSolution,vertiportsCandidates,deserializedTripItems);
        List<Integer> bestSolution = new ArrayList<>(currentSolution);
        double bestEnergy= currentEnergey;
        double currentTemperature = INITIAL_TEMPERATURE;
        int notChangeCount=0;
        for (int iteration=0;iteration<10000;iteration++ ){
            List<Integer> newSolution= generateNewSolution(currentSolution);
            double newEnergy= calculateFitness(newSolution,vertiportsCandidates,deserializedTripItems);
            double deltaEnergy= newEnergy-currentEnergey;
            Random random=new Random();
            if(deltaEnergy>0 || Math.exp(deltaEnergy/currentTemperature)>random.nextDouble()){
                currentSolution=new ArrayList<>(newSolution);
                currentEnergey=newEnergy;
            }
            // update the best solution
            if (currentEnergey>bestEnergy){
                bestSolution=new ArrayList<>(currentSolution);
                bestEnergy=currentEnergey;
                notChangeCount=0;
            }
            else {
                notChangeCount++;
            }


            // update the temperature
            currentTemperature=currentTemperature*ANNEALING_RATE;
            log.info("Iteration: "+iteration+" Current Temperature: "+currentTemperature+" Current Energy: "+currentEnergey+" Best Energy: "+bestEnergy);

            // if the best solution is not updated for more than 2000 iterations, break
            if(notChangeCount>2000){
                break;
            }
        }
        log.info("Best Solution: "+ bestSolution + " Best Energy: "+ bestEnergy);


    }

    public static List<Integer> generateRandomSolution() {
        // generate a Integer List with a length of 200, 74 of them are 1, the rest are 0
        List<Integer> solution = new ArrayList<>();
        Random random = new Random();
        while (solution.size() < NUM_OF_SELECTED_VERTIPORTS) {
            int randomNum = random.nextInt(NUM_OF_CANDIDATE_VERTIPORTS);
            if(!solution.contains(randomNum)) {
            solution.add(randomNum);}
        }
        return solution;
    }

    public static double calculateFitness(List<Integer> chosenVertiportID,List<Vertiport> vertiportCandidates, List<TripItem> deserializedTripItems) {
        // 实现适应度函数的具体逻辑

        double sumVertiportConstructionCost=0.0;
        double savedGeneralizedCost=0.0;

        List<TripItem> uamAvailableTrips = new ArrayList<>();
        // calculate the sum of vertiport construction cost
        for (Integer vertiportID:chosenVertiportID){
            sumVertiportConstructionCost=sumVertiportConstructionCost+vertiportCandidates.get(vertiportID).constructionCost;
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

       return savedGeneralizedCost;


}
public static List<Integer> generateNewSolution (List<Integer> currentSolution){
        List<Integer> newSolution=new ArrayList<>(currentSolution);
        List<Integer> notChosenVertiport = new ArrayList<>();
        for (int i=0;i<200;i++){
            if (!newSolution.contains(i)){
                notChosenVertiport.add(i);
            }
        }
        // randomly substitute one element in newSolution with one element in notChosenVertiport
    Random random=new Random();
        int index1= random.nextInt(newSolution.size());
        int index2= random.nextInt(notChosenVertiport.size());

        newSolution.set(index1,notChosenVertiport.get(index2));

    return newSolution;
}
    public static double calculateEuciDistance(Coord coord1, Coord coord2) {
        double euciDistance = Math.sqrt(Math.pow(coord1.getX() - coord2.getX(), 2) + Math.pow(coord1.getY() - coord2.getY(), 2));
        return euciDistance;
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


}
