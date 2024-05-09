package net.bhl.matsim.uam.optimization;

import net.bhl.matsim.uam.optimization.utils.TripItemForOptimization;
import org.apache.log4j.Logger;
import org.matsim.api.core.v01.Coord;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import org.matsim.utils.MemoryObserver;


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
    public static final double CAR_EMISSION_FACTOR=0.42; // kg/km
    public static final double PT_EMISSION_FACTOR=0.1; // kg/km
    public static final double UAM_EMISSION_FACTOR=0.1; // kg/km
    public static final double CARBON_EQUIVALENCE_FACTOR=2.48; // Euro/kgCO2
    public static final boolean CONSIDER_CARBON= false; // Euro/kgCO2
    public static void main(String[] args) throws IOException {
        MemoryObserver.start(600);
        // Provide the file via program arguments
        if (args.length > 0) {
            serializedTripItemFile = args[0];
            vertiportCandidateFile = args[1];
        }
        // Get the object TripItemForOptimization from the serialized file

        log.info("Loading the vertiport candidates...");
        VertiportReader vertiportReader = new VertiportReader();
        List<Vertiport> vertiportsCandidates = vertiportReader.getVertiportsWithNeighbors(vertiportCandidateFile);
        log.info("Finished loading the vertiport candidates.");
        log.info("Loading the trips...");
        List<TripItemForOptimization> deserializedTripItemForOptimizations = deserializeTripItems(serializedTripItemFile);
        log.info("Finished loading the trips.");
        Random random = new Random(42);
        List<Integer> currentSolutionID = generateRandomSolution(random,vertiportsCandidates,NUM_OF_SELECTED_VERTIPORTS);
        double currentEnergey = calculateFitness(currentSolutionID, vertiportsCandidates,deserializedTripItemForOptimizations);
        List<Integer> bestSolutionID = new ArrayList<>(currentSolutionID);
        double bestEnergy= currentEnergey;
        double currentTemperature = INITIAL_TEMPERATURE;
        int notChangeCount=0;
        int saturatedVertiportCount=0;
        double maxSaturationRate=0;
        for (Integer vertiportID:currentSolutionID){
            if (vertiportsCandidates.get(vertiportID).maxSaturationRate>1){
                saturatedVertiportCount++;
            }
            if (vertiportsCandidates.get(vertiportID).maxSaturationRate>maxSaturationRate){
                maxSaturationRate=vertiportsCandidates.get(vertiportID).maxSaturationRate;
            }
        }
        for (int iteration=0;iteration<10000;iteration++ ){
            List<Integer> newSolutionID= generateNewSolution(random,currentSolutionID,vertiportsCandidates,CONSIDER_CARBON);
            // set the saturation rate of all vertiports to 0
            for (Vertiport vertiport:vertiportsCandidates){
                for (int i=0;i<24;i++){
                    vertiport.saturationRates.put(i,0.0);
                }
            }
            double newEnergy= calculateFitness(newSolutionID,vertiportsCandidates, deserializedTripItemForOptimizations);
            double deltaEnergy= newEnergy-currentEnergey;


            if(deltaEnergy>0 || Math.exp(deltaEnergy/currentTemperature)>random.nextDouble()){
                currentSolutionID=new ArrayList<>(newSolutionID);
                currentEnergey=newEnergy;
                for (Integer vertiportID:currentSolutionID){
                    if (vertiportsCandidates.get(vertiportID).maxSaturationRate>1){
                        saturatedVertiportCount++;
                    }
                    if (vertiportsCandidates.get(vertiportID).maxSaturationRate>maxSaturationRate){
                        maxSaturationRate=vertiportsCandidates.get(vertiportID).maxSaturationRate;
                    }
                }

            }
            // update the best solution
            if (currentEnergey>bestEnergy){
                bestSolutionID=new ArrayList<>(currentSolutionID);
                bestEnergy=currentEnergey;
                notChangeCount=0;
            }
            else {
                notChangeCount++;
            }


            // update the temperature
            currentTemperature=currentTemperature*ANNEALING_RATE;
            log.info("Iteration: "+iteration+" Current Temperature: "+currentTemperature+" Current Energy: "+currentEnergey+" Best Energy: "+bestEnergy+ " Saturated Vertiport Count: "+saturatedVertiportCount+" Max Saturation Rate: "+maxSaturationRate);

            // if the best solution is not updated for more than 2000 iterations, break
            if(notChangeCount>2000){
                break;
            }
            saturatedVertiportCount=0;
            maxSaturationRate=0;
        }

        log.info("Best Solution: "+ bestSolutionID + " Best Energy: "+ bestEnergy);


    }


    public static List<Integer> generateRandomSolution(Random random,List<Vertiport> VertiportCandidates, int numOfSelectedVertiports) {

      List<Integer> selectedVertiportsID = new ArrayList<>();
      while (selectedVertiportsID.size() < numOfSelectedVertiports) {
          int vertiportID = VertiportCandidates.get(random.nextInt(VertiportCandidates.size())).ID;
          if (!selectedVertiportsID.contains(vertiportID)) {
              selectedVertiportsID.add(vertiportID);
          }
      }
      return selectedVertiportsID;
    }

    public static double calculateFitness(List<Integer> chosenVertiportID, List<Vertiport> vertiportsCandidates,List<TripItemForOptimization> deserializedTripItemForOptimizations) {
        // 实现适应度函数的具体逻辑

        double sumVertiportConstructionCost=0.0;
        double savedGeneralizedCost=0.0;
// decode the chosenVertiportID to the vertiport object
        List<Vertiport> chosenVertiport = new ArrayList<>();
        for (Integer vertiportID:chosenVertiportID){
            chosenVertiport.add(vertiportsCandidates.get(vertiportID));
        }
        List<TripItemForOptimization> uamAvailableTrips = new ArrayList<>();
//        // calculate the sum of vertiport construction cost
//        for (Integer vertiportID:chosenVertiportID){
//            sumVertiportConstructionCost=sumVertiportConstructionCost+vertiportCandidates.get(vertiportID).constructionCost;
//        }

        for (TripItemForOptimization tripItemForOptimization : deserializedTripItemForOptimizations)
        {
            tripItemForOptimization.isUAMAvailable = false;
            tripItemForOptimization.uamUtility=-9999;
            List<Vertiport> originNeighbourVertiports = findAvailableNeighbourVertiports(chosenVertiportID, tripItemForOptimization.originNeighborVertiportCandidates);
            List<Vertiport> destinationNeighbourVertiports = findAvailableNeighbourVertiports(chosenVertiportID, tripItemForOptimization.destinationNeighborVertiportCandidates);
            if (originNeighbourVertiports.size() > 0 && destinationNeighbourVertiports.size() > 0) {
                if(originNeighbourVertiports.size()>1 || destinationNeighbourVertiports.size()>1){
                    tripItemForOptimization.isUAMAvailable = true;
                    tripItemForOptimization.originNeighborVertiports = new ArrayList<>();
                    tripItemForOptimization.destinationNeighborVertiports = new ArrayList<>();
                    tripItemForOptimization.originNeighborVertiports.addAll(originNeighbourVertiports);
                    tripItemForOptimization.destinationNeighborVertiports.addAll(destinationNeighbourVertiports);
                    uamAvailableTrips.add(tripItemForOptimization);
                }
                else if (originNeighbourVertiports.get(0).ID != destinationNeighbourVertiports.get(0).ID) {
                    tripItemForOptimization.isUAMAvailable = true;
                    tripItemForOptimization.originNeighborVertiports = new ArrayList<>();
                    tripItemForOptimization.destinationNeighborVertiports = new ArrayList<>();
                    tripItemForOptimization.originNeighborVertiports.addAll(originNeighbourVertiports);
                    tripItemForOptimization.destinationNeighborVertiports.addAll(destinationNeighbourVertiports);
                    uamAvailableTrips.add(tripItemForOptimization);
                }
 }
        }

        for (TripItemForOptimization tripItemForOptimization : uamAvailableTrips) {
            if (tripItemForOptimization.tripID==1146924){
                System.out.println("tripID: "+tripItemForOptimization.tripID);
            }
            for (Vertiport vertiport : tripItemForOptimization.originNeighborVertiports) {
                tripItemForOptimization.originNeighborVertiportsTimeAndDistance.put(vertiport, tripItemForOptimization.originNeighborVertiportCandidatesTimeAndDistance.get(vertiport));
            }
            for (Vertiport vertiport : tripItemForOptimization.destinationNeighborVertiports) {
                tripItemForOptimization.destinationNeighborVertiportsTimeAndDistance.put(vertiport, tripItemForOptimization.destinationNeighborVertiportCandidatesTimeAndDistance.get(vertiport));
            }
            double objectiveFunctionBefore;
            double objectiveFunctionAfter;
            double depatureTime= tripItemForOptimization.departureTime;
            // Initialize the Vertiport Allocation
            // Find the vertiport pair for a trip with the lowest uam generalized cost, the access and ergress vertiport could not be the same
            double lowestUAMGeneralizedCost = Double.MAX_VALUE;
            Vertiport originVertiport = null;
            Vertiport destinationVertiport = null;
            for (Vertiport origin : tripItemForOptimization.originNeighborVertiports) {
                for (Vertiport destination :  tripItemForOptimization.destinationNeighborVertiports) {
                    if (origin.ID != destination.ID) {
                        double accessTime = tripItemForOptimization.originNeighborVertiportsTimeAndDistance.get(origin).get("travelTime");
                        double egressTime = tripItemForOptimization.destinationNeighborVertiportsTimeAndDistance.get(destination).get("travelTime");
                        double accessDistance = tripItemForOptimization.originNeighborVertiportsTimeAndDistance.get(origin).get("distance");
                        double egressDistance = tripItemForOptimization.destinationNeighborVertiportsTimeAndDistance.get(destination).get("distance");
                        double accessCost =0;
                        double egressCost =0;
                        double accessEmisson =0;
                        double egressEmisson =0;
                        double flightDistance= calculateEuciDistance(origin.coord,destination.coord);
                        double flightTime=flightDistance/flightSpeed+takeOffLandingTime;
                        double flightCost=6.1+ flightDistance/1000*0.6;
                        double flightEmission=flightDistance/1000*UAM_EMISSION_FACTOR;
                        double uamTravelTime=accessTime+egressTime+flightTime+UAM_PROCESS_TIME;

                        if (tripItemForOptimization.accessMode.equals("car") ){
                            accessCost=accessDistance/1000*0.42;
                            accessEmisson=accessDistance/1000*CAR_EMISSION_FACTOR;
                        }
                        if (tripItemForOptimization.egressMode.equals("car") ){
                            egressCost=egressDistance/1000*0.42;
                            egressEmisson=egressDistance/1000*CAR_EMISSION_FACTOR;
                        }
                        double UAMCost=accessCost+egressCost+flightCost;
                        double UAMGeneralizedCost=UAMCost+uamTravelTime* tripItemForOptimization.VOT;
                        if (CONSIDER_CARBON){
                            UAMGeneralizedCost=UAMGeneralizedCost+CARBON_EQUIVALENCE_FACTOR*(accessEmisson+egressEmisson+flightEmission);
                        }
                        if (UAMGeneralizedCost < lowestUAMGeneralizedCost) {
                            tripItemForOptimization.uamTravelTime=uamTravelTime;
                            tripItemForOptimization.UAMCost=UAMCost;
                            tripItemForOptimization.UAMUtilityVar=-2.48*UAMCost/100-4.28*flightTime/6000-6.79*(uamTravelTime-flightTime)/6000;
                            tripItemForOptimization.uamUtility= tripItemForOptimization.UAMUtilityFix+ tripItemForOptimization.UAMUtilityVar;
                            tripItemForOptimization.UAMGeneralizedCost=UAMGeneralizedCost;
                            tripItemForOptimization.accessVertiport = origin;
                            tripItemForOptimization.egressVertiport = destination;
                            tripItemForOptimization.accessTime=accessTime;
                            tripItemForOptimization.egressTime=egressTime;
                            tripItemForOptimization.flightTime=flightTime;
                            lowestUAMGeneralizedCost = UAMGeneralizedCost;
                        }
                    }
                }
            }
            // determine the probability of mode choice of each trip
            tripItemForOptimization.uamProbability=calculateModeProbability(tripItemForOptimization.uamUtility, tripItemForOptimization.carUtility, tripItemForOptimization.ptUtility).get(0);
            tripItemForOptimization.carProbability=calculateModeProbability(tripItemForOptimization.uamUtility, tripItemForOptimization.carUtility, tripItemForOptimization.ptUtility).get(1);
            tripItemForOptimization.ptProbability=calculateModeProbability(tripItemForOptimization.uamUtility, tripItemForOptimization.carUtility, tripItemForOptimization.ptUtility).get(2);
            Double [] modeProbability={tripItemForOptimization.carProbability, tripItemForOptimization.ptProbability, tripItemForOptimization.uamProbability};
            // update the saturation rate of the access and egress vertiport

            int arriveVertiportHour=(int)Math.floor((tripItemForOptimization.departureTime+ tripItemForOptimization.accessTime)/3600);
            if (arriveVertiportHour>=24){
                arriveVertiportHour=arriveVertiportHour-24;
            }
            int leaveVertiportHour=(int)Math.floor((tripItemForOptimization.departureTime+ tripItemForOptimization.accessTime+UAM_PROCESS_TIME+ tripItemForOptimization.flightTime)/3600);
            if (leaveVertiportHour>=24){
                leaveVertiportHour=leaveVertiportHour-24;
            }
            vertiportsCandidates.get(tripItemForOptimization.accessVertiport.ID).saturationRates.put(arriveVertiportHour,  vertiportsCandidates.get(tripItemForOptimization.accessVertiport.ID).saturationRates.get(arriveVertiportHour)+ tripItemForOptimization.uamProbability/ tripItemForOptimization.accessVertiport.capacity);
            vertiportsCandidates.get(tripItemForOptimization.egressVertiport.ID).saturationRates.put(leaveVertiportHour,  vertiportsCandidates.get(tripItemForOptimization.egressVertiport.ID).saturationRates.get(leaveVertiportHour)+ tripItemForOptimization.uamProbability/ tripItemForOptimization.egressVertiport.capacity);

            ModeDecider modeDecider=new ModeDecider(modeProbability);
            Double [] modeSamples=modeDecider.sample(100);
            if (!CONSIDER_CARBON)
            {objectiveFunctionBefore= tripItemForOptimization.currentGeneralizedCost;
            objectiveFunctionAfter= modeSamples[0]* tripItemForOptimization.carGeneralizedCost+modeSamples[1]* tripItemForOptimization.ptGeneralizedCost+modeSamples[2]* tripItemForOptimization.UAMGeneralizedCost;
            }
            else {
                objectiveFunctionBefore= tripItemForOptimization.currentGeneralizedCost+tripItemForOptimization.currentEmission*CARBON_EQUIVALENCE_FACTOR;;
                objectiveFunctionAfter= modeSamples[0]*(tripItemForOptimization.carGeneralizedCost+ tripItemForOptimization.carEmission*CARBON_EQUIVALENCE_FACTOR)+modeSamples[1]*(tripItemForOptimization.ptGeneralizedCost+ tripItemForOptimization.ptEmission*CARBON_EQUIVALENCE_FACTOR)+modeSamples[2]* tripItemForOptimization.UAMGeneralizedCost;
            }


            double savedGeneralizedCostOneTrip=objectiveFunctionBefore-objectiveFunctionAfter;

            if (savedGeneralizedCostOneTrip<0){
                savedGeneralizedCostOneTrip=0;
            }
            if (tripItemForOptimization.tripPurpose.startsWith("H")){
                savedGeneralizedCostOneTrip=savedGeneralizedCostOneTrip*2;
            }
            savedGeneralizedCost=savedGeneralizedCost+savedGeneralizedCostOneTrip;




        }


       return savedGeneralizedCost;


}
public static List<Integer> generateNewSolution (Random random, List<Integer> currentSolutionID, List<Vertiport> vertiportsCandidates, boolean considerCarbon){
        List<Integer> newSolutionID=new ArrayList<>(currentSolutionID);
        List<Integer> notChosenVertiportID = new ArrayList<>();
        List<Integer> SaturationVertiportsID = new ArrayList<>();
        List<Vertiport> newSolution = new ArrayList<>();
        for (Integer vertiportID:newSolutionID){
            newSolution.add(vertiportsCandidates.get(vertiportID));
    }
        List<Integer> NotSaturationVertiportsID = new ArrayList<>();
        for (Vertiport vertiport:newSolution){
            // get the maxSaturationRate for each vertiport
            double maxSaturationRate=0;
            for (int i=0;i<24;i++){
                if (vertiport.saturationRates.get(i)>vertiport.maxSaturationRate){
                    vertiport.maxSaturationRate=vertiport.saturationRates.get(i);
                }
            }

            if (vertiport.maxSaturationRate>1){
                SaturationVertiportsID.add(vertiport.ID);
            }
            else {
                NotSaturationVertiportsID.add(vertiport.ID);
            }
        }
        for (Vertiport vertiport:vertiportsCandidates){
            if (!currentSolutionID.contains(vertiport.ID)){
                notChosenVertiportID.add(vertiport.ID);
            }
        }
        if (SaturationVertiportsID.size()>0 && considerCarbon)
            // get the vertiport with highest saturationRate in SaturationVertiports
        {
            Integer vertiportWithHighestSaturationRateID=SaturationVertiportsID.get(0);
            for (Integer vertiportID:SaturationVertiportsID){
                if ( vertiportsCandidates.get(vertiportID).maxSaturationRate>vertiportsCandidates.get(vertiportWithHighestSaturationRateID).maxSaturationRate){
                    vertiportWithHighestSaturationRateID=vertiportID;
                }
                }
                Vertiport vertiportWithHighestSaturationRate=vertiportsCandidates.get(vertiportWithHighestSaturationRateID);
            // randomly select a vertiport from the neighbors of the vertiportWithHighestSaturationRate and make sure it is not in the current solution
                List<Vertiport> neighbors=vertiportWithHighestSaturationRate.neighbors;
                List<Integer> neighborsID=new ArrayList<>();
                for (Vertiport vertiport:neighbors){
                    neighborsID.add(vertiport.ID);
                }
            // if the neighbors are all in the current solution, remove the vertiportWithHighestSaturationRate from the SaturationVertiports, repeat the process
            while (currentSolutionID.containsAll(neighborsID)||neighborsID.size()==0){
                SaturationVertiportsID.remove(vertiportWithHighestSaturationRateID);
                if ( SaturationVertiportsID.size()>0){
                    vertiportWithHighestSaturationRateID=SaturationVertiportsID.get(0);
                    for (Integer vertiportID:SaturationVertiportsID){
                        if ( vertiportsCandidates.get(vertiportID).maxSaturationRate>vertiportsCandidates.get(vertiportWithHighestSaturationRateID).maxSaturationRate){
                            vertiportWithHighestSaturationRateID=vertiportID;
                        }
                    }
                    vertiportWithHighestSaturationRate=vertiportsCandidates.get(vertiportWithHighestSaturationRateID);
                    neighbors=vertiportWithHighestSaturationRate.neighbors;
                    neighborsID=new ArrayList<>();
                    for (Vertiport vertiport:neighbors){
                        neighborsID.add(vertiport.ID);
                    }
                }
                else {
                    // randomly select a vertiport from the notChosenVertiport and remove a vertiport from the current solution
                    Integer newVertiportID=notChosenVertiportID.get(random.nextInt(notChosenVertiportID.size()));
                    Integer removedVertiportID=currentSolutionID.get(random.nextInt(currentSolutionID.size()));
                    newSolutionID.remove(removedVertiportID);
                    newSolutionID.add(newVertiportID);
                    return newSolutionID;
                }
            }
            Integer newVertiportID=neighborsID.get(random.nextInt(neighborsID.size()));
            while (currentSolutionID.contains(newVertiportID)){
                newVertiportID=neighborsID.get(random.nextInt(neighborsID.size()));
        }
            Integer removedVertiportID= NotSaturationVertiportsID.get(random.nextInt(NotSaturationVertiportsID.size()));
            newSolutionID.remove(removedVertiportID);
            newSolutionID.add(newVertiportID);
        }
        else {
            // randomly select a vertiport from the notChosenVertiport and remove a vertiport from the current solution
            Integer newVertiportID=notChosenVertiportID.get(random.nextInt(notChosenVertiportID.size()));
            Integer removedVertiportID=currentSolutionID.get(random.nextInt(currentSolutionID.size()));
            newSolutionID.remove(removedVertiportID);
            newSolutionID.add(newVertiportID);
        }
    return newSolutionID;
}
    public static double calculateEuciDistance(Coord coord1, Coord coord2) {
        double euciDistance = Math.sqrt(Math.pow(coord1.getX() - coord2.getX(), 2) + Math.pow(coord1.getY() - coord2.getY(), 2));
        return euciDistance;
    }
    public static List<TripItemForOptimization> deserializeTripItems(String fileName) {
        List<TripItemForOptimization> tripItemForOptimizations = new ArrayList<>();

        try  {
            FileInputStream fileIn = new FileInputStream(fileName);
            ObjectInputStream objectIn = new ObjectInputStream(fileIn)  ;
            tripItemForOptimizations = (List<TripItemForOptimization>) objectIn.readObject();


        } catch (Exception e) {
            e.printStackTrace();
        }

        return tripItemForOptimizations;
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
