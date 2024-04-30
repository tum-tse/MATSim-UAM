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
    public static final boolean CONSIDER_CARBON= true; // Euro/kgCO2
    public static void main(String[] args) throws IOException {
        MemoryObserver.start(600);
        // Provide the file via program arguments
        if (args.length > 0) {
            serializedTripItemFile = args[0];
            vertiportCandidateFile = args[1];
        }
        // Get the object TripItem from the serialized file

        log.info("Loading the vertiport candidates...");
        VertiportReader vertiportReader = new VertiportReader();
        List<Vertiport> vertiportsCandidates = vertiportReader.getVertiportsWithNeighbors(vertiportCandidateFile);
        log.info("Finished loading the vertiport candidates.");
        log.info("Loading the trips...");
        List<TripItem> deserializedTripItems = deserializeTripItems(serializedTripItemFile);
        log.info("Finished loading the trips.");
        Random random = new Random(42);
        List<Vertiport> currentSolution = generateRandomSolution(random,vertiportsCandidates,NUM_OF_SELECTED_VERTIPORTS);
        double currentEnergey = calculateFitness(currentSolution,vertiportsCandidates,deserializedTripItems);
        List<Vertiport> bestSolution = new ArrayList<>(currentSolution);
        double bestEnergy= currentEnergey;
        double currentTemperature = INITIAL_TEMPERATURE;
        int notChangeCount=0;
        for (int iteration=0;iteration<10000;iteration++ ){
            List<Vertiport> newSolution= generateNewSolution(random,currentSolution,vertiportsCandidates);
            // set the saturation rate of all vertiports to 0
            for (Vertiport vertiport:vertiportsCandidates){
                for (int i=0;i<24;i++){
                    vertiport.saturationRates.put(i,0.0);
                }
            }
            double newEnergy= calculateFitness(newSolution,vertiportsCandidates,deserializedTripItems);
            double deltaEnergy= newEnergy-currentEnergey;
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

    public static List<Vertiport> generateRandomSolution(Random random,List<Vertiport> VertiportCandidates, int numOfSelectedVertiports) {

      List<Vertiport> selectedVertiports = new ArrayList<>();
      while (selectedVertiports.size() < numOfSelectedVertiports) {
          Vertiport vertiport = VertiportCandidates.get(random.nextInt(VertiportCandidates.size()));
          if (!selectedVertiports.contains(vertiport)) {
              selectedVertiports.add(vertiport);
          }
      }
      return selectedVertiports;
    }

    public static double calculateFitness(List<Vertiport> chosenVertiport,List<Vertiport> vertiportCandidates, List<TripItem> deserializedTripItems) {
        // 实现适应度函数的具体逻辑

        double sumVertiportConstructionCost=0.0;
        double savedGeneralizedCost=0.0;

        List<TripItem> uamAvailableTrips = new ArrayList<>();
//        // calculate the sum of vertiport construction cost
//        for (Integer vertiportID:chosenVertiportID){
//            sumVertiportConstructionCost=sumVertiportConstructionCost+vertiportCandidates.get(vertiportID).constructionCost;
//        }

        for (TripItem tripItem : deserializedTripItems)
        {
            tripItem.isUAMAvailable = false;
            tripItem.uamUtility=-9999;
            List<Vertiport> originNeighbourVertiports = findAvailableNeighbourVertiports(chosenVertiport, tripItem.originNeighborVertiportCandidates);
            List<Vertiport> destinationNeighbourVertiports = findAvailableNeighbourVertiports(chosenVertiport,tripItem.destinationNeighborVertiportCandidates);
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
            double objectiveFunctionBefore;
            double objectiveFunctionAfter;
            double depatureTime=tripItem.departureTime;
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
                        double accessEmisson =0;
                        double egressEmisson =0;
                        double flightDistance= calculateEuciDistance(origin.coord,destination.coord);
                        double flightTime=flightDistance/flightSpeed+takeOffLandingTime;
                        double flightCost=6.1+ flightDistance/1000*0.6;
                        double flightEmission=flightDistance/1000*UAM_EMISSION_FACTOR;
                        double uamTravelTime=accessTime+egressTime+flightTime+UAM_PROCESS_TIME;

                        if (tripItem.accessMode.equals("car") ){
                            accessCost=accessDistance/1000*0.42;
                            accessEmisson=accessDistance/1000*CAR_EMISSION_FACTOR;
                        }
                        if (tripItem.egressMode.equals("car") ){
                            egressCost=egressDistance/1000*0.42;
                            egressEmisson=egressDistance/1000*CAR_EMISSION_FACTOR;
                        }
                        double UAMCost=accessCost+egressCost+flightCost;
                        double UAMGeneralizedCost=UAMCost+uamTravelTime*tripItem.VOT;
                        if (CONSIDER_CARBON){
                            UAMGeneralizedCost=UAMGeneralizedCost+CARBON_EQUIVALENCE_FACTOR*(accessEmisson+egressEmisson+flightEmission);
                        }
                        if (UAMGeneralizedCost < lowestUAMGeneralizedCost) {
                            tripItem.uamTravelTime=uamTravelTime;
                            tripItem.UAMCost=UAMCost;
                            tripItem.UAMUtilityVar=-2.48*UAMCost/100-4.28*flightTime/6000-6.79*(uamTravelTime-flightTime)/6000;
                            tripItem.uamUtility=tripItem.UAMUtilityFix+tripItem.UAMUtilityVar;
                            tripItem.UAMGeneralizedCost=UAMGeneralizedCost;
                            tripItem.accessVertiport = origin;
                            tripItem.egressVertiport = destination;
                            tripItem.accessTime=accessTime;
                            tripItem.egressTime=egressTime;
                            tripItem.flightTime=flightTime;
                            lowestUAMGeneralizedCost = UAMGeneralizedCost;
                        }
                    }
                }
            }
            // determine the probability of mode choice of each trip
            tripItem.uamProbability=calculateModeProbability(tripItem.uamUtility,tripItem.carUtility,tripItem.ptUtility).get(0);
            tripItem.carProbability=calculateModeProbability(tripItem.uamUtility,tripItem.carUtility,tripItem.ptUtility).get(1);
            tripItem.ptProbability=calculateModeProbability(tripItem.uamUtility,tripItem.carUtility,tripItem.ptUtility).get(2);
            Double [] modeProbability={tripItem.carProbability,tripItem.ptProbability,tripItem.uamProbability};
            // update the saturation rate of the access and egress vertiport

            int arriveVertiportHour=(int)Math.floor((tripItem.departureTime+tripItem.accessTime)/3600);
            int leaveVertiportHour=(int)Math.floor((tripItem.departureTime+tripItem.accessTime+UAM_PROCESS_TIME+tripItem.flightTime)/3600);
            tripItem.accessVertiport.saturationRates.put(arriveVertiportHour,tripItem.accessVertiport.saturationRates.get(arriveVertiportHour)+tripItem.uamProbability/tripItem.accessVertiport.capacity);
            tripItem.egressVertiport.saturationRates.put(leaveVertiportHour,tripItem.egressVertiport.saturationRates.get(leaveVertiportHour)+tripItem.uamProbability/tripItem.egressVertiport.capacity);

            ModeDecider modeDecider=new ModeDecider(modeProbability);
            Double [] modeSamples=modeDecider.sample(100);
            if (!CONSIDER_CARBON)
            {objectiveFunctionBefore= tripItem.currentGeneralizedCost+tripItem.currentEmission*CARBON_EQUIVALENCE_FACTOR;
            objectiveFunctionAfter= modeSamples[0]*tripItem.carGeneralizedCost+modeSamples[1]*tripItem.ptGeneralizedCost+modeSamples[2]*tripItem.UAMGeneralizedCost;
            }
            else {
                objectiveFunctionBefore= tripItem.currentGeneralizedCost;
                objectiveFunctionAfter= modeSamples[0]*(tripItem.carGeneralizedCost+tripItem.carEmission*CARBON_EQUIVALENCE_FACTOR)+modeSamples[1]*(tripItem.ptGeneralizedCost+tripItem.ptEmission*CARBON_EQUIVALENCE_FACTOR)+modeSamples[2]*tripItem.UAMGeneralizedCost;
            }


            double savedGeneralizedCostOneTrip=objectiveFunctionBefore-objectiveFunctionAfter;

            if (savedGeneralizedCostOneTrip<0){
                savedGeneralizedCostOneTrip=0;
            }
            if (tripItem.tripPurpose.startsWith("H")){
                savedGeneralizedCostOneTrip=savedGeneralizedCostOneTrip*2;
            }
            savedGeneralizedCost=savedGeneralizedCost+savedGeneralizedCostOneTrip;




        }


       return savedGeneralizedCost;


}
public static List<Vertiport> generateNewSolution (Random random, List<Vertiport> currentSolution, List<Vertiport> vertiportCandidates){
        List<Vertiport> newSolution=new ArrayList<>(currentSolution);
        List<Vertiport> notChosenVertiport = new ArrayList<>();
        List<Vertiport> SaturationVertiports = new ArrayList<>();
        List<Vertiport> NotSaturationVertiports = new ArrayList<>();
        for (Vertiport vertiport:currentSolution){
            // get the maxSaturationRate for each vertiport
            double maxSaturationRate=0;
            for (int i=0;i<24;i++){
                if (vertiport.saturationRates.get(i)>vertiport.maxSaturationRate){
                    vertiport.maxSaturationRate=vertiport.saturationRates.get(i);
                }
            }

            if (vertiport.maxSaturationRate>1){
                SaturationVertiports.add(vertiport);
            }
            else {
                NotSaturationVertiports.add(vertiport);
            }
        }
        for (Vertiport vertiport:vertiportCandidates){
            if (!currentSolution.contains(vertiport)){
                notChosenVertiport.add(vertiport);
            }
        }
        if (SaturationVertiports.size()>0)
            // get the vertiport with highest saturationRate in SaturationVertiports
        {
            Vertiport vertiportWithHighestSaturationRate=SaturationVertiports.get(0);
            for (Vertiport vertiport:SaturationVertiports){
                if (vertiport.maxSaturationRate>vertiportWithHighestSaturationRate.maxSaturationRate){
                    vertiportWithHighestSaturationRate=vertiport;
                }
            }
            // randomly select a vertiport from the neighbors of the vertiportWithHighestSaturationRate and make sure it is not in the current solution
            List<Vertiport> neighbors=vertiportWithHighestSaturationRate.neighbors;
            // if the neighbors are all in the current solution, remove the vertiportWithHighestSaturationRate from the SaturationVertiports, repeat the process
            while (currentSolution.containsAll(neighbors)){
                SaturationVertiports.remove(vertiportWithHighestSaturationRate);
                if (SaturationVertiports.size()>0){
                    vertiportWithHighestSaturationRate=SaturationVertiports.get(0);
                    for (Vertiport vertiport:SaturationVertiports){
                        if (vertiport.maxSaturationRate>vertiportWithHighestSaturationRate.maxSaturationRate){
                            vertiportWithHighestSaturationRate=vertiport;
                        }
                    }
                    neighbors=vertiportWithHighestSaturationRate.neighbors;
                }
                else {
                    // randomly select a vertiport from the notChosenVertiport and remove a vertiport from the current solution
                    Vertiport newVertiport=notChosenVertiport.get(random.nextInt(notChosenVertiport.size()));
                    Vertiport removedVertiport=currentSolution.get(random.nextInt(currentSolution.size()));
                    newSolution.remove(removedVertiport);
                    newSolution.add(newVertiport);
                    return newSolution;
                }
            }
            Vertiport newVertiport=neighbors.get(random.nextInt(neighbors.size()));
            while (currentSolution.contains(newVertiport)){
                newVertiport=neighbors.get(random.nextInt(neighbors.size()));
        }
            Vertiport removedVertiport= NotSaturationVertiports.get(random.nextInt(NotSaturationVertiports.size()));
            newSolution.remove(removedVertiport);
            newSolution.add(newVertiport);
            return newSolution;
        }
        else {
            // randomly select a vertiport from the notChosenVertiport and remove a vertiport from the current solution
            Vertiport newVertiport=notChosenVertiport.get(random.nextInt(notChosenVertiport.size()));
            Vertiport removedVertiport=currentSolution.get(random.nextInt(currentSolution.size()));
            newSolution.remove(removedVertiport);
            newSolution.add(newVertiport);
            return newSolution;
        }
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
    public static List<Vertiport> findAvailableNeighbourVertiports(List<Vertiport> chosenVertiport, List<Vertiport> neighborVertiportCandidates) {
        List<Vertiport> duplicates = new ArrayList<>();

        for (Vertiport vertiport : neighborVertiportCandidates) {
            if (chosenVertiport.contains(vertiport)) {
                duplicates.add(vertiport);
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
