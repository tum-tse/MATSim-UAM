package net.bhl.matsim.uam.optimization;

import net.bhl.matsim.uam.optimization.utils.ScenarioSpecific;
import net.bhl.matsim.uam.optimization.utils.TripItemForOptimization;
import org.apache.log4j.Logger;
import org.matsim.api.core.v01.Coord;
import org.matsim.utils.MemoryObserver;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class SimulatedAnnealingForPartA {
    // This class provides the simulated annealing algorithm to solve the vertiport siting problem.

    public static final Logger log = Logger.getLogger(SimulatedAnnealingForPartA.class);

    public static double INITIAL_TEMPERATURE;
    public static double FINAL_TEMPERATURE;
    public static double ANNEALING_RATE;
    public static double MAX_ITERATION;
    public static double MAX_NOT_CHANGE_COUNT;
    public static String serializedTripItemFile;
    public static String vertiportCandidateFile;
    public static double flightSpeed; // m/s
    public static double UAM_PROCESS_TIME; // s
    public static double takeOffLandingTime; // s
    public static int NUM_OF_SELECTED_VERTIPORTS;
    private static final int MEMORY_CHECK_INTERVAL=600;
    private static double UAM_FIX_COST;
    private static double UAM_KM_COST;
    private static double CAR_COST;
    private static boolean considerReturnTrip;
    private static  int sampleSize;
    private static long[] RANDOM_SEEDS;
    private static int num_of_run;
    public static double CAR_EMISSION_FACTOR; // kg/km
    public static double PT_EMISSION_FACTOR; // kg/km
    public static double UAM_EMISSION_FACTOR; // kg/km
    public static double CARBON_EQUIVALENCE_FACTOR; // Euro/kgCO2
    public static boolean CONSIDER_CARBON; // Euro/kgCO2
    private static int SIMULATION_HOURS;
    public static String scenarioName;
    // Number of parallel tasks

    public static void main(String[] args) throws IOException, InterruptedException {
        MemoryObserver.start(MEMORY_CHECK_INTERVAL);
        // Provide the file via program arguments
        if (args.length > 0) {
            serializedTripItemFile = args[0];
            vertiportCandidateFile = args[1];
            sampleSize = Integer.parseInt(args[2]);
            num_of_run = Integer.parseInt(args[3]);
            scenarioName = args[4];
            RANDOM_SEEDS = new long[num_of_run];
            for (int i = 0; i < num_of_run; i++) {
                RANDOM_SEEDS[i] = Long.parseLong(args[5 + i]);
            }
        }
        // Build the scenario of Munich
        ScenarioSpecific scenarioSpecific = new ScenarioSpecific(scenarioName);
        scenarioSpecific.buildScenario();
        NUM_OF_SELECTED_VERTIPORTS = scenarioSpecific.num_of_selected_vertiports;
        UAM_FIX_COST = scenarioSpecific.uam_fix_cost;
        UAM_KM_COST = scenarioSpecific.uam_km_cost;
        CAR_COST = scenarioSpecific.car_km_cost;
        flightSpeed = scenarioSpecific.flight_speed;
        UAM_PROCESS_TIME = scenarioSpecific.uam_process_time;
        takeOffLandingTime = scenarioSpecific.uam_take_off_landing_time;
        considerReturnTrip = scenarioSpecific.consider_return_trip;
        CAR_EMISSION_FACTOR = scenarioSpecific.car_emission_factor;


        // Get the object TripItemForOptimization from the serialized file

        log.info("Loading the vertiport candidates...");
        VertiportReader vertiportReader = new VertiportReader();
        List<Vertiport> vertiportsCandidates = vertiportReader.getVertiportsWithNeighbors(vertiportCandidateFile);
        log.info("Finished loading the vertiport candidates.");
        log.info("Loading the trips...");
        List<TripItemForOptimization> deserializedTripItemForOptimizations = deserializeTripItems(serializedTripItemFile);
        log.info("Finished loading the trips.");
        // Initialize the UAM enabled trips and UAM utility
        log.info("Start the simulated annealing algorithm...");


        // run the simulated annealing 10 times
        for (int index_of_run=0;index_of_run<num_of_run;index_of_run++) {

            for (TripItemForOptimization tripItemForOptimization : deserializedTripItemForOptimizations) {
                tripItemForOptimization.isUAMAvailable = false;
                tripItemForOptimization.uamUtility = Double.NEGATIVE_INFINITY;
                tripItemForOptimization.savedGeneralizedCost = 0;
            }

            log.info("Run: " + index_of_run+1);
            // record the start time
            long startTime = System.currentTimeMillis();
            Random random = new Random(RANDOM_SEEDS[index_of_run]);
            List<Integer> currentSolutionID = generateRandomSolution(random, vertiportsCandidates, NUM_OF_SELECTED_VERTIPORTS);
            double currentEnergey = calculateFitness(currentSolutionID,  vertiportsCandidates,currentSolutionID, deserializedTripItemForOptimizations,random,scenarioSpecific);
            for (TripItemForOptimization tripItemForOptimization : deserializedTripItemForOptimizations) {
                if (!tripItemForOptimization.tempSavedGeneralizedCosts.isEmpty()) {
                    tripItemForOptimization.savedGeneralizedCost = tripItemForOptimization.tempSavedGeneralizedCosts.get(0);
                }
            }
            List<Integer> bestSolutionID = new ArrayList<>(currentSolutionID);
            double bestEnergy = currentEnergey;
            double currentTemperature = INITIAL_TEMPERATURE;
            int notChangeCount = 0;
            int saturatedVertiportCount = 0;
            double maxSaturationRate = 0;

            if (CONSIDER_CARBON) {
                for (Integer vertiportID : currentSolutionID) {
                    if (vertiportsCandidates.get(vertiportID).maxSaturationRate > 1) {
                        saturatedVertiportCount++;
                    }
                    if (vertiportsCandidates.get(vertiportID).maxSaturationRate > maxSaturationRate) {
                        maxSaturationRate = vertiportsCandidates.get(vertiportID).maxSaturationRate;
                    }
                }
            }
            for (int iteration = 0; iteration < MAX_ITERATION; iteration++) {
                List<List<Integer>> newSolutionList = generateNewSolution(random, currentSolutionID, vertiportsCandidates);
                List<Integer> newSolutionID = newSolutionList.get(0);
                List<Integer> selectionDifference = newSolutionList.get(1);

                // set the saturation rate of all vertiports to 0
                if (CONSIDER_CARBON) {
                    for (Vertiport vertiport : vertiportsCandidates) {
                        for (int i = 0; i < SIMULATION_HOURS; i++) {
                            vertiport.saturationRates.put(i, 0.0);
                        }
                    }
                }
                double newEnergy = calculateFitness(newSolutionID, vertiportsCandidates,selectionDifference ,deserializedTripItemForOptimizations,random,scenarioSpecific);
                double deltaEnergy = newEnergy - currentEnergey;

                if (deltaEnergy > 0 || Math.exp(deltaEnergy / currentTemperature) > random.nextDouble()) {
                    currentSolutionID = new ArrayList<>(newSolutionID);
                    currentEnergey = newEnergy;
                    // update the savedGeneralizedCost for each trip
                    for (TripItemForOptimization tripItemForOptimization : deserializedTripItemForOptimizations) {
                        if (!tripItemForOptimization.tempSavedGeneralizedCosts.isEmpty()) {
                            tripItemForOptimization.savedGeneralizedCost = tripItemForOptimization.tempSavedGeneralizedCosts.get(0);
                        }
                    }
                    if (CONSIDER_CARBON) {
                        for (Integer vertiportID : currentSolutionID) {
                            if (vertiportsCandidates.get(vertiportID).maxSaturationRate > 1) {
                                saturatedVertiportCount++;
                            }
                            if (vertiportsCandidates.get(vertiportID).maxSaturationRate > maxSaturationRate) {
                                maxSaturationRate = vertiportsCandidates.get(vertiportID).maxSaturationRate;
                            }
                        }
                    }
                }
                // update the best solution
                if (currentEnergey > bestEnergy) {
                    bestSolutionID = new ArrayList<>(currentSolutionID);
                    bestEnergy = currentEnergey;
                    notChangeCount = 0;
                } else {
                    notChangeCount++;
                }


                // update the temperature
                if (currentTemperature > FINAL_TEMPERATURE) {
                    currentTemperature = currentTemperature * ANNEALING_RATE;
                }
                if (CONSIDER_CARBON) {
                    log.info("Iteration: " + iteration + " Current Temperature: " + currentTemperature + " Current Energy: " + currentEnergey + " Best Energy: " + bestEnergy + " Saturated Vertiport Count: " + saturatedVertiportCount + " Max Saturation Rate: " + maxSaturationRate);
                } else {
                    log.info("Iteration: " + iteration + " Current Temperature: " + currentTemperature + " Current Energy: " + currentEnergey + " Best Energy: " + bestEnergy);
                }
                // if the best solution is not updated for more than 1000 iterations, break
                if (notChangeCount > MAX_NOT_CHANGE_COUNT) {
                    break;
                }
                saturatedVertiportCount = 0;
                maxSaturationRate = 0;
            }

            log.info("Best Solution: " + bestSolutionID + " Best Energy: " + bestEnergy);
            // record the end time
            long endTime = System.currentTimeMillis();
            log.info("Run: " + index_of_run+1 + " Time: " + (endTime - startTime) / 1000 + "s");
        }
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

    public static double calculateFitness(List<Integer> chosenVertiportID, List<Vertiport> vertiportsCandidates, List<Integer> selectionDifference, List<TripItemForOptimization> deserializedTripItemForOptimizations, Random random, ScenarioSpecific scenarioSpecific) {
        List<TripItemForOptimization> uamAvailableTrips = new ArrayList<>();
        List<TripItemForOptimization> differenceTrips = new ArrayList<>();
        double totalSavedGeneralizedCost = 0;
        double totalConstructionCost = 0;
        double totalObjectiveFunctionValue = 0;

        for (TripItemForOptimization tripItemForOptimization : deserializedTripItemForOptimizations) {
            tripItemForOptimization.tempSavedGeneralizedCosts.clear();
            boolean added = false;
            for (Vertiport vertiport : tripItemForOptimization.originNeighborVertiportCandidates) {
                if (selectionDifference.contains(vertiport.ID)) {
                    differenceTrips.add(tripItemForOptimization);
                    added = true;
                    break;
                }
            }
            if (added) {
                continue;
            }
            for (Vertiport vertiport : tripItemForOptimization.destinationNeighborVertiportCandidates) {
                if (selectionDifference.contains(vertiport.ID)) {
                    differenceTrips.add(tripItemForOptimization);
                    break;
                }
            }
        }

        for (TripItemForOptimization tripItemForOptimization : differenceTrips) {
            List<Vertiport> originNeighbourVertiports = findAvailableNeighbourVertiports(chosenVertiportID, tripItemForOptimization.originNeighborVertiportCandidates);
            List<Vertiport> destinationNeighbourVertiports = findAvailableNeighbourVertiports(chosenVertiportID, tripItemForOptimization.destinationNeighborVertiportCandidates);
            if (!originNeighbourVertiports.isEmpty() && !destinationNeighbourVertiports.isEmpty()) {
                tripItemForOptimization.isUAMAvailable = true;
                tripItemForOptimization.originNeighborVertiports = originNeighbourVertiports;
                tripItemForOptimization.destinationNeighborVertiports = destinationNeighbourVertiports;
                uamAvailableTrips.add(tripItemForOptimization);
            } else {
                tripItemForOptimization.isUAMAvailable = false;
                tripItemForOptimization.UAMUtilityVar = Double.NEGATIVE_INFINITY;
                tripItemForOptimization.tempSavedGeneralizedCosts.add(0.0);
            }
        }



        for (TripItemForOptimization tripItemForOptimization : uamAvailableTrips) {
            calculateTripSavedCost(tripItemForOptimization, vertiportsCandidates, random, scenarioSpecific);
        }
        for (TripItemForOptimization tripItemForOptimization : deserializedTripItemForOptimizations) {
            if (differenceTrips.contains(tripItemForOptimization)) {
                totalSavedGeneralizedCost += tripItemForOptimization.tempSavedGeneralizedCosts.get(0);
            }
            else {
                totalSavedGeneralizedCost += tripItemForOptimization.savedGeneralizedCost;
            }
        }

        if (CONSIDER_CARBON) {
            for (int vertiportID:chosenVertiportID){
                totalConstructionCost+=vertiportsCandidates.get(vertiportID).constructionCost;
            }
            totalObjectiveFunctionValue=scenarioSpecific.beta_savedCost*totalSavedGeneralizedCost+scenarioSpecific.beta_constructionCost*totalConstructionCost;
        }
        else {
            totalObjectiveFunctionValue=totalSavedGeneralizedCost;
        }
        return totalObjectiveFunctionValue;
    }

    public static void calculateTripSavedCost(TripItemForOptimization tripItemForOptimization, List<Vertiport> vertiportsCandidates, Random random,ScenarioSpecific scenarioSpecific) {
        for (Vertiport vertiport : tripItemForOptimization.originNeighborVertiports) {
            tripItemForOptimization.originNeighborVertiportsTimeAndDistance.put(vertiport, tripItemForOptimization.originNeighborVertiportCandidatesTimeAndDistance.get(vertiport));
        }
        for (Vertiport vertiport : tripItemForOptimization.destinationNeighborVertiports) {
            tripItemForOptimization.destinationNeighborVertiportsTimeAndDistance.put(vertiport, tripItemForOptimization.destinationNeighborVertiportCandidatesTimeAndDistance.get(vertiport));
        }
        double objectiveFunctionBefore;
        double objectiveFunctionAfter;
        // Initialize the Vertiport Allocation
        // Find the vertiport pair for a trip with the lowest uam generalized cost, the access and ergress vertiport could not be the same
        double lowestUAMGeneralizedCost = Double.MAX_VALUE;
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
                    double flightEmission=flightDistance/1000*UAM_EMISSION_FACTOR; // convert m to km
                    double flightCost=UAM_FIX_COST+ calculateEuciDistance(origin.coord,destination.coord)/1000*UAM_KM_COST;
                    double uamTravelTime=accessTime+egressTime+flightTime+UAM_PROCESS_TIME;
                    if (tripItemForOptimization.accessMode.equals("car") ){
                        accessCost=accessDistance/1000*CAR_COST;
                    }
                    if (tripItemForOptimization.egressMode.equals("car") ){
                        egressCost=egressDistance/1000*CAR_COST;
                    }
                    double UAMCost=accessCost+egressCost+flightCost;
                    double UAMGeneralizedCost=UAMCost+uamTravelTime* tripItemForOptimization.VOT;
                    if (CONSIDER_CARBON){
                        UAMGeneralizedCost=UAMGeneralizedCost+CARBON_EQUIVALENCE_FACTOR*(accessEmisson+egressEmisson+flightEmission);
                    }
                    if (UAMGeneralizedCost < lowestUAMGeneralizedCost) {
                        tripItemForOptimization.uamTravelTime=uamTravelTime;
                        tripItemForOptimization.UAMCost=UAMCost;
                        tripItemForOptimization.UAMGeneralizedCost=UAMGeneralizedCost;
                        tripItemForOptimization.UAMUtilityVar= scenarioSpecific.calculateUAMUtilityVAR(flightTime,uamTravelTime-flightTime,UAMCost);
                        tripItemForOptimization.uamUtility= tripItemForOptimization.UAMUtilityFix+ tripItemForOptimization.UAMUtilityVar;

                        tripItemForOptimization.accessVertiport = origin;
                        tripItemForOptimization.egressVertiport = destination;
                        lowestUAMGeneralizedCost = UAMGeneralizedCost;
                    }
                }
            }
        }
        // determine the probability of mode choice of each trip
        if (CONSIDER_CARBON) {
            int arriveVertiportHour = (int) Math.floor((tripItemForOptimization.departureTime + tripItemForOptimization.accessTime) / 3600);
            int leaveVertiportHour = (int) Math.floor((tripItemForOptimization.departureTime + tripItemForOptimization.accessTime + UAM_PROCESS_TIME + tripItemForOptimization.flightTime) / 3600);
            vertiportsCandidates.get(tripItemForOptimization.accessVertiport.ID).saturationRates.put(arriveVertiportHour, vertiportsCandidates.get(tripItemForOptimization.accessVertiport.ID).saturationRates.get(arriveVertiportHour) + tripItemForOptimization.uamProbability / tripItemForOptimization.accessVertiport.totalCapacity);
            vertiportsCandidates.get(tripItemForOptimization.egressVertiport.ID).saturationRates.put(leaveVertiportHour, vertiportsCandidates.get(tripItemForOptimization.egressVertiport.ID).saturationRates.get(leaveVertiportHour) + tripItemForOptimization.uamProbability / tripItemForOptimization.egressVertiport.totalCapacity);
        }
        // Create a double list to store the probability of each mode
        ModeDecider modeDecider=new ModeDecider(tripItemForOptimization.uamUtility,tripItemForOptimization.carUtility,tripItemForOptimization.ptUtility,random);
        Double [] modeSamples=modeDecider.sample(sampleSize);
        if (!CONSIDER_CARBON)
        {objectiveFunctionBefore= tripItemForOptimization.currentGeneralizedCost;
            objectiveFunctionAfter= modeSamples[1]* tripItemForOptimization.carGeneralizedCost+modeSamples[2]* tripItemForOptimization.ptGeneralizedCost+modeSamples[0]* tripItemForOptimization.UAMGeneralizedCost;
        }
        else {
            objectiveFunctionBefore= tripItemForOptimization.currentGeneralizedCost+tripItemForOptimization.currentEmission*CARBON_EQUIVALENCE_FACTOR;
            objectiveFunctionAfter= modeSamples[1]*(tripItemForOptimization.carGeneralizedCost+ tripItemForOptimization.carEmission*CARBON_EQUIVALENCE_FACTOR)+modeSamples[2]*(tripItemForOptimization.ptGeneralizedCost+ tripItemForOptimization.ptEmission*CARBON_EQUIVALENCE_FACTOR)+modeSamples[0]* tripItemForOptimization.UAMGeneralizedCost;
        }
        double savedGeneralizedCostOneTrip=objectiveFunctionBefore-objectiveFunctionAfter;
        if (savedGeneralizedCostOneTrip<0){
            savedGeneralizedCostOneTrip=0;
        }
        if (tripItemForOptimization.tripPurpose.startsWith("H") && considerReturnTrip){
            savedGeneralizedCostOneTrip=savedGeneralizedCostOneTrip*2;
        }

        // Include your logic here that was previously inside the loop over all trips
        tripItemForOptimization.tempSavedGeneralizedCosts.add(savedGeneralizedCostOneTrip);
    }

public static List<List<Integer>> generateNewSolution (Random random, List<Integer> currentSolutionID, List<Vertiport> vertiportsCandidates) {
    List<Integer> newSolutionID = new ArrayList<>(currentSolutionID);
    List<Integer> differenceID = new ArrayList<>();
    List<Integer> notChosenVertiportID = new ArrayList<>();
    List<Integer> SaturationVertiportsID = new ArrayList<>();
    List<List<Integer>> result = new ArrayList<>();
    List<Vertiport> newSolution = new ArrayList<>();

    List<Integer> NotSaturationVertiportsID = new ArrayList<>();
    if (CONSIDER_CARBON) {
        for (Integer vertiportID : newSolutionID) {
            newSolution.add(vertiportsCandidates.get(vertiportID));
        }
        for (Vertiport vertiport : newSolution) {
            // get the maxSaturationRate for each vertiport
            vertiport.maxSaturationRate = 0;
            for (int i = 0; i < SIMULATION_HOURS; i++) {
                if (vertiport.saturationRates.get(i) > vertiport.maxSaturationRate) {
                    vertiport.maxSaturationRate = vertiport.saturationRates.get(i);
                }
            }

            if (vertiport.maxSaturationRate > 1) {
                SaturationVertiportsID.add(vertiport.ID);
            } else {
                NotSaturationVertiportsID.add(vertiport.ID);
            }
        }
        for (Vertiport vertiport : vertiportsCandidates) {
            if (!currentSolutionID.contains(vertiport.ID)) {
                notChosenVertiportID.add(vertiport.ID);
            }
        }
        if (!SaturationVertiportsID.isEmpty()) {
            Integer vertiportWithHighestSaturationRateID = SaturationVertiportsID.get(0);
            for (Integer vertiportID : SaturationVertiportsID) {
                if (vertiportsCandidates.get(vertiportID).maxSaturationRate > vertiportsCandidates.get(vertiportWithHighestSaturationRateID).maxSaturationRate) {
                    vertiportWithHighestSaturationRateID = vertiportID;
                }
            }
            Vertiport vertiportWithHighestSaturationRate = vertiportsCandidates.get(vertiportWithHighestSaturationRateID);
            // randomly select a vertiport from the neighbors of the vertiportWithHighestSaturationRate and make sure it is not in the current solution
            List<Vertiport> neighbors = vertiportWithHighestSaturationRate.neighbors;
            List<Integer> neighborsID = new ArrayList<>();
            for (Vertiport vertiport : neighbors) {
                neighborsID.add(vertiport.ID);
            }
            // if the neighbors are all in the current solution, remove the vertiportWithHighestSaturationRate from the SaturationVertiports, repeat the process
            while (currentSolutionID.containsAll(neighborsID) || neighborsID.isEmpty()) {
                SaturationVertiportsID.remove(vertiportWithHighestSaturationRateID);
                if (!SaturationVertiportsID.isEmpty()) {
                    vertiportWithHighestSaturationRateID = SaturationVertiportsID.get(0);
                    for (Integer vertiportID : SaturationVertiportsID) {
                        if (vertiportsCandidates.get(vertiportID).maxSaturationRate > vertiportsCandidates.get(vertiportWithHighestSaturationRateID).maxSaturationRate) {
                            vertiportWithHighestSaturationRateID = vertiportID;
                        }
                    }
                    vertiportWithHighestSaturationRate = vertiportsCandidates.get(vertiportWithHighestSaturationRateID);
                    neighbors = vertiportWithHighestSaturationRate.neighbors;
                    neighborsID = new ArrayList<>();
                    for (Vertiport vertiport : neighbors) {
                        neighborsID.add(vertiport.ID);
                    }
                }
                else {
                    // randomly select a vertiport from the notChosenVertiport and remove a vertiport from the current solution
                    Integer newVertiportID = notChosenVertiportID.get(random.nextInt(notChosenVertiportID.size()));
                    Integer removedVertiportID = currentSolutionID.get(random.nextInt(currentSolutionID.size()));
                    differenceID.add(removedVertiportID);
                    differenceID.add(newVertiportID);
                    newSolutionID.remove(removedVertiportID);
                    newSolutionID.add(newVertiportID);
                    result.add(newSolutionID);
                    result.add(differenceID);
                    return result;
                }
                Integer newVertiportID = neighborsID.get(random.nextInt(neighborsID.size()));
                while (currentSolutionID.contains(newVertiportID)) {
                    newVertiportID = neighborsID.get(random.nextInt(neighborsID.size()));
                }
                Integer removedVertiportID = NotSaturationVertiportsID.get(random.nextInt(NotSaturationVertiportsID.size()));
                newSolutionID.remove(removedVertiportID);
                newSolutionID.add(newVertiportID);
                differenceID.add(removedVertiportID);
                differenceID.add(newVertiportID);
                result.add(newSolutionID);
                result.add(differenceID);
            }
        }
        else{
            Integer newVertiportID = notChosenVertiportID.get(random.nextInt(notChosenVertiportID.size()));
            Integer removedVertiportID = currentSolutionID.get(random.nextInt(currentSolutionID.size()));
            differenceID.add(removedVertiportID);
            differenceID.add(newVertiportID);
            newSolutionID.remove(removedVertiportID);
            newSolutionID.add(newVertiportID);
            result.add(newSolutionID);
            result.add(differenceID);
        }
    }
    else{
    for (Vertiport vertiport : vertiportsCandidates) {
        if (!currentSolutionID.contains(vertiport.ID)) {
            notChosenVertiportID.add(vertiport.ID);
        }
    }


        // randomly select a vertiport from the notChosenVertiport and remove a vertiport from the current solution
        Integer newVertiportID = notChosenVertiportID.get(random.nextInt(notChosenVertiportID.size()));
        Integer removedVertiportID = currentSolutionID.get(random.nextInt(currentSolutionID.size()));
        differenceID.add(removedVertiportID);
        differenceID.add(newVertiportID);
        newSolutionID.remove(removedVertiportID);
        newSolutionID.add(newVertiportID);
        result.add(newSolutionID);
        result.add(differenceID);
    }
        return result;
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




}
