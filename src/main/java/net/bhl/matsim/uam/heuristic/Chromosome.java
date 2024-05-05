package net.bhl.matsim.uam.heuristic;

import net.bhl.matsim.uam.optimization.utils.TripItemForOptimization;
import net.bhl.matsim.uam.optimization.Vertiport;

import java.util.ArrayList;
import java.util.List;

import static net.bhl.matsim.uam.optimization.VertiportOptimizerGenetic.*;

public class Chromosome {

        public List<Integer> genes;
        public double fitness;
public Chromosome(List<Integer> genes) {
            this.genes = genes;
        }

public List<Integer> encodeChromosome(){
    List<Integer> selectedVertiportID = new ArrayList<>();
    for (int i=0;i<this.genes.size();i++){
        if (this.genes.get(i)==1){
            selectedVertiportID.add(i);
        }
    }
    return selectedVertiportID;
}

    public void calculateFitness(List<Vertiport> vertiportCandidates, List<TripItemForOptimization> deserializedTripItemForOptimizations) {
        // 实现适应度函数的具体逻辑
        double sumGeneralizedCost=0.0;
        double sumVertiportConstructionCost=0.0;
        double savedGeneralizedCost=0.0;
        List<Integer> chosenVertiportID=encodeChromosome();
        List<TripItemForOptimization> uamAvailableTrips = new ArrayList<>();
        // calculate the sum of vertiport construction cost
        for (Integer vertiportID:chosenVertiportID){
            sumVertiportConstructionCost=sumVertiportConstructionCost+vertiportCandidates.get(vertiportID).constructionCost;
        }

        for (TripItemForOptimization tripItemForOptimization : deserializedTripItemForOptimizations)
        {
            tripItemForOptimization.isUAMAvailable = false;
            tripItemForOptimization.uamUtility=-9999;
            List<Vertiport> originNeighbourVertiports = findAvailableNeighbourVertiports(chosenVertiportID, tripItemForOptimization.originNeighborVertiportCandidates);
            List<Vertiport> destinationNeighbourVertiports = findAvailableNeighbourVertiports(chosenVertiportID, tripItemForOptimization.destinationNeighborVertiportCandidates);
            if (originNeighbourVertiports.size() > 0 && destinationNeighbourVertiports.size() > 0) {
                tripItemForOptimization.isUAMAvailable = true;
                tripItemForOptimization.originNeighborVertiports = originNeighbourVertiports;
                tripItemForOptimization.destinationNeighborVertiports = destinationNeighbourVertiports;
                uamAvailableTrips.add(tripItemForOptimization);}
        }

        for (TripItemForOptimization tripItemForOptimization : uamAvailableTrips) {

            for (Vertiport vertiport : tripItemForOptimization.originNeighborVertiports) {
                tripItemForOptimization.originNeighborVertiportsTimeAndDistance.put(vertiport, tripItemForOptimization.originNeighborVertiportCandidatesTimeAndDistance.get(vertiport));
            }
            for (Vertiport vertiport : tripItemForOptimization.destinationNeighborVertiports) {
                tripItemForOptimization.destinationNeighborVertiportsTimeAndDistance.put(vertiport, tripItemForOptimization.destinationNeighborVertiportCandidatesTimeAndDistance.get(vertiport));
            }

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
                        double flightDistance= calculateEuciDistance(origin.coord,destination.coord);
                        double flightTime=flightDistance/flightSpeed+takeOffLandingTime;
                        double flightCost=6.1+ calculateEuciDistance(origin.coord,destination.coord)/1000*0.6;
                        double uamTravelTime=accessTime+egressTime+flightTime+UAM_PROCESS_TIME;
                        if (tripItemForOptimization.accessMode.equals("car") ){
                            accessCost=accessDistance/1000*0.42;
                        }
                        if (tripItemForOptimization.egressMode.equals("car") ){
                            egressCost=egressDistance/1000*0.42;
                        }
                        double UAMCost=accessCost+egressCost+flightCost;
                        double UAMGeneralizedCost=UAMCost+uamTravelTime* tripItemForOptimization.VOT;
                        if (UAMGeneralizedCost < lowestUAMGeneralizedCost) {
                            tripItemForOptimization.uamTravelTime=uamTravelTime;
                            tripItemForOptimization.UAMCost=UAMCost;
                            tripItemForOptimization.UAMUtilityVar=-2.48*UAMCost/100-4.28*flightTime/6000-6.79*(uamTravelTime-flightTime)/6000;
                            tripItemForOptimization.uamUtility= tripItemForOptimization.UAMUtilityFix+ tripItemForOptimization.UAMUtilityVar;
                            tripItemForOptimization.UAMGeneralizedCost=UAMGeneralizedCost;
                            tripItemForOptimization.accessVertiport = origin;
                            tripItemForOptimization.egressVertiport = destination;
                            lowestUAMGeneralizedCost = UAMGeneralizedCost;
                        }
                    }
                }
            }
            // determine the probability of mode choice of each trip
            tripItemForOptimization.uamProbability=calculateModeProbability(tripItemForOptimization.uamUtility, tripItemForOptimization.carUtility, tripItemForOptimization.ptUtility).get(0);
            tripItemForOptimization.carProbability=calculateModeProbability(tripItemForOptimization.uamUtility, tripItemForOptimization.carUtility, tripItemForOptimization.ptUtility).get(1);
            tripItemForOptimization.ptProbability=calculateModeProbability(tripItemForOptimization.uamUtility, tripItemForOptimization.carUtility, tripItemForOptimization.ptUtility).get(2);
            double generalizedCostOneTripBefore= tripItemForOptimization.carGeneralizedCost*calculateModeProbability(-9999, tripItemForOptimization.carUtility, tripItemForOptimization.ptUtility).get(1)+ tripItemForOptimization.ptGeneralizedCost*calculateModeProbability(-9999, tripItemForOptimization.carUtility, tripItemForOptimization.ptUtility).get(2);
            double generalizedCostOneTripAfter= tripItemForOptimization.UAMGeneralizedCost* tripItemForOptimization.uamProbability+ tripItemForOptimization.carGeneralizedCost* tripItemForOptimization.carProbability+ tripItemForOptimization.ptGeneralizedCost* tripItemForOptimization.ptProbability;
            double savedGeneralizedCostOneTrip=generalizedCostOneTripBefore-generalizedCostOneTripAfter;
            if (savedGeneralizedCostOneTrip<0){
                savedGeneralizedCostOneTrip=0;
            }
            if (tripItemForOptimization.tripPurpose.startsWith("H")){
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




        this.fitness=savedGeneralizedCost;

    }
        }

