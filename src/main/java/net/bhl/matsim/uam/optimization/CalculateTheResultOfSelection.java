package net.bhl.matsim.uam.optimization;

import com.opencsv.CSVWriter;
import net.bhl.matsim.uam.analysis.traveltimes.utils.TripItem;
import org.apache.log4j.Logger;
import org.matsim.api.core.v01.Coord;

import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.List;

public class CalculateTheResultOfSelection {
    public static final double flightSpeed= 350/3.6; // m/s
    public static final double UAM_PROCESS_TIME= 10*60; // s
    public static final double takeOffLandingTime= 60; // s
    private static String vertiportCandidateFile;
    private static String servedTripsFileGRD_F;
    private static String servedTripsFileGRD_B;
    private static String servedTripsFileGA;
    private static String servedTripsFileSA;
    private static String servedTripsFileGRD_F_U;
    private static String servedTripsAllCandidates;
    private static int[] GRD_F;
    private static int[] GRD_B;
    private static int[] GA;
    private static int[] SA;
    private static int[] GRD_F_U;
    private static int[] AllCandidates;
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
            servedTripsFileGRD_F = args[2];
            servedTripsFileGRD_B = args[3];
            servedTripsFileGA = args[4];
            servedTripsFileSA = args[5];
            servedTripsFileGRD_F_U = args[6];
            servedTripsAllCandidates = args[7];

        }
        // Allcandidates is 0-199
AllCandidates=new int[200];
        for (int i=0;i<200;i++){
            AllCandidates[i]=i;
        }
        GRD_F = new int[]{6, 150, 124, 154, 94, 156, 53, 16, 135, 159, 188, 31, 74, 82, 133, 76, 39, 166, 44, 65, 151, 174, 22, 24, 7, 122, 87, 179, 152, 111, 199, 55, 198, 175, 2, 20, 168, 46, 119, 88, 141, 52, 45, 189, 48, 192, 176, 49, 158, 84, 50, 32, 100, 109, 169, 21, 17, 99, 81, 61, 178, 145, 92, 149, 132, 42, 144, 153, 193, 3, 93, 41, 160, 165};
        GRD_B= new int[]{150, 94, 135, 6, 87, 74, 122, 16, 53, 46, 188, 39, 119, 166, 154, 124, 99, 24, 175, 84, 82, 49, 141, 17, 192, 50, 168, 179, 155, 158, 30, 57, 169, 65, 55, 48, 149, 116, 13, 151, 105, 64, 51, 26, 156, 176, 133, 132, 45, 2, 109, 41, 31, 44, 22, 178, 145, 143, 165, 153, 196, 42, 197, 189, 193, 159, 104, 100, 76, 144, 7, 88, 152, 32};
        GA= new int[]{6, 7, 16, 17, 19, 20, 23, 24, 26, 30, 31, 32, 36, 39, 41, 44, 45, 46, 48, 49, 50, 51, 53, 57, 61, 64, 65, 74, 78, 79, 80, 87, 88, 93, 94, 96, 99, 100, 104, 105, 107, 109, 116, 119, 122, 126, 132, 133, 137, 138, 141, 143, 144, 146, 149, 150, 152, 154, 155, 156, 158, 159, 160, 164, 165, 169, 175, 176, 178, 179, 186, 190, 192, 193};
        SA= new int[]{111, 164, 80, 39, 107, 104, 152, 3, 32, 74, 36, 84, 149, 159, 132, 192, 150, 65, 176, 158, 188, 57, 165, 145, 30, 64, 88, 51, 175, 160, 138, 20, 186, 133, 44, 169, 99, 93, 52, 48, 105, 196, 122, 146, 16, 119, 53, 46, 96, 31, 178, 24, 61, 155, 168, 7, 154, 6, 126, 17, 79, 23, 87, 179, 78, 144, 26, 100, 116, 141, 41, 156, 50, 143};
        GRD_F_U= new int[]{6, 150, 156, 16, 53, 154, 159, 122, 179, 87, 74, 175, 7, 149, 168, 39, 133, 46, 141, 119, 188, 24, 65, 48, 31, 44, 158, 99, 57, 155, 189, 178, 145, 176, 55, 196, 52, 50, 51, 32, 26, 132, 153, 144, 152, 143, 64, 169, 105, 104, 107, 164, 116, 93, 78, 84, 20, 96, 80, 79, 186, 41, 138, 23, 100, 135, 160, 30, 3, 137, 126, 146, 21, 117};
// cast the int[] to List<Integer>
        List<Integer> AllCandidates_List = new ArrayList<>();
        for (int i = 0; i < AllCandidates.length; i++) {
            AllCandidates_List.add(AllCandidates[i]);
        }
        List<Integer> GRD_F_List = new ArrayList<>();
        for (int i = 0; i < GRD_F.length; i++) {
            GRD_F_List.add(GRD_F[i]);
        }
        List<Integer> GRD_B_List = new ArrayList<>();
        for (int i = 0; i < GRD_B.length; i++) {
            GRD_B_List.add(GRD_B[i]);
        }
        List<Integer> GA_List = new ArrayList<>();
        for (int i = 0; i < GA.length; i++) {
            GA_List.add(GA[i]);
        }
        List<Integer> SA_List = new ArrayList<>();
        for (int i = 0; i < SA.length; i++) {
            SA_List.add(SA[i]);
        }
        List<Integer> GRD_F_U_List = new ArrayList<>();
        for (int i = 0; i < GRD_F_U.length; i++) {
            GRD_F_U_List.add(GRD_F_U[i]);
        }

        // Get the object TripItem from the serialized file

        log.info("Loading the vertiport candidates...");
        VertiportReader vertiportReader = new VertiportReader();
        List<Vertiport> vertiportsCandidates = vertiportReader.getVertiports(vertiportCandidateFile);
        // Test
        // Only the first 50 vertiports are used for testing
        //      vertiportsCandidates = vertiportsCandidates.subList(0, 50);
        log.info("Finished loading the vertiport candidates.");
        log.info("Loading the trips...");
        List<TripItem> deserializedTripItems = deserializeTripItems(fileName);
        log.info("Finished loading the trips.");
        System.out.println(calculateSelectionScore(vertiportsCandidates, AllCandidates_List, deserializedTripItems, servedTripsAllCandidates));
        System.out.println(calculateSelectionScore(vertiportsCandidates, GRD_F_List, deserializedTripItems, servedTripsFileGRD_F));
        System.out.println(calculateSelectionScore(vertiportsCandidates, GRD_B_List, deserializedTripItems, servedTripsFileGRD_B));
        System.out.println(calculateSelectionScore(vertiportsCandidates, GA_List, deserializedTripItems, servedTripsFileGA));
        System.out.println(calculateSelectionScore(vertiportsCandidates, SA_List, deserializedTripItems, servedTripsFileSA));
        System.out.println(calculateSelectionScore(vertiportsCandidates, GRD_F_U_List, deserializedTripItems, servedTripsFileGRD_F_U));
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

    public static double calculateSelectionScore(List<Vertiport> vertiportCandidate, List<Integer> chosenVertiportID, List<TripItem> deserializedTripItems, String servedTripsIDFile) throws IOException {
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

        System.out.println("The number of served trips is "+uamAvailableTrips.size());

        List<Integer> uamAvailableTripsID = new ArrayList<>();
        for (TripItem tripItem : uamAvailableTrips) {
            uamAvailableTripsID.add(tripItem.tripID);
        }
        FileWriter fileWriter = new FileWriter(servedTripsIDFile);
        CSVWriter csvWriter = new CSVWriter(fileWriter);
        csvWriter.writeNext(new String[]{"tripID"});
        for (Integer tripID : uamAvailableTripsID) {
            csvWriter.writeNext(new String[]{tripID.toString()});
        }
        csvWriter.close();
        fileWriter.close();


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
