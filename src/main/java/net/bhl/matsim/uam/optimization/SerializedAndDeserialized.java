package net.bhl.matsim.uam.optimization;

import net.bhl.matsim.uam.optimization.utils.TripItemForOptimization;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

public class SerializedAndDeserialized implements Serializable{
    public static String fileName = "C:\\Users\\24524\\OneDrive - TUM\\MasterArbeit\\demand\\motorized_trips_for_run.dat";
    public static String outputTripFile = "C:\\Users\\24524\\OneDrive - TUM\\MasterArbeit\\demand\\tripItemsForRunNew.dat";
    public static void main(String[] args) {
       // deserialize the trips
        List<TripItemForOptimization> deserializedTripItemForOptimizations = deserializeTripItems(fileName);
        


        // Assuming 'trips' is a List<Trip> containing all the trip data
        try (FileOutputStream fileOut = new FileOutputStream(outputTripFile);
             ObjectOutputStream out = new ObjectOutputStream(fileOut)) {
            out.writeObject(deserializedTripItemForOptimizations);
        } catch (IOException e) {
            e.printStackTrace();
        }
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
}
