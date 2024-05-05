package net.bhl.matsim.uam.optimization;

import net.bhl.matsim.uam.optimization.utils.TripItemForOptimization;
import org.apache.log4j.Logger;

import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.List;

public class test {
    private static final Logger log = Logger.getLogger(test.class);
    public static void main(String[] args) {
        String serializedTripItemFile = "E:\\OneDrive - TUM\\MasterArbeit\\demand\\TLongTripsWithUtilityForRunFilippos1.dat";
        log.info("Loading the trips...");
        List<TripItemForOptimization> deserializedTripItemForOptimizations = deserializeTripItems(serializedTripItemFile);
        log.info("Finished loading the trips.");
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
