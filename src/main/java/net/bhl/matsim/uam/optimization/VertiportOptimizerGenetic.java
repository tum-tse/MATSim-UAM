package net.bhl.matsim.uam.optimization;

import com.opencsv.CSVWriter;
import net.bhl.matsim.uam.heuristic.Chromosome;
import net.bhl.matsim.uam.heuristic.Population;
import net.bhl.matsim.uam.optimization.utils.TripItemForOptimization;
import org.apache.log4j.Logger;
import org.matsim.api.core.v01.Coord;

import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.*;

public class VertiportOptimizerGenetic {
    public static final Logger log = Logger.getLogger(VertiportOptimizerGenetic.class);
    public static final double flightSpeed= 350/3.6; // m/s
    public static final double UAM_PROCESS_TIME= 10*60; // s
    public static final double takeOffLandingTime= 60; // s
    private static String vertiportCandidateFile;
    public static final double CURRENT_TOTAL_GENERALIZED_COST= 7825440.944198524;
    public static String fileName;
    public static String outputfile;

    public static double calculateEuciDistance(Coord coord1, Coord coord2) {
        double euciDistance = Math.sqrt(Math.pow(coord1.getX() - coord2.getX(), 2) + Math.pow(coord1.getY() - coord2.getY(), 2));
        return euciDistance;
    }
    public static void main(String[] args) throws IOException {
        // Provide the file via program arguments
        if (args.length > 0) {
            fileName = args[0];
            vertiportCandidateFile=args[1];
            outputfile=args[2];
        }
        // Get the object TripItemForOptimization from the serialized file
        log.info("Loading the vertiport candidates...");

        VertiportReader vertiportReader= new VertiportReader();
        List<Vertiport> vertiportsCandidates = vertiportReader.getVertiports(vertiportCandidateFile);
        log.info("Finished loading the vertiport candidates.");
        log.info("Loading the trips...");
        List<TripItemForOptimization> deserializedTripItemForOptimizations = deserializeTripItems(fileName);
        log.info("Finished loading the trips.");
        List<Vertiport> vertiports = new ArrayList<>();
       double vertiportCost = 0.0;
       double sumGeneralizedCost=0.0;

       // implement the population
        Population population = new Population(50,vertiportsCandidates, deserializedTripItemForOptimizations);
        // initialize the population
        population.initializePopulation();
        // calculate the fitness of all chromosomes in the population
        population.calculateFitnessOfAllChromosomes();
        // get the best chromosome in the population
        Chromosome bestChromosome = population.getbestChromosome(); // 初始化最优染色体
        double bestFitness = bestChromosome.fitness;

        // create a Map<Integer,Map<List<Integer>,Double>> to store the best fitness of each generation
        Map<Integer,List<Double>> bestFitnessList= new HashMap<>(); // key: generation (Integer), value: best fitness in this generation (Double), best fitness in all generations (Double)

        bestFitnessList.put(0,Arrays.asList(bestFitness,bestFitness));
        int nonchangeGeneration=0;
        for (int generation=1; generation <10000; generation++){
            // start update the new generation of population
            population.select();
            population.crossover();
            population.mutate();
            // end update the new generation of population

            // calculate the fitness of all chromosomes in the population
            population.calculateFitnessOfAllChromosomes();
            // get the best chromosome in the population
            Chromosome bestChromosomeInGeneration = population.getbestChromosome(); // 初始化最优染色体
            double bestFitnessInGeneration = bestChromosomeInGeneration.fitness;
            // put the best chromosome and best fitness in the Map<Integer,Map<List<Integer>,Double>> bestFitnessList
            if (bestFitnessInGeneration > bestFitness) {
                bestChromosome = bestChromosomeInGeneration; // 更新最优染色体
                bestFitness = bestFitnessInGeneration;
                nonchangeGeneration=0;
            }
            else
            {
                nonchangeGeneration++;
            }
            bestFitnessList.put(generation,Arrays.asList(bestFitnessInGeneration,bestFitness));
            // plot the change of best fitness across generations

            log.info("Generation " + generation + ": Best Fitness In This Generation = " + bestFitnessInGeneration);
            if (nonchangeGeneration>1000){
                break;
            }
        }
        // write the Map<Integer,Map<List<Integer>,Double>> bestFitnessList to a csv file
          writeCSV(bestFitnessList,outputfile);




        log.info("Best Chromosome: " + bestChromosome.genes);
        log.info("Best Fitness: " + bestFitness);
        log.info("Best Vertiport ID: " + bestChromosome.encodeChromosome());

       // check if the chosenVertiport is in the origin and destination neighbour vertiport of each trip in the deserializedTripItemForOptimizations


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
    public static void writeCSV(Map<Integer,List<Double>> bestFitnessList, String filePath) {
        try {
            CSVWriter writer = new CSVWriter(new FileWriter(filePath));
            // Integer is the generation, List<Double> contains the best fitness in this generation and the best fitness in all generations
            String[] header = {"Generation", "Best Fitness In This Generation", "Best Fitness In All Generations"};
            writer.writeNext(header);
            for (Map.Entry<Integer, List<Double>> entry : bestFitnessList.entrySet()) {
                Integer generation = entry.getKey();
                List<Double> bestFitnessInThisGenerationAndAllGenerations = entry.getValue();
                String[] data = {generation.toString(), bestFitnessInThisGenerationAndAllGenerations.get(0).toString(), bestFitnessInThisGenerationAndAllGenerations.get(1).toString()};
                writer.writeNext(data);
            }
            writer.close();


        } catch (IOException e) {
            e.printStackTrace();
        }
    }



}