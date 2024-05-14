package net.bhl.matsim.uam.optimization;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;
public class ModeDecider {
    public ModeDecider(double UAM_Utlility, double carUtility, double ptUtility, Random random) {
        this.UAM_Utlility = UAM_Utlility;
        this.carUtility = carUtility;
        this.ptUtility = ptUtility;
        this.random = random;
    }

    public double UAM_Utlility;
    public double carUtility;
    public double ptUtility;
    public Random random;



    // return the number of samples of each mode in a Integer array
    public Double [] sample(int samplesize) {
        Double [] probabilities = calculateModeProbability(UAM_Utlility, carUtility, ptUtility).toArray(new Double[3]);
        double [] samples = new double[probabilities.length];
        for(int i=0; i< samples.length; i++) {
            samples[i] = 0.0;
        }

        double sum = 0;
        for(int i=0; i< samplesize; i++) {
            double r = this.random.nextDouble();
            for (int j = 0; j < probabilities.length; j++) {
                sum += probabilities[j];
                if (r < sum) {
                    samples[j] += 1;
                    sum = 0;
                    break;
                }
            }
        }
        return new Double[]{ (samples[0] / samplesize), samples[1] / samplesize, samples[2] / samplesize};

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