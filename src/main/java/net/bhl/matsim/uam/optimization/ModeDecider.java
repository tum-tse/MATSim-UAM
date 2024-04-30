package net.bhl.matsim.uam.optimization;
import java.util.concurrent.ThreadLocalRandom;
public class ModeDecider {
    public Double [] probabilities;

    public ModeDecider(Double [] probabilities) {
        this.probabilities = probabilities;

    }


    // return the number of samples of each mode in a Integer array
    public Double [] sample(int samplesize) {
        Double [] samples = new Double[probabilities.length];
        double r = ThreadLocalRandom.current().nextDouble();
        double sum = 0;
        for(int i=0; i< samplesize; i++) {
            for (int j = 0; j < probabilities.length; j++) {
                sum += probabilities[j];
                if (r < sum) {
                    samples[j] += 1;
                    break;
                }
            }
        }
        return new Double[]{ (samples[0] / samplesize), samples[1] / samplesize, samples[2] / samplesize};

    }




}