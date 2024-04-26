package net.bhl.matsim.uam.optimization;

public class ModeDecider {
    public ModeDecider(double carProbability, double ptProbability, double uamProbability) {
        this.carProbability = carProbability;
        this.ptProbability = ptProbability;
        this.uamProbability = uamProbability;
    }

    public double carProbability;
public double ptProbability;
public double uamProbability;


    public String decideMode() {

Double random=Math.random();
// Decide the mode by Mote carlo sampling based on the probability
if (random<this.carProbability){
 return "car";
}
else if (random<this.carProbability+this.ptProbability){
return "pt";
}
else{
    return "uam";
}
    }


}