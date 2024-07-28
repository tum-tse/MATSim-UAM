package net.bhl.matsim.uam.optimization.utils;

public class ProcessorCheck {

        public static void main(String[] args) {
            int processors = Runtime.getRuntime().availableProcessors();
            System.out.println("Number of processors available to the Java Virtual Machine: " + processors);
        }

}
