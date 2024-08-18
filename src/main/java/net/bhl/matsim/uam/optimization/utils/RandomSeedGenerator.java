package net.bhl.matsim.uam.optimization.utils;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;
public class RandomSeedGenerator {




            public static void main(String[] args) {
                // Create a primary random generator for seed generation
                Random seedGenerator = new Random(); // You can optionally set a seed here
                // Array to hold the seeds
                long[] seeds = new long[10];

                // Generate 10 unique seeds
                for (int i = 0; i < seeds.length; i++) {
                    seeds[i] = seedGenerator.nextLong();
                    System.out.print(seeds[i]+",");
                }
            }
    }


