package net.bhl.matsim.uam.optimization.syntheticInstances;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Locale;
import java.util.Random;

public class GeneratePopulationTrips {
    public static void main(String[] args) {
        int[] sizes = {10000, 25000, 50000}; // Network Size：10km, 25km, 50km
        int[] densities = {1000, 2500, 5000}; // Population Density: people/km2
        String[] distributionTypes = {"central", "dual", "uniform"}; // Destination Distribution Type
        double percent = 0.05; // Filter size

        for (int size : sizes) {
            for (int density : densities) {
                for (String type : distributionTypes) {
                    generateAndSavePopulationTrips(size, (int) (density* percent), type);
                }
            }
        }
    }

    private static void generateAndSavePopulationTrips(int size, int density, String distributionType) {
        Locale.setDefault(Locale.US); // Make sure to use "." as decimal separator
        Random random = new Random();
        int totalPopulation = (density * (size / 1000) * (size / 1000)); // 计算总人口
        String fileName = String.format("trips_%dkm_%d_%s.csv", size / 1000, density*20, distributionType);

        try (FileWriter writer = new FileWriter(fileName)) {
            writer.append("originX,originY,destinationX,destinationY,departureTime\n"); // 写入表头
            for (int i = 1; i <= totalPopulation; i++) {
                // get a random origin and destination within the network size
                double originX = random.nextDouble() * size;
                double originY = random.nextDouble() * size;
                double[] destination = getDestination(size, distributionType, random);
                int departureTime = 21600 + random.nextInt(10801); // 随机生成出发时间

                writer.append(String.format("%.2f,%.2f,%.2f,%.2f,%d\n", originX, originY, destination[0], destination[1], departureTime));
            }
        } catch (IOException e) {
            System.out.println("Error writing to file: " + e.getMessage());
        }
    }

    private static double[] getDestination(int size, String distributionType, Random random) {
        double centerX = size / 2.0;
        double centerY = size / 2.0;
        double spread = size / 4.0; // 用于控制分散区域的大小

        switch (distributionType) {
            case "central":
                // 终点集中于中心区域
                return new double[]{centerX + random.nextGaussian() * spread, centerY + random.nextGaussian() * spread};
            case "dual":
                // 两个中心，东西部
                boolean isEast = random.nextBoolean();
                double eastX = 3 * size / 4.0;
                double westX = size / 4.0;
                return new double[]{isEast ? eastX + random.nextGaussian() * spread : westX + random.nextGaussian() * spread, centerY + random.nextGaussian() * spread};
            case "uniform":
                // 终点均匀分散
                return new double[]{random.nextDouble() * size, random.nextDouble() * size};
            default:
                return new double[]{0, 0}; // 默认返回中心点
        }
    }
}

