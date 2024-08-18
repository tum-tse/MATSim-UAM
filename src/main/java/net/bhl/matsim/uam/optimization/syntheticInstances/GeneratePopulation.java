package net.bhl.matsim.uam.optimization.syntheticInstances;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Locale;
import java.util.Random;

public class GeneratePopulation {
    public static void main(String[] args) {
        int[] sizes = {10000, 25000, 50000}; // 网格大小：10km, 25km, 50km
        int[] densities = {1000,2500,5000}; // 人口密度：人/km2

        for (int size : sizes) {
            for (int density : densities) {
                generateAndSavePopulation(size, density);
            }
        }
    }

    private static void generateAndSavePopulation(int size, int density) {
        Locale.setDefault(Locale.US); // 确保使用点作为小数分隔符
        Random random = new Random();
        int totalPopulation = density * (size / 1000) * (size / 1000) /20; // 计算总人口
        String fileName = String.format("population_%dkm_%d_density.csv", size / 1000, density);

        try (FileWriter writer = new FileWriter(fileName)) {
            writer.append("ID,X,Y\n"); // 写入表头
            for (int i = 1; i <= totalPopulation; i++) {
                double x = random.nextDouble() * size;
                double y = random.nextDouble() * size;
                writer.append(String.format("%d,%.2f,%.2f\n", i, x, y));
            }
        } catch (IOException e) {
            System.out.println("Error writing to file: " + e.getMessage());
        }
    }
}

