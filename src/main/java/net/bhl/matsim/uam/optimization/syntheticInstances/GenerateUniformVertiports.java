package net.bhl.matsim.uam.optimization.syntheticInstances;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Locale;

public class GenerateUniformVertiports {
    public static void main(String[] args) {
        Locale.setDefault(Locale.US);
        int[] sizes = {10000, 25000, 50000}; // Network Size
        int[] counts = {10,62,250}; // vertiports Number

        for (int size : sizes) {
            for (int count : counts) {
                generateAndSaveVertiports(size, count);
            }
        }
    }

    private static void generateAndSaveVertiports(int size, int count) {
        int numRows = (int) Math.sqrt(count); // 确定每行的vertiport数目
        double spacing = (double) size / (numRows - 1); // 计算间距
        String fileName = "uniform_vertiports_" + size + "_" + count + ".csv";

        try (FileWriter writer = new FileWriter(fileName)) {
            writer.append("ID,X,Y,Capacity\n"); // 写入表头
            int id = 1;
            for (int i = 0; i < numRows; i++) {
                for (int j = 0; j < numRows; j++) {
                    double x = i * spacing;
                    double y = j * spacing;
                    int capacity = 10;
                    writer.append(String.format("%d,%.2f,%.2f,%d\n", id++, x, y, capacity));
                }
            }
        } catch (IOException e) {
            System.out.println("Error writing to file: " + e.getMessage());
        }
    }
}

