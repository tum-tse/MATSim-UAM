package net.bhl.matsim.uam.optimization;

import org.matsim.api.core.v01.Coord;
import org.matsim.contrib.util.CSVReaders;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class VertiportReader {
    public final double NEIGHBOR_DISTANCE = 1000.0;
    public List<Vertiport> getVertiports(String vertiportCandidateInput) throws IOException {
        List<Vertiport> vertiports = new ArrayList<>();
        List<String[]> rows = CSVReaders.readCSV(vertiportCandidateInput);
        for (String[] row : rows.subList(1, rows.size())) {
            int j = 0;
            Vertiport vertiport = new Vertiport();
            vertiport.ID = Integer.parseInt(row[j++]);
            vertiport.coord = new Coord(Double.parseDouble(row[j++]), Double.parseDouble(row[j++]));
            vertiport.capacity = Integer.parseInt(row[j]);
            vertiports.add(vertiport);
        }
        return vertiports;
    }

    public List<Vertiport> getVertiportsWithNeighbors(String vertiportCandidateInput) throws IOException {
        List<Vertiport> vertiports = new ArrayList<>();
        List<String[]> rows = CSVReaders.readCSV(vertiportCandidateInput);
        for (String[] row : rows.subList(1, rows.size())) {
            int j = 0;
            Vertiport vertiport = new Vertiport();
            vertiport.ID = Integer.parseInt(row[j++]);
            vertiport.coord = new Coord(Double.parseDouble(row[j++]), Double.parseDouble(row[j++]));
            vertiport.capacity = Integer.parseInt(row[j++]);
            vertiport.constructionCost = Double.parseDouble(row[j]);
            for (int i = 0; i < 24; i++) {
                vertiport.saturationRates.put(i, 0.0);
            }
            vertiport.neighbors = new ArrayList<>();
            vertiports.add(vertiport);
        }
        for (int i = 0; i < vertiports.size(); i++) {
            for (int j = 0; j < vertiports.size(); j++) {
                if (i != j) {
                    if (calculateEuciDistance(vertiports.get(i).coord, vertiports.get(j).coord) < NEIGHBOR_DISTANCE) {
                        vertiports.get(i).neighbors.add(vertiports.get(j));
                    }
                }
            }
        }
        return vertiports;
    }
    public double calculateEuciDistance(Coord coord1, Coord coord2) {
        double euciDistance = Math.sqrt(Math.pow(coord1.getX() - coord2.getX(), 2) + Math.pow(coord1.getY() - coord2.getY(), 2));
        return euciDistance;
    }
}
