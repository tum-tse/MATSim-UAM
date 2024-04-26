package net.bhl.matsim.uam.optimization;

import org.matsim.api.core.v01.Coord;
import org.matsim.contrib.util.CSVReaders;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class VertiportReader {
    public static List<Vertiport> getVertiports(String vertiportCandidateInput) throws IOException {
        List<Vertiport> vertiports = new ArrayList<>();
        List<String[]> rows = CSVReaders.readCSV(vertiportCandidateInput);
        for (String[] row : rows.subList(1, rows.size())) {
            int j = 0;
            Vertiport vertiport = new Vertiport();
            vertiport.ID = Integer.parseInt(row[j++]);
            vertiport.coord = new Coord(Double.parseDouble(row[j++]), Double.parseDouble(row[j++]));
            vertiport.constructionCost = Double.parseDouble(row[j]);
            vertiports.add(vertiport);

        }
        return vertiports;
    }
}
