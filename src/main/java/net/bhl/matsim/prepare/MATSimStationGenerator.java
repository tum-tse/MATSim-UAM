package net.bhl.matsim.prepare;

import org.matsim.api.core.v01.Coord;
import org.matsim.api.core.v01.population.Activity;
import org.matsim.api.core.v01.population.Person;
import org.matsim.api.core.v01.population.Plan;
import org.matsim.api.core.v01.population.Population;
import org.matsim.core.population.PopulationUtils;

import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class MATSimStationGenerator {

    // Inner class to represent a 2D point
    public static class Point {
        public double x, y;

        public Point(double x, double y) {
            this.x = x;
            this.y = y;
        }

        public double distanceTo(Point other) {
            return Math.sqrt(Math.pow(this.x - other.x, 2) + Math.pow(this.y - other.y, 2));
        }

        @Override
        public String toString() {
            return String.format("(%.2f, %.2f)", x, y);
        }
    }

    // Inner class to represent a cluster
    public static class Cluster {
        public Point centroid;
        public List<Point> points;

        public Cluster(Point centroid) {
            this.centroid = centroid;
            this.points = new ArrayList<>();
        }

        public void addPoint(Point point) {
            points.add(point);
        }

        public void updateCentroid() {
            if (points.isEmpty()) return;

            double sumX = 0, sumY = 0;
            for (Point point : points) {
                sumX += point.x;
                sumY += point.y;
            }
            centroid.x = sumX / points.size();
            centroid.y = sumY / points.size();
        }

        public void clearPoints() {
            points.clear();
        }
    }

    /**
     * Parses MATSim population file and extracts start activity locations
     * @param populationFilePath Path to the MATSim population XML file
     * @return List of start activity locations
     */
    public static List<Point> extractStartActivityLocations(String populationFilePath) {
        List<Point> startLocations = new ArrayList<>();

        try {
            // Read population using MATSim utilities
            Population population = PopulationUtils.readPopulation(populationFilePath);

            for (Person person : population.getPersons().values()) {
                Plan selectedPlan = person.getSelectedPlan();
                if (selectedPlan != null && !selectedPlan.getPlanElements().isEmpty()) {

                    // Get the first activity (start activity)
                    if (selectedPlan.getPlanElements().get(0) instanceof Activity) {
                        Activity startActivity = (Activity) selectedPlan.getPlanElements().get(0);
                        Coord coord = startActivity.getCoord();

                        if (coord != null) {
                            startLocations.add(new Point(coord.getX(), coord.getY()));
                        }
                    }
                }
            }

            System.out.println("Extracted " + startLocations.size() + " start activity locations");

        } catch (Exception e) {
            System.err.println("Error reading population file: " + e.getMessage());
            e.printStackTrace();
        }

        return startLocations;
    }

    /**
     * Performs K-means++ clustering on the given points
     * @param points List of points to cluster
     * @param k Number of clusters
     * @param maxIterations Maximum number of iterations
     * @return List of cluster centroids
     */
    public static List<Point> performKMeansPlusPlus(List<Point> points, int k, int maxIterations) {
        if (points.size() < k) {
            System.err.println("Number of points is less than number of clusters!");
            return new ArrayList<>();
        }

        Random random = new Random();
        List<Cluster> clusters = new ArrayList<>();

        // Step 1: Choose first centroid randomly
        Point firstCentroid = points.get(random.nextInt(points.size()));
        clusters.add(new Cluster(new Point(firstCentroid.x, firstCentroid.y)));

        // Step 2: Choose remaining centroids using K-means++ method
        for (int i = 1; i < k; i++) {
            double[] distances = new double[points.size()];
            double totalDistance = 0;

            // Calculate squared distances to nearest centroid for each point
            for (int j = 0; j < points.size(); j++) {
                Point point = points.get(j);
                double minDistance = Double.MAX_VALUE;

                for (Cluster cluster : clusters) {
                    double distance = point.distanceTo(cluster.centroid);
                    minDistance = Math.min(minDistance, distance);
                }

                distances[j] = minDistance * minDistance; // Squared distance
                totalDistance += distances[j];
            }

            // Choose next centroid with probability proportional to squared distance
            double randomValue = random.nextDouble() * totalDistance;
            double cumulativeDistance = 0;

            for (int j = 0; j < points.size(); j++) {
                cumulativeDistance += distances[j];
                if (cumulativeDistance >= randomValue) {
                    Point newCentroid = points.get(j);
                    clusters.add(new Cluster(new Point(newCentroid.x, newCentroid.y)));
                    break;
                }
            }
        }

        // Step 3: Perform K-means iterations
        for (int iteration = 0; iteration < maxIterations; iteration++) {
            // Clear all clusters
            for (Cluster cluster : clusters) {
                cluster.clearPoints();
            }

            // Assign each point to nearest cluster
            for (Point point : points) {
                double minDistance = Double.MAX_VALUE;
                Cluster nearestCluster = null;

                for (Cluster cluster : clusters) {
                    double distance = point.distanceTo(cluster.centroid);
                    if (distance < minDistance) {
                        minDistance = distance;
                        nearestCluster = cluster;
                    }
                }

                if (nearestCluster != null) {
                    nearestCluster.addPoint(point);
                }
            }

            // Update centroids
            boolean converged = true;
            for (Cluster cluster : clusters) {
                Point oldCentroid = new Point(cluster.centroid.x, cluster.centroid.y);
                cluster.updateCentroid();

                if (oldCentroid.distanceTo(cluster.centroid) > 1e-6) {
                    converged = false;
                }
            }

            if (converged) {
                System.out.println("K-means converged after " + (iteration + 1) + " iterations");
                break;
            }
        }

        // Extract centroids
        List<Point> centroids = new ArrayList<>();
        for (Cluster cluster : clusters) {
            centroids.add(cluster.centroid);
            System.out.println("Cluster centroid: " + cluster.centroid +
                    " with " + cluster.points.size() + " points");
        }

        return centroids;
    }

    /**
     * Generates stations.csv file from cluster centroids
     * @param centroids List of cluster centroids
     * @param outputFilePath Path for output CSV file
     */
    public static void generateStationsCSV(List<Point> centroids, String outputFilePath) {
        try (FileWriter writer = new FileWriter(outputFilePath)) {
            // Write header
            writer.append("station_id,station_name,x,y,z,vtol_z,ground_access_capacity,")
                    .append("ground_access_freespeed,flight_access_capacity,flight_access_freespeed,")
                    .append("preflighttime,postflighttime,defaultwaittime,numberOfChargers,chargingSpeed\n");

            // Write station data
            for (int i = 0; i < centroids.size(); i++) {
                Point centroid = centroids.get(i);
                writer.append(String.format("%d,station_%d,%.0f,%.0f,0,600,1000,50,500,100,300,180,600,99999,99999\n",
                        i + 1, i + 1, centroid.x, centroid.y));
            }

            System.out.println("Generated " + centroids.size() + " stations in " + outputFilePath);

        } catch (IOException e) {
            System.err.println("Error writing stations CSV: " + e.getMessage());
            e.printStackTrace();
        }
    }

    /**
     * Main method to execute the complete workflow
     */
    public static void main(String[] args) {
        // Configuration parameters
        String populationFilePath = "examples/munich-scenario/input/matsimPlans_5percent.xml.gz"; // Update with your file path
        String outputFilePath = "examples/munich-scenario/uam/stations.csv";
        int numberOfClusters = 10; // Adjust as needed
        int maxIterations = 100;

        // Parse command line arguments if provided
        if (args.length >= 1) {
            populationFilePath = args[0];
        }
        if (args.length >= 2) {
            numberOfClusters = Integer.parseInt(args[1]);
        }
        if (args.length >= 3) {
            outputFilePath = args[2];
        }

        System.out.println("Starting MATSim Station Generator...");
        System.out.println("Population file: " + populationFilePath);
        System.out.println("Number of clusters: " + numberOfClusters);
        System.out.println("Output file: " + outputFilePath);

        // Step 1: Extract start activity locations
        List<Point> startLocations = extractStartActivityLocations(populationFilePath);

        if (startLocations.isEmpty()) {
            System.err.println("No start locations found!");
            return;
        }

        // Step 2: Perform K-means++ clustering
        List<Point> centroids = performKMeansPlusPlus(startLocations, numberOfClusters, maxIterations);

        if (centroids.isEmpty()) {
            System.err.println("Clustering failed!");
            return;
        }

        // Step 3: Generate stations CSV
        generateStationsCSV(centroids, outputFilePath);

        System.out.println("Process completed successfully!");
    }

    /**
     * Utility method to run with custom parameters
     */
    public static void generateStations(String populationFilePath, int numberOfClusters, String outputFilePath) {
        List<Point> startLocations = extractStartActivityLocations(populationFilePath);
        if (!startLocations.isEmpty()) {
            List<Point> centroids = performKMeansPlusPlus(startLocations, numberOfClusters, 100);
            if (!centroids.isEmpty()) {
                generateStationsCSV(centroids, outputFilePath);
            }
        }
    }
}
