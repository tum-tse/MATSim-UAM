package net.bhl.matsim.uam.optimization.syntheticInstances;

import org.matsim.api.core.v01.Coord;
import org.matsim.api.core.v01.Id;
import org.matsim.api.core.v01.network.Network;
import org.matsim.api.core.v01.network.Node;
import org.matsim.core.network.NetworkUtils;
import org.matsim.core.network.io.NetworkWriter;
import org.matsim.core.utils.geometry.CoordUtils;

public class NetworkGenerator {

    private Network network;

    public static void main(String[] args) {
        NetworkGenerator generator = new NetworkGenerator();

        // Generate raster networks
        generator.createGridNetwork(10000, 10000, 750, "raster10.xml");
        generator.createGridNetwork(25000, 25000, 750, "raster25.xml");
        generator.createGridNetwork(50000, 50000, 750, "raster50.xml");

        // Generate radial-ring networks
        generator.createRadialRingNetwork(10000, 750, 16, "radialRing10.xml");
        generator.createRadialRingNetwork(25000, 750, 16, "radialRing25.xml");
        generator.createRadialRingNetwork(50000, 750, 16, "radialRing50.xml");

        // Generate orthogonal-diagonal networks
        generator.createOrthogonalDiagonalNetwork(10000, 10000, 750, "orthoDiag10.xml");
        generator.createOrthogonalDiagonalNetwork(25000, 25000, 750, "orthoDiag25.xml");
        generator.createOrthogonalDiagonalNetwork(50000, 50000, 750, "orthoDiag50.xml");
    }

    private void createGridNetwork(double width, double height, double spacing, String filename) {
        this.network = NetworkUtils.createNetwork();
        int numX = (int) (width / spacing);
        int numY = (int) (height / spacing);
        // Create grid nodes
        for (int i = 0; i <= numX; i++) {
            for (int j = 0; j <= numY; j++) {
                Coord coord = new Coord(i * spacing, j * spacing);
                Node node = NetworkUtils.createAndAddNode(network, Id.createNodeId(i + "_" + j), coord);
                // Connect nodes bidirectionally
                if (i > 0) {
                    Node prevNodeX = network.getNodes().get(Id.createNodeId((i - 1) + "_" + j));
                    createBidirectionalLink(prevNodeX, node, spacing);
                }
                if (j > 0) {
                    Node prevNodeY = network.getNodes().get(Id.createNodeId(i + "_" + (j - 1)));
                    createBidirectionalLink(prevNodeY, node, spacing);
                }
            }
        }
        saveNetwork(filename);
    }

    private void createRadialRingNetwork(double diameter, double spacing, int numRadials, String filename) {
        this.network = NetworkUtils.createNetwork();
        double radius = diameter / 2;
        Coord center = new Coord(radius, radius);
        Node centralNode = NetworkUtils.createAndAddNode(network, Id.createNodeId("center"), center);

        // Create radial roads
        Node[][] ringNodes = new Node[(int)(radius / spacing) + 1][numRadials];
        for (int i = 0; i < numRadials; i++) {
            double angle = 2 * Math.PI * i / numRadials;
            for (int r = 0; r * spacing <= radius; r++) {
                double currentRadius = r * spacing;
                Coord coord = new Coord(currentRadius * Math.cos(angle) + radius, currentRadius * Math.sin(angle) + radius);
                Node node = NetworkUtils.createAndAddNode(network, Id.createNodeId("radial_" + i + "_ring_" + r), coord);
                ringNodes[r][i] = node;
                if (r == 0) {
                    createBidirectionalLink(centralNode, node, spacing);
                } else {
                    createBidirectionalLink(ringNodes[r-1][i], node, spacing);
                }
            }
        }

        // Connect radial roads with ring roads
        for (int r = 1; r * spacing <= radius; r++) {
            for (int i = 0; i < numRadials; i++) {
                int nextIndex = (i + 1) % numRadials;
                createBidirectionalLink(ringNodes[r][i], ringNodes[r][nextIndex], spacing);
            }
        }

        saveNetwork(filename);
    }

    private void createOrthogonalDiagonalNetwork(double width, double height, double spacing, String filename) {
        createGridNetwork(width, height, spacing, filename); // Create base bidirectional raster network
        int numX = (int) (width / spacing);
        int numY = (int) (height / spacing);
        // Add diagonal links in two directions
        for (int i = 0; i < numX; i++) {
            for (int j = 0; j < numY; j++) {
                Node node = network.getNodes().get(Id.createNodeId(i + "_" + j));
                if (node == null) continue;
                // Connect diagonally in both directions
                if (i < numX && j < numY) {
                    Node diagNode1 = network.getNodes().get(Id.createNodeId((i + 1) + "_" + (j + 1)));
                    Node diagNode2 = network.getNodes().get(Id.createNodeId((i + 1) + "_" + (j - 1)));
                    if (diagNode1 != null) {
                        createBidirectionalLink(node, diagNode1, Math.sqrt(2) * spacing);
                    }
                    if (diagNode2 != null && j > 0) {
                        createBidirectionalLink(node, diagNode2, Math.sqrt(2) * spacing);
                    }
                }
            }
        }
        saveNetwork(filename);
    }
    private void createBidirectionalLink(Node fromNode, Node toNode, double length) {
        NetworkUtils.createAndAddLink(network, Id.createLinkId("from_" + fromNode.getId() + "_to_" + toNode.getId()), fromNode, toNode, length, 50 / 3.6, 999999, 1);
        NetworkUtils.createAndAddLink(network, Id.createLinkId("from_" + toNode.getId() + "_to_" + fromNode.getId()), toNode, fromNode, length, 50 / 3.6, 999999, 1);
    }
    private void saveNetwork(String filename) {
        new NetworkWriter(this.network).write(filename);
    }
}
