package net.bhl.matsim.uam.optimization.utils;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.DocumentBuilder;
import org.w3c.dom.Document;
import org.w3c.dom.NodeList;
import org.w3c.dom.Node;
import org.w3c.dom.Element;
import java.io.File;
public class ClusteringConfiguration {
    public ClusteringConfiguration(String configuration_file) {
        this.configuration_file = configuration_file;
    }

    public String configuration_file;
    public int epsilon;
    public int min_pts;
    public int num_of_clusters;
    public int num_of_cluster_iterations;
    public double tolerance;
    public double neighbor_distance;
    public int num_cluster_in_existing_vertiports;
    public String clustering_output_file_candidates;
    public String clustering_output_file_existing;
    public void buildScenario() {
        // Read clusteringParameters
        try {
            File inputFile = new File(this.configuration_file);
            DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
            DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
            Document doc = dBuilder.parse(inputFile);
            doc.getDocumentElement().normalize();
        Element element = (Element) doc.getElementsByTagName("clusteringParameters").item(0);
        this.epsilon = Integer.parseInt(getTagValue("epsilon", element));
        this.min_pts = Integer.parseInt(getTagValue("min_pts", element));
        this.num_of_clusters = Integer.parseInt(getTagValue("num_of_clusters", element));
        this.num_of_cluster_iterations = Integer.parseInt(getTagValue("num_of_cluster_iterations", element));
        this.tolerance = Double.parseDouble(getTagValue("tolerance", element));
        this.neighbor_distance = Double.parseDouble(getTagValue("neighbor_distance", element));
        this.num_cluster_in_existing_vertiports = Integer.parseInt(getTagValue("num_of_clusters_existing", element));
        this.clustering_output_file_candidates = getTagValue("clustering_output_file_candidates", element);
        this.clustering_output_file_existing = getTagValue("clustering_output_file_existing", element);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }


    private static String getTagValue(String tag, Element element) {
        NodeList nodeList = element.getElementsByTagName(tag).item(0).getChildNodes();
        Node node = (Node) nodeList.item(0);
        return node.getNodeValue();
    }

}
