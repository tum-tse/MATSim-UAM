package net.bhl.matsim.uam.optimization.utils;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import java.io.File;

public class OptimizationConfiguration {
    public OptimizationConfiguration(String configuration_file) {
        this.configuration_file = configuration_file;
    }

    public String configuration_file;
    public double initial_temperature;
    public double final_temperature;
    public double annealing_rate;
    public int max_iteration;
    public int max_not_change_count;

    public void buildScenario() {
        // Read clusteringParameters
        try {
            File inputFile = new File(this.configuration_file);
            DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
            DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
            Document doc = dBuilder.parse(inputFile);
            doc.getDocumentElement().normalize();
            Element element = (Element) doc.getElementsByTagName("heuristicParameters").item(0);
            this.initial_temperature = Double.parseDouble(getTagValue("initial_temperature", element));
            this.final_temperature = Double.parseDouble(getTagValue("final_temperature", element));
            this.annealing_rate = Double.parseDouble(getTagValue("annealing_rate", element));
            this.max_iteration = Integer.parseInt(getTagValue("max_iteration", element));
            this.max_not_change_count = Integer.parseInt(getTagValue("max_not_change_count", element));

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
