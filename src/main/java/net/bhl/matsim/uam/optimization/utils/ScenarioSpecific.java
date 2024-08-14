package net.bhl.matsim.uam.optimization.utils;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.DocumentBuilder;
import org.w3c.dom.Document;
import org.w3c.dom.NodeList;
import org.w3c.dom.Node;
import org.w3c.dom.Element;
import java.io.File;
public class ScenarioSpecific {
    public ScenarioSpecific(String scenarioConfigurationFile) {
        this.scenarioConfigurationFile = scenarioConfigurationFile;
    }
    public  String scenarioConfigurationFile;
    public  boolean incrementaSiting;
    public  String existingVertiportFile;
    public long random_seed_PS; // random seed for population sampling
    public double population_sampling_factor;
    public long random_seed_RT; // random seed for remaining tasks
    public int num_of_selected_vertiports;
    public double flight_speed;
    public double uam_process_time;
    public double uam_take_off_landing_time;
    public double uam_cruise_altitude;
    public double uam_fix_cost;
    public double uam_km_cost;
    public double car_km_cost;
    public double uam_utility_cost_parameter;
    public double uam_utility_flight_time_parameter;
    public double uam_utility_waiting_time_parameter;
    public boolean consider_return_trip;
    public double car_emission_factor;
    public double uam_emission_factor_horizontal;
    public double uam_emission_factor_vertical;
    public double pt_emission_factor;
    public double carbon_equivalent_cost;
    public int simulation_hours;
    public boolean consider_pt; // whether to consider public transport as access and egress mode
    public double search_radius;
    public double pt_cost;
    public double beta_savedCost; // the parameter of saved generalised cost in objective function
    public double beta_constructionCost; // the parameter of construction cost in objective function
    public double beta_savedEmission; // the parameter of saved emission in objective function
    public double scale_factor_utility;
    public double neighbouring_distance;
    public int MonteCarlosampleSize;
    public double car_utility_mean;
    public double car_utility_sigma;
    public double pt_utility_mean;
    public double pt_utility_sigma;
    public String outputVertiportFile;
    public String outputTripBasedIndicatorFile;
    public String outputVertiportBasedIndicatorFile;
    public double wait_area_demand_factor;
    public int threshold_GRD_Unit_Selection;


    public void buildScenario() {
        try {
            File inputFile = new File(scenarioConfigurationFile);
            DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
            DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
            Document doc = dBuilder.parse(inputFile);
            doc.getDocumentElement().normalize();

            // Read scenarioSpecific parameters
            Element scenarioSpecific = (Element) doc.getElementsByTagName("scenarioSpecific").item(0);

            this.uam_fix_cost = Double.parseDouble(getTagValue("uam_fix_cost", scenarioSpecific));
            this.uam_km_cost = Double.parseDouble(getTagValue("uam_km_cost", scenarioSpecific));
            this.num_of_selected_vertiports = Integer.parseInt(getTagValue("num_of_selected_vertiports", scenarioSpecific));
            this.flight_speed = Double.parseDouble(getTagValue("flight_speed", scenarioSpecific));
            this.uam_process_time = Double.parseDouble(getTagValue("uam_process_time", scenarioSpecific));
            this.uam_cruise_altitude = Double.parseDouble(getTagValue("uam_cruise_altitude", scenarioSpecific));
            this.uam_take_off_landing_time = Double.parseDouble(getTagValue("uam_cruise_altitude", scenarioSpecific))/Double.parseDouble(getTagValue("uam_vertical_speed", scenarioSpecific))*2;
            this.car_km_cost = Double.parseDouble(getTagValue("car_km_cost", scenarioSpecific));
            this.uam_utility_cost_parameter = Double.parseDouble(getTagValue("uam_utility_cost_parameter", scenarioSpecific));
            this.uam_utility_flight_time_parameter = Double.parseDouble(getTagValue("uam_utility_flight_time_parameter", scenarioSpecific));
            this.uam_utility_waiting_time_parameter = Double.parseDouble(getTagValue("uam_utility_waiting_time_parameter", scenarioSpecific));
            this.consider_return_trip = Boolean.parseBoolean(getTagValue("consider_return_trip", scenarioSpecific));
            this.carbon_equivalent_cost = Double.parseDouble(getTagValue("carbon_equivalent_cost", scenarioSpecific));
            this.car_emission_factor = Double.parseDouble(getTagValue("car_emission_factor", scenarioSpecific));
            this.pt_emission_factor = Double.parseDouble(getTagValue("pt_emission_factor", scenarioSpecific));
            this.uam_emission_factor_horizontal = Double.parseDouble(getTagValue("uam_emission_factor", scenarioSpecific));
            this.uam_emission_factor_vertical = Double.parseDouble(getTagValue("uam_emission_factor_vertical", scenarioSpecific));
            this.simulation_hours = Integer.parseInt(getTagValue("simulation_hours", scenarioSpecific));
            this.consider_pt = Boolean.parseBoolean(getTagValue("consider_pt", scenarioSpecific));
            this.search_radius = Integer.parseInt(getTagValue("search_radius", scenarioSpecific));
            this.pt_cost = Double.parseDouble(getTagValue("pt_cost", scenarioSpecific));
            this.beta_savedCost = Double.parseDouble(getTagValue("beta_savedCost", scenarioSpecific));
            this.beta_constructionCost = Double.parseDouble(getTagValue("beta_constructionCost", scenarioSpecific));
            this.beta_savedEmission = Double.parseDouble(getTagValue("beta_savedEmission", scenarioSpecific));
            this.scale_factor_utility = Double.parseDouble(getTagValue("scale_factor_utility", scenarioSpecific));
            this.neighbouring_distance = Double.parseDouble(getTagValue("neighbouring_distance", scenarioSpecific));
            this.random_seed_RT = Long.parseLong(getTagValue("random_seed_RT", scenarioSpecific));
            this.random_seed_PS = Long.parseLong(getTagValue("random_seed_PS", scenarioSpecific));
            this.population_sampling_factor = Double.parseDouble(getTagValue("population_sampling_factor", scenarioSpecific));
            this.MonteCarlosampleSize = Integer.parseInt(getTagValue("Monte_carlo_sample_size", scenarioSpecific));
            this.car_utility_mean = Double.parseDouble(getTagValue("car_utility_mean", scenarioSpecific));
            this.car_utility_sigma = Double.parseDouble(getTagValue("car_utility_sigma", scenarioSpecific));
            this.pt_utility_mean = Double.parseDouble(getTagValue("pt_utility_mean", scenarioSpecific));
            this.pt_utility_sigma = Double.parseDouble(getTagValue("pt_utility_sigma", scenarioSpecific));
            this.existingVertiportFile = getTagValue("existing_vertiports_file", scenarioSpecific);
            this.incrementaSiting = Boolean.parseBoolean(getTagValue("incremental_siting", scenarioSpecific));
            this.outputTripBasedIndicatorFile = getTagValue("output_trip_indicators_file", scenarioSpecific);
            this.outputVertiportFile = getTagValue("output_vertiports_file", scenarioSpecific);
            this.outputVertiportBasedIndicatorFile = getTagValue("output_vertiport_indicators_file", scenarioSpecific);
            this.wait_area_demand_factor = Double.parseDouble(getTagValue("wait_area_demand_factor", scenarioSpecific));
            this.threshold_GRD_Unit_Selection = Integer.parseInt(getTagValue("threshold_GRD_Unit_Selection", scenarioSpecific));
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static String getTagValue(String tag, Element element) {
        NodeList nodeList = element.getElementsByTagName(tag).item(0).getChildNodes();
        Node node = (Node) nodeList.item(0);
        return node.getNodeValue();
    }



public double calculateUAMUtilityVAR(double flightTime, double waitingTime, double cost) {

    return (uam_utility_cost_parameter * cost + uam_utility_flight_time_parameter * flightTime / 60 + uam_utility_waiting_time_parameter * waitingTime / 60) / scale_factor_utility;

}
        }
