package net.bhl.matsim.uam.optimization.utils;

public class ScenarioSpecific {
    public ScenarioSpecific(String scenarioName) {
        ScenarioSpecific.scenarioName = scenarioName;
    }

    public static String scenarioName;
    public int num_of_selected_vertiports;
    public double flight_speed;
    public double uam_process_time;
    public double uam_take_off_landing_time;
    public double uam_fix_cost;
    public double uam_km_cost;
    public double car_km_cost;
    public double uam_utility_cost_parameter;
    public double uam_utility_flight_time_parameter;
    public double uam_utility_waiting_time_parameter;
    public double uam_utility_generalized_cost_parameter;
    public boolean consider_return_trip;
    public double car_emission_factor;
    public double uam_emission_factor;
    public double pt_emission_factor;
    public double carbon_equivalent_cost;
    public int simulation_hours;
    public boolean consider_pt;
    public double search_radius;
    public double pt_cost;
    public double beta_savedCost; // the parameter of saved generalised cost in objective function
    public double beta_constructionCost; // the parameter of construction cost in objective function
    public double beta_savedEmission; // the parameter of saved emission in objective function
    public void buildScenario(){
        if (scenarioName.equals("Munich_A")) {
            this.uam_fix_cost = 6.1;
            this.uam_km_cost = 0.6;
            this.num_of_selected_vertiports = 74;
            this.flight_speed = 350 / 3.6;
            this.uam_process_time = 10 * 60;
            this.uam_take_off_landing_time = 1 * 60;
            this.car_km_cost = 0.42;
            this.uam_utility_cost_parameter = -2.48;
            this.uam_utility_flight_time_parameter = -4.28;
            this.uam_utility_waiting_time_parameter = -6.79;
            this.consider_return_trip = true;
            this.carbon_equivalent_cost = 2.48;
            this.car_emission_factor = 0.4376;
            this.pt_emission_factor = 0.1;
            this.uam_emission_factor = 0.1;
            this.simulation_hours = 36;
            this.consider_pt = false;
            this.search_radius = 5000;
            this.pt_cost = 0.51;
            this.beta_savedCost = 1;
            this.beta_constructionCost = 1;
            this.beta_savedEmission = 1;
        }
        if (scenarioName.equals("Synthetic_10")) {
            this.uam_fix_cost = 6.1;
            this.uam_km_cost = 0.6;
            this.num_of_selected_vertiports = 5;
            this.flight_speed = 350 / 3.6;
            this.uam_process_time = 10 * 60;
            this.uam_take_off_landing_time = 60;
            this.car_km_cost = 0.42;
            this.uam_utility_generalized_cost_parameter = -0.1;
            this.consider_return_trip = false;
            this.carbon_equivalent_cost = 0;
            this.car_emission_factor = 0;
            this.pt_emission_factor = 0;
            this.uam_emission_factor = 0;
            this.simulation_hours = 36;
            this.consider_pt = false;
            this.search_radius = 5000;
            this.pt_cost = 2;
        }
        if (scenarioName.equals("Synthetic_25")) {
            this.uam_fix_cost = 6.1;
            this.uam_km_cost = 0.6;
            this.num_of_selected_vertiports = 20;
            this.flight_speed = 350 / 3.6;
            this.uam_process_time = 10 * 60;
            this.uam_take_off_landing_time = 60;
            this.car_km_cost = 0.42;
            this.uam_utility_generalized_cost_parameter = -0.1;
            this.consider_return_trip = false;
            this.carbon_equivalent_cost = 0;
            this.car_emission_factor = 0;
            this.pt_emission_factor = 0;
            this.uam_emission_factor = 0;
            this.simulation_hours = 36;
            this.consider_pt = false;
            this.search_radius = 5000;
            this.pt_cost = 2;
        }
        if (scenarioName.equals("Synthetic_50")) {
            this.uam_fix_cost = 6.1;
            this.uam_km_cost = 0.6;
            this.num_of_selected_vertiports = 80;
            this.flight_speed = 350 / 3.6;
            this.uam_process_time = 10 * 60;
            this.uam_take_off_landing_time = 60;
            this.car_km_cost = 0.42;
            this.uam_utility_generalized_cost_parameter = -0.1;
            this.consider_return_trip = false;
            this.carbon_equivalent_cost = 0;
            this.car_emission_factor = 0;
            this.pt_emission_factor = 0;
            this.uam_emission_factor = 0;
            this.simulation_hours = 36;
            this.consider_pt = false;
            this.search_radius = 5000;
            this.pt_cost = 2;
        }
}

public double calculateUAMUtilityVAR(double flightTime, double waitingTime, double cost,double generalizedCost) {
        if (scenarioName.equals("Munich_A")) {
            return (uam_utility_cost_parameter * cost + uam_utility_flight_time_parameter * flightTime/60 + uam_utility_waiting_time_parameter * waitingTime/60)/100;
        }
        if (scenarioName.equals("Synthetic_10")) {
            return generalizedCost*uam_utility_generalized_cost_parameter;
        }
        if (scenarioName.equals("Synthetic_25")) {
            return generalizedCost*uam_utility_generalized_cost_parameter;
        }
        if (scenarioName.equals("Synthetic_50")) {
            return generalizedCost*uam_utility_generalized_cost_parameter;
        }
           return 0; // should not reach here
}
        }
