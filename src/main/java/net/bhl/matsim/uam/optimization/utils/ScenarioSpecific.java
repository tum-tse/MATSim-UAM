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
    public boolean consider_return_trip;
    public double car_emission_factor;
    public double uam_emission_factor;
    public double pt_emission_factor;
    public double carbon_equivalent_cost;
    public int simulation_hours;
    public boolean consider_pt;
    public double search_radius;
    public double pt_cost;
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
            this.car_emission_factor = 0.42;
            this.pt_emission_factor = 0.1;
            this.uam_emission_factor = 0.1;
            this.simulation_hours = 36;
            this.consider_pt = false;
            this.search_radius = 5000;
            this.pt_cost = 0.81;
        }
}
        }
