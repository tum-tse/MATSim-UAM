package net.bhl.matsim.uam.optimization.utils;

public class ScenarioSpecific {
    public ScenarioSpecific(String scenarioName) {
        ScenarioSpecific.scenarioName = scenarioName;
    }

    public static String scenarioName;
    public static int num_of_selected_vertiports;
    public static double flight_speed;
    public static double uam_process_time;
    public static double uam_take_off_landing_time;
    public static double uam_fix_cost;
    public static double uam_km_cost;
    public static double car_km_cost;
    public static double uam_utility_cost_parameter;
    public static double uam_utility_flight_time_parameter;
    public static double uam_utility_waiting_time_parameter;
    public static boolean consider_return_trip;
    public static double car_emission_factor;
    public static double uam_emission_factor;
    public static double pt_emission_factor;
    public static double carbon_equivalent_cost;
    public static int simulation_hours;
    public void buildScenario(){
        if (scenarioName.equals("Munich_A")) {
            uam_fix_cost = 6.1;
            uam_km_cost = 0.6;
            num_of_selected_vertiports = 74;
            flight_speed = 350 / 3.6;
            uam_process_time = 10 * 60;
            uam_take_off_landing_time = 1 * 60;
            car_km_cost = 0.42;
            uam_utility_cost_parameter = -2.48;
            uam_utility_flight_time_parameter = -4.28;
            uam_utility_waiting_time_parameter = -6.79;
            consider_return_trip = true;
            carbon_equivalent_cost = 2.48;
            car_emission_factor = 0.42;
            pt_emission_factor = 0.1;
            uam_emission_factor = 0.1;
            simulation_hours = 36;

        }
}
        }
