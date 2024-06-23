package net.bhl.matsim.uam.optimization.syntheticInstances;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class RunPreCalculateInstance {
    public static void main(String[] args) {
        // Define multiple sets of arguments
        List<String[]> argsList = new ArrayList<>();
        argsList.add(new String[]{
                "D:\\OneDrive - TUM\\MasterArbeit\\OptimizingUAMVertiport\\syntheticInstancesInput\\Trips\\10_modified\\modified_trips_10km_orth_1000_central_car.csv",
                "D:\\OneDrive - TUM\\MasterArbeit\\OptimizingUAMVertiport\\syntheticInstancesInput\\NetworkConfig\\10\\config_10_OrthoDiag.xml",
                "D:\\OneDrive - TUM\\MasterArbeit\\OptimizingUAMVertiport\\syntheticInstancesInput\\VertiportsCandidates\\uniform_vertiports_10000_50.csv",
                "D:\\OneDrive - TUM\\MasterArbeit\\OptimizingUAMVertiport\\syntheticInstancesInput\\Trips\\10_modified\\modified_trips_10km_orth_1000_central_car.dat"
        });
        // More parameter sets can be added as needed

        // Execute the target program for each set of parameters
        for (String[] runArgs : argsList) {
            runPreCalculateAccessEgressCost(runArgs);
        }
    }

    private static void runPreCalculateAccessEgressCost(String[] args) {
        try {
            // Construct the command
            String[] command = createCommand(args);

            // Print the command being executed
            System.out.print("Executing: ");
            for (String cmd : command) {
                System.out.print(cmd + " ");
            }
            System.out.println();

            // Create the process
            ProcessBuilder processBuilder = new ProcessBuilder(command);
            processBuilder.inheritIO(); // Redirect child process output and error streams to the main process
            Process process = processBuilder.start();

            // Wait for the process to finish
            int exitCode = process.waitFor();
            System.out.println("Process exited with code: " + exitCode);
        } catch (IOException | InterruptedException e) {
            e.printStackTrace();
        }
    }

    private static String[] createCommand(String[] args) {
        String[] command = new String[args.length + 4];
        command[0] = "java";
        command[1] = "-cp";
        command[2] = "target/classes"; // Adjust this if your classes are in a different directory
        command[3] = "net.bhl.matsim.uam.optimization.PreCalculateAccessEgressCost";
        System.arraycopy(args, 0, command, 4, args.length);
        return command;
    }
}
