package net.bhl.matsim.uam.analysis.traveltimes;

import ch.ethz.matsim.av.plcpc.DefaultParallelLeastCostPathCalculator;
import ch.sbb.matsim.routing.pt.raptor.*;
import net.bhl.matsim.uam.config.UAMConfigGroup;
import net.bhl.matsim.uam.data.UAMRoute;
import net.bhl.matsim.uam.data.UAMStationConnectionGraph;
import net.bhl.matsim.uam.dispatcher.UAMManager;
import net.bhl.matsim.uam.infrastructure.UAMStations;
import net.bhl.matsim.uam.infrastructure.readers.UAMXMLReader;
import net.bhl.matsim.uam.modechoice.CustomCarDisutility;
import net.bhl.matsim.uam.qsim.UAMLinkSpeedCalculator;
import net.bhl.matsim.uam.router.UAMModes;
import net.bhl.matsim.uam.router.strategy.*;
import net.bhl.matsim.uam.router.strategy.UAMStrategy.UAMStrategyType;
import org.apache.log4j.Logger;
import org.matsim.api.core.v01.Coord;
import org.matsim.api.core.v01.Scenario;
import org.matsim.api.core.v01.TransportMode;
import org.matsim.api.core.v01.network.Link;
import org.matsim.api.core.v01.network.Network;
import org.matsim.contrib.dvrp.run.DvrpConfigGroup;
import org.matsim.core.config.Config;
import org.matsim.core.config.ConfigUtils;
import org.matsim.core.config.groups.PlanCalcScoreConfigGroup;
import org.matsim.core.controler.AbstractModule;
import org.matsim.core.controler.Injector;
import org.matsim.core.mobsim.qsim.qnetsimengine.linkspeedcalculator.DefaultLinkSpeedCalculator;
import org.matsim.core.network.NetworkUtils;
import org.matsim.core.network.algorithms.TransportModeNetworkFilter;
import org.matsim.core.router.AStarLandmarksFactory;
import org.matsim.core.router.DijkstraFactory;
import org.matsim.core.router.LinkWrapperFacility;
import org.matsim.core.router.util.*;
import org.matsim.core.scenario.ScenarioUtils;
import org.matsim.core.trafficmonitoring.TravelTimeCalculator;
import org.matsim.core.utils.misc.Time;
import org.matsim.facilities.Facility;

import java.io.*;
import java.util.*;

/**
 * This script generates csv file containing estimated travel times by UAM for
 * trips. The trips file must contain departure time and origin and destination
 * coordinates for the trips. Necessary inputs are in the following order:
 * -Network file; -Transit Schedule file; -Transit Vehicles file; -Trips file;
 * -strategy name(minTraveltime, minDistance, minAccessTravelTime,
 * minAccessDistance) -output file; ;
 *
 * @author Aitanm (Aitan Militao), RRothfeld (Raoul Rothfeld)
 */

public class RunCalculateUAMTravelTimes {
	private static final Logger log = Logger.getLogger(RunCalculateUAMTravelTimes.class);

	public static void main(String[] args) throws Exception {
		System.out.println(
				"ARGS: base-network.xml* uam.xml* transitScheduleFile.xml* transitVehiclesFile.xml* tripsCoordinateFile.csv* strategy-name* outputfile-name* total-process-times-mins search-radius-km access-modes");
		System.out.println("(* required)");

		// READ THE INPUTS
		// ARGS
		int j = 0;
		String networkInput = args[j++];
		String uamVehicles = args[j++];
		String transitScheduleInput = args[j++];
		String transitVehiclesInput = args[j++];
		String tripsInput = args[j++];
		String strategyName = args[j++];
		String outputPath = args[j++];

		double processTime = 0;
		if (args.length == j + 1)
			processTime = Double.parseDouble(args[j++]) * 60;

		double searchRadius = 50000;
		if (args.length == j + 1)
			searchRadius = Double.parseDouble(args[j++]) * 1000;

		String accessModes = "car,pt,bike,walk";
		if (args.length == j + 1)
			accessModes = args[j];

		// READ NETWORK
		Config config = ConfigUtils.createConfig(new UAMConfigGroup(), new DvrpConfigGroup());
		config.network().setInputFile(networkInput);
		// READ TRANSIT SCHEDULE
		config.transit().setTransitScheduleFile(transitScheduleInput);
		// READ TRANSIT VEHICLES FILE
		config.transit().setVehiclesFile(transitVehiclesInput);

		// CONFIGURING PT PARAMETERS
		config.transitRouter().setSearchRadius(2500);
		config.planCalcScore().setPerforming_utils_hr(0);
		config.planCalcScore().setMarginalUtlOfWaitingPt_utils_hr(-2.3);
		config.planCalcScore().setUtilityOfLineSwitch(-0.17);
		config.transitRouter().setExtensionRadius(500);

		PlanCalcScoreConfigGroup.ModeParams accessWalk = new PlanCalcScoreConfigGroup.ModeParams(TransportMode.access_walk);
		accessWalk.setMarginalUtilityOfTraveling(-4.0);
		config.planCalcScore().addModeParams(accessWalk);
		PlanCalcScoreConfigGroup.ModeParams egressWalk = new PlanCalcScoreConfigGroup.ModeParams(TransportMode.egress_walk);
		egressWalk.setMarginalUtilityOfTraveling(-4.0);
		config.planCalcScore().addModeParams(egressWalk);

		config.planCalcScore().getOrCreateModeParams(TransportMode.pt).setMarginalUtilityOfTraveling(-1.32);
		config.planCalcScore().getOrCreateModeParams(TransportMode.walk).setMarginalUtilityOfTraveling(-6.46);
		config.transit().setUseTransit(true);

		// setting walk and bike teleportation parameters equal to the test scenario
		config.plansCalcRoute().getModeRoutingParams().get(TransportMode.walk).setTeleportedModeSpeed(1.2);
		config.plansCalcRoute().getModeRoutingParams().get(TransportMode.walk).setBeelineDistanceFactor(1.3);
		config.plansCalcRoute().getModeRoutingParams().get(TransportMode.bike).setTeleportedModeSpeed(3.1);
		config.plansCalcRoute().getModeRoutingParams().get(TransportMode.bike).setBeelineDistanceFactor(1.4);

		// set available modes
		((UAMConfigGroup) config.getModules().get(UAMModes.UAM_MODE)).setAvailableAccessModes(accessModes);
		((UAMConfigGroup) config.getModules().get(UAMModes.UAM_MODE)).setSearchRadius("" + searchRadius);

		// Build scenario
		Scenario scenario = ScenarioUtils.createScenario(config);
		ScenarioUtils.loadScenario(scenario);
		Network network = scenario.getNetwork();

		// CREATE CAR/UAM NETWORK
		TransportModeNetworkFilter filter = new TransportModeNetworkFilter(network);
		Set<String> modes = new HashSet<>();
		modes.add(UAMModes.UAM_MODE);
		Network networkUAM = NetworkUtils.createNetwork();
		filter.filter(networkUAM, modes);
		Network networkCar = NetworkUtils.createNetwork();
		Set<String> modesCar = new HashSet<>();
		modesCar.add(TransportMode.car);
		filter.filter(networkCar, modesCar);

		// LEAST COST Parallel PATH CALCULATOR - uam
		TravelTimeCalculator tcc2 = TravelTimeCalculator.create(network, config.travelTimeCalculator());
		TravelTime travelTime = tcc2.getLinkTravelTimes();
		TravelDisutility travelDisutility = TravelDisutilityUtils
				.createFreespeedTravelTimeAndDisutility(config.planCalcScore());
		DefaultParallelLeastCostPathCalculator pathCalculator = DefaultParallelLeastCostPathCalculator.create(
				Runtime.getRuntime().availableProcessors(), new DijkstraFactory(), networkUAM, travelDisutility,
				travelTime);

		// SETUP UAM MANAGER AND STATIONCONENCTIONUTILITIES
		UAMXMLReader uamReader = new UAMXMLReader(networkUAM);
		uamReader.readFile(uamVehicles);
		UAMManager uamManager = new UAMManager(network);
		uamManager.setStations(new UAMStations(uamReader.getStations(), network));
		uamManager.setVehicles(uamReader.getVehicles());

		UAMStationConnectionGraph stationConnectionutilities = new UAMStationConnectionGraph(uamManager, null,
				pathCalculator);

		// PLCPCCAR
		com.google.inject.Injector injector = Injector.createInjector(config, new AbstractModule() {
			@Override
			public void install() {
				bind(LeastCostPathCalculatorFactory.class).to(AStarLandmarksFactory.class);
				// should not be necessary: created in the controler
				// install(new ScenarioByInstanceModule(scenario));
			}
		});

		DefaultLinkSpeedCalculator delegate = new DefaultLinkSpeedCalculator();
		UAMLinkSpeedCalculator speedCalculator = new UAMLinkSpeedCalculator(uamReader.getMapVehicleVerticalSpeeds(),
				uamReader.getMapVehicleHorizontalSpeeds(), delegate);
		CustomCarDisutility customCarDisutility = new CustomCarDisutility(speedCalculator, travelTime);
		LeastCostPathCalculatorFactory pathCalculatorFactory = injector
				.getInstance(LeastCostPathCalculatorFactory.class); // AStarLandmarksFactory
		LeastCostPathCalculator plcpccar = pathCalculatorFactory.createPathCalculator(networkCar, customCarDisutility,
				travelTime);

		// SET PUBLIC TRANSPORT ROUTER
		RaptorStaticConfig raptorStaticConfig = RaptorUtils.createStaticConfig(config);
		raptorStaticConfig.setBeelineWalkSpeed(0.9230769);
		raptorStaticConfig.setMinimalTransferTime(0);
		raptorStaticConfig.setBeelineWalkConnectionDistance(250);
		raptorStaticConfig.setMarginalUtilityOfTravelTimeAccessWalk_utl_s(-0.0017944);
		raptorStaticConfig.setMarginalUtilityOfTravelTimeEgressWalk_utl_s(-0.0017944);
		raptorStaticConfig.setMarginalUtilityOfTravelTimeWalk_utl_s(-0.0017944);
		SwissRailRaptorData data = SwissRailRaptorData.create(scenario.getTransitSchedule(), raptorStaticConfig,
				network);
		DefaultRaptorParametersForPerson parametersForPerson = new DefaultRaptorParametersForPerson(config);
		SwissRailRaptor transitRouter = new SwissRailRaptor(data, parametersForPerson,
				new LeastCostRaptorRouteSelector(), new DefaultRaptorIntermodalAccessEgress());

		// Create StrategyUtils
		// Using default parameters from the uamConfig group : walkDistance, etc.
		UAMStrategyUtils strategyUtils = new UAMStrategyUtils(uamManager.getStations(),
				(UAMConfigGroup) config.getModules().get(UAMModes.UAM_MODE), scenario, stationConnectionutilities,
				networkCar, transitRouter, pathCalculator, plcpccar, null);
		UAMStrategy strategy = null;
		switch (UAMStrategyType.valueOf(strategyName.toUpperCase())) {
			case MINTRAVELTIME:
				strategy = new UAMMinTravelTimeStrategy(strategyUtils);
				break;
			case MINACCESSTRAVELTIME:
				strategy = new UAMMinAccessTravelTimeStrategy(strategyUtils);
				break;
			case MINDISTANCE:
				strategy = new UAMMinDistanceStrategy(strategyUtils);
				break;
			case MINACCESSDISTANCE:
				strategy = new UAMMinAccessDistanceStrategy(strategyUtils);
				break;
			default:
				log.warn("Strategy not available, please provide an available strategy.");
				System.exit(-1);
		}

		// READ TRIPS INPUT
		BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(tripsInput)));
		String line;
		List<String> header = null;

		List<TripItem> trips = new ArrayList<>();
		while ((line = reader.readLine()) != null) {
			List<String> row = Arrays.asList(line.split(","));

			if (header == null) {
				header = row;
			} else {
				double originX = Double.parseDouble(row.get(header.indexOf("origin_x")));
				double originY = Double.parseDouble(row.get(header.indexOf("origin_y")));
				double destX = Double.parseDouble(row.get(header.indexOf("destination_x")));
				double destY = Double.parseDouble(row.get(header.indexOf("destination_y")));
				double departureTime = Time.parseTime(row.get(header.indexOf("trip_time")));

				Coord originCood = new Coord(originX, originY);
				Coord destinationCoord = new Coord(destX, destY);

				TripItem trip = new TripItem();
				trip.origin = originCood;
				trip.destination = destinationCoord;
				trip.departureTime = departureTime;
				trip.processTime = processTime;

				trips.add(trip);
			}
		}
		reader.close();

		// Calculate travel times
		log.info("Calculating travel times...");
		int counter = 1;
		for (TripItem trip : trips) {
			if (trips.size() < 100 || counter % (trips.size() / 100) == 0)
				log.info("Calculation completion: " + counter + "/" + trips.size() + " ("
						+ String.format("%.0f", (double) counter / trips.size() * 100) + "%).");
			Link from = NetworkUtils.getNearestLink(network, trip.origin);
			Link to = NetworkUtils.getNearestLink(network, trip.destination);
			Facility<?> fromFacility = new LinkWrapperFacility(from);
			Facility<?> toFacility = new LinkWrapperFacility(to);

			UAMRoute uamRoute = strategy.getRoute(null, fromFacility, toFacility, trip.departureTime);

			trip.accessTime = strategyUtils.getAccessTime(fromFacility, trip.departureTime,
					uamRoute.bestOriginStation, uamRoute.accessMode);
			trip.flightTime = strategyUtils.getFlightTime(uamRoute.bestOriginStation, uamRoute.bestDestinationStation);
			trip.egressTime = strategyUtils.getEgressTime(toFacility, trip.departureTime,
					uamRoute.bestDestinationStation, uamRoute.egressMode);

			trip.travelTime = trip.accessTime + trip.flightTime + trip.egressTime + trip.processTime;
			counter++;
		}

		pathCalculator.close();

		// Writes output file
		log.info("Writing travel times file...");
		write(outputPath, trips);
		log.info("...done.");

	}

	public static void write(String outputPath, List<TripItem> trips) throws IOException {
		BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outputPath)));

		writer.write(formatHeader() + "\n");
		for (TripItem trip : trips) {
			writer.write(String.join(",",
					new String[]{String.valueOf(trip.origin.getX()), String.valueOf(trip.origin.getY()),
							String.valueOf(trip.destination.getX()), String.valueOf(trip.destination.getY()),
							String.valueOf(trip.departureTime), String.valueOf(trip.travelTime),
							String.valueOf(trip.accessTime), String.valueOf(trip.flightTime),
							String.valueOf(trip.egressTime), String.valueOf(trip.processTime)})
					+ "\n");
		}

		writer.flush();
		writer.close();
	}

	private static String formatHeader() {
		return String.join(",", new String[]{"origin_x", "origin_y", "destination_x", "destination_y",
				"departure_time", "travel_time", "access_time", "flight_time", "egress_time", "process_time"});
	}

	static class TripItem {
		public Coord origin;
		public Coord destination;
		public double departureTime;
		public double travelTime;
		public double distance;
		public double accessTime;
		public double flightTime;
		public double egressTime;
		public double processTime;
	}

}
