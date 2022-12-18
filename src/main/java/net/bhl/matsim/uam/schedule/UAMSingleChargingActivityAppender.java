package net.bhl.matsim.uam.schedule;

import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ExecutionException;

import org.matsim.api.core.v01.Id;
import org.matsim.api.core.v01.network.Link;
import org.matsim.contrib.dvrp.path.VrpPathWithTravelData;
import org.matsim.contrib.dvrp.path.VrpPaths;
import org.matsim.contrib.dvrp.schedule.DriveTask;
import org.matsim.contrib.dvrp.schedule.Schedule;
import org.matsim.contrib.dvrp.schedule.Schedules;
import org.matsim.contrib.dvrp.schedule.StayTask;
import org.matsim.contrib.dvrp.schedule.Task;
import org.matsim.core.router.util.LeastCostPathCalculator;
import org.matsim.core.router.util.LeastCostPathCalculator.Path;
import org.matsim.core.router.util.TravelTime;

import com.google.inject.Inject;
import com.google.inject.name.Named;

import net.bhl.matsim.uam.infrastructure.UAMStation;
import net.bhl.matsim.uam.infrastructure.UAMStations;
import net.bhl.matsim.uam.infrastructure.UAMVehicle;
import net.bhl.matsim.uam.run.UAMConstants;

/**
 * This class adds tasks for each vehicle schedule based on the requests
 *
 * @author balacmi (Milos Balac), RRothfeld (Raoul Rothfeld)
 */
public class UAMSingleChargingActivityAppender {
	@Inject
	@Named(UAMConstants.uam)
	private LeastCostPathCalculator uamPathCalculator;

	@Inject
	@Named(UAMConstants.uam)
	private TravelTime travelTime;

	@Inject
	private UAMStations uamStations;

	private List<AppendTask> tasks = new LinkedList<>();

	/**
	 * @param request UAM request
	 * @param vehicle UAM Vehicle
	 * @param now     simulation step now
	 *                <p>
	 *                This method generates the paths for pickup and drop-off and
	 *                create a new AppendTask containing this information. The new
	 *                AppendTask is added to the AppendTask list.
	 */
	public void schedule(UAMVehicle vehicle, double now, Id<UAMStation> chargingStationId) {
		tasks.add(new AppendTask(vehicle, now, chargingStationId)); // adds the task to the tasks list
	}

	// There is a difference between an AppendTask and a Task that implements the
	// Task interface.

	// Uses the AppendTasks from the list containing the paths to generate the Tasks
	// (Tasks that implements the Task interface) and add them in order to the
	// vehicle schedule
	public void schedule(AppendTask task) throws ExecutionException, InterruptedException {
		UAMVehicle vehicle = task.vehicle;
		double now = task.time;
		Schedule schedule = vehicle.getSchedule();

		StayTask stayTask = (StayTask) Schedules.getLastTask(schedule);
		Link chargingLink = this.uamStations.getUAMStations().get(task.chargingStationId).getLocationLink();
		// either we are already there and then we just add charging task
		boolean requiresPickupFlight = !stayTask.getLink().getId().equals(chargingLink.getId());
		double startTime = stayTask.getStatus() == Task.TaskStatus.STARTED ? now : stayTask.getBeginTime();
		double scheduleEndTime = schedule.getEndTime();

		if (!requiresPickupFlight) {
			// we are already there
			UAMChargingTask chargingTask = new UAMChargingTask(UAMTaskType.CHARGING, startTime, startTime + 60.0,
					chargingLink, this.uamStations.getNearestUAMStation(chargingLink).getId());
			schedule.addTask(chargingTask);
		} else {
			Path pickupPath = uamPathCalculator.calcLeastCostPath(stayTask.getLink().getToNode(),
					chargingLink.getFromNode(), startTime, null, null);

			VrpPathWithTravelData chargePathWithTravelData = VrpPaths.createPath(stayTask.getLink(), chargingLink,
					startTime, pickupPath, travelTime);
			DriveTask chargingFlyTask = new DriveTask(UAMTaskType.FLY, chargePathWithTravelData);
			schedule.addTask(chargingFlyTask);
			now = chargingFlyTask.getEndTime();
			UAMChargingTask chargingTask = new UAMChargingTask(UAMTaskType.CHARGING, chargingFlyTask.getEndTime(),
					chargingFlyTask.getEndTime() + 60.0, chargingLink,
					this.uamStations.getNearestUAMStation(chargingLink).getId());
			schedule.addTask(chargingTask);
		}

		schedule.addTask(new StayTask(UAMTaskType.STAY, now + 60.0, scheduleEndTime, chargingLink));
	}

	public void update() {
		// TODO: This can be made more efficient if one knows which ones have
		// just been added and which ones are still
		// to be processed. Depends mainly on if "update" is called before new
		// tasks are submitted or after ...
		try {
			for (AppendTask task : tasks)
				schedule(task);
		} catch (ExecutionException | InterruptedException e) {
			throw new RuntimeException(e);
		}

		tasks.clear();
	}

	private class AppendTask {
		final public UAMVehicle vehicle;
		final public Id<UAMStation> chargingStationId;

		final public double time;

		public AppendTask(UAMVehicle vehicle, double time, Id<UAMStation> chargingStationId) {
			this.vehicle = vehicle;
			this.time = time;
			this.chargingStationId = chargingStationId;
		}
	}
}
