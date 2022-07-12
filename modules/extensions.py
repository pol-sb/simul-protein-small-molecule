import openmm
import openmm.unit as unit
from openmm import XmlSerializer, app

from modules.analyse import *


def extend_thermostat(
    args,
    logger,
    top,
    system,
    platform,
    platform_props,
    check_point,
    folder,
    name,
    temp,
    log_report_interval,
    chk_reporter_flag,
):
    # Getting the temperature and n_steps from the argument parser.
    ext_T = args.extend_thermostat[0]
    ext_steps = args.extend_thermostat[1]

    # Printing information
    logger.info(
        f"\nMain simulation done.\nApplying thermostat at {ext_T}K and extending"
        f" simulation for {int(ext_steps)} steps."
    )

    # Updating the integrator with the new thermostat temperature.
    integrator = openmm.openmm.LangevinIntegrator(
        ext_T * unit.kelvin, 0.01 / unit.picosecond, 0.005 * unit.picosecond
    )

    simulation_ext = app.simulation.Simulation(
        top,
        system,
        integrator,
        platform,
        platformProperties=platform_props,
    )

    simulation_ext.loadCheckpoint(check_point)

    simulation_ext.reporters.append(
        app.statedatareporter.StateDataReporter(
            f"{folder}{name}_{temp}.log",
            reportInterval=log_report_interval,
            potentialEnergy=True,
            temperature=True,
            step=True,
            speed=True,
            volume=True,
            elapsedTime=True,
            separator="\t",
            append=True,
        )
    )

    if not chk_reporter_flag:
        # Adding a checkpoint reporter manually, as the step method does not
        # have an option to automatically add the reporter.
        simulation_ext.reporters.append(
            app.CheckpointReporter(
                file=check_point,
                reportInterval=ext_steps * 0.05,
            )
        )

        # Adding a flag when the checkpoint file reporter is added to the simulation
        # to avoid adding it more than once.
        chk_reporter_flag = True

    # Running the simulations for the given timesteps.
    simulation_ext.step(int(ext_steps))

    # Printing information
    logger.info(f"\nSimulation extension done.")

    # Save post extension final system position.
    positions = simulation_ext.context.getState(getPositions=True).getPositions()
    app.PDBFile.writeFile(
        simulation_ext.topology,
        positions,
        open(f"{folder}final_system_state.pdb", "w"),
    )

    # Saves checkpoint file after extension
    simulation_ext.saveCheckpoint(check_point)
