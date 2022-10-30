import time
from ast import literal_eval

import matplotlib.pyplot as plt
import openmm
import openmm.unit as unit
import pandas as pd
import scipy as sci
import simulate as sim
from openmm import XmlSerializer, app

import modules.small_molecule as smol
import modules.utils as ut
from modules.analyse import *


def minimize_montecarlo(args, real_path, logger, folder):
    if args.potential:
        minimize_potential_en(args, real_path, logger, folder)
    elif args.dispersion:
        minimize_dispersion(args, real_path, logger, folder)
    else:
        raise NotImplementedError()


def setup_montecarlo(args, real_path, logger, folder):

    logger.info("Performing Monte Carlo minimization.")
    print("")

    residues = pd.read_csv(f"{real_path}/data/residues.csv").set_index(
        "three", drop=False
    )

    try:
        proteins = pd.read_pickle(f"{real_path}/data/proteins.pkl")
    except ValueError:
        proteins = pd.read_csv(
            f"{real_path}/data/proteins.csv", converters={"fasta": literal_eval}
        )
        prot_new = proteins.rename(columns={"Unnamed: 0": "Name"})
        proteins = prot_new.set_index("Name")

    logger.info(f"Working with protein {args.name[0]} at {args.temp[0]} K.")

    t0 = time.time()

    # Adding additional arguments to the argument parser.
    vars(args)["proteins"] = proteins
    vars(args)["residues"] = residues

    residues = args.residues
    name = args.name[0]
    prot = proteins.loc[args.name[0]]
    temp = args.temp[0]
    sm_mol = args.small_molec
    sim_time = [args.time, args.nsteps]
    verbosity = args.verbosity
    plat_cpu = args.cpu
    plat_gpu = args.gpu
    check_collision = args.check_collision
    lambd_override = args.lambd
    sigma_override = args.sigma
    mass_override = args.mass

    # Number of MC runs where energy minimization is attempted.
    n_mc_runs = args.runs

    residues = residues.set_index("one")

    # folder = "test_mc_minimize_" + name + f"/{temp}/"

    # Generates the parameters for the LJ interaction º
    lj_eps, fasta, types, MWs = genParamsLJ(residues, name, prot)

    # Generates the parameters for the Debye-Huckel long range interaction,
    yukawa_eps, yukawa_kappa = genParamsDH(residues, name, prot, temp)

    N = len(fasta)

    # Preparing the base openmm.system that will contain the simulation object and parameters.
    system, L, Lz, marg, Nstep = sim.prepare_system(N)

    # Placing IDP chains in the simulation box.
    xy, top, pos, n_chains = sim.place_idp(L, Lz, N, marg, fasta, n_idp=100)

    # Storing the topology into a trajectory with one frame
    in_traj = md.Trajectory(
        np.array(pos).reshape(n_chains * N, 3),
        top,
        0,
        [L, L, Lz],
        [90, 90, 90],
    )

    # Attempting to create directories in which to save the topology
    ut.create_dirs(args)

    if args.resume:

        logger.info(f"Attempting to resume a simulation...")
        print("")

        # Checking if there is a checkpoint file
        chk_path, chk_file = ut.find_checkpoint()

        if chk_path != "":
            logger.info(f"Checkpoint file found in '{chk_path}.'")

        # print('check_point: ', check_point)
        check_point = chk_path + "/" + chk_file

        # check_point = chk_path + f"/{name}_{temp}_{sm_mol[0]}_restart.chk"
        # except TypeError:
        # check_point = chk_path + f"/{name}_{temp}_NODRG_restart.chk"
    else:
        logger.info(f"Starting new simulation...")
        print("")
        check_point = ""

    # Saving a .pdb file with the current configuration.
    logger.info(f"Storing files in {folder}")
    in_traj.save_pdb(f"{folder}/top.pdb")
    pdb = app.pdbfile.PDBFile(f"{folder}/top.pdb")

    for _ in range(n_chains):
        system.addParticle((residues.loc[prot.fasta[0]].MW + 2) * unit.amu)
        for a in prot.fasta[1:-1]:
            system.addParticle(residues.loc[a].MW * unit.amu)
        system.addParticle((residues.loc[prot.fasta[-1]].MW + 16) * unit.amu)

    logger.debug(in_traj.xyz.shape)
    logger.debug(top)
    n_parts_old = system.getNumParticles()

    if not os.path.isfile(check_point) and not args.resume:
        # This block of code will be executed if
        # there is not a checkpoint file

        if sm_mol:
            smol_name = sm_mol[0]
            smol_conc = float(sm_mol[1])
            comp_dist = float(sm_mol[2])

            drug_comp = smol_name.split("-")

            print("")
            logger.info("Adding small molecules to the system...")

            # Adding small particles to the system
            in_traj, top, system, n_drugs, drg_param = smol.add_drugs(
                system=system,
                in_traj=in_traj,
                in_top=top,
                conc=smol_conc,
                dist_threshold=2,
                drug_components=drug_comp,
                directory=folder,
                comp_dist=comp_dist,
                verbosity=verbosity,
                residues=residues,
                col_chk_flag=check_collision,
                lambd_override=lambd_override,
                mass_override=mass_override,
                sigma_override=sigma_override,
            )

            pdb = app.pdbfile.PDBFile(folder + "/sm_drg_traj.pdb")
            logger.debug(f"Number of particles: {system.getNumParticles()}")
            logger.debug(f"Number of drugs: {n_drugs}")

        else:
            logger.info("No small molecule given. Proceeding with only protein.")
            drg_param = "None"
            sigma = "None"

    else:
        # This code will get executed if there is a checkpoint file and
        # will load the previously used small molecule parameters and
        # configuration

        logger.info("\nReading small molecules from stored files...")
        top_ats = pd.read_csv(chk_path + "/sm_drg_ats.csv")

        # TODO: Add a way of getting the comp_dist used in a simulation
        # when resuming

        # This was before: n_drugs = (len(top_ats) - n_parts_old) // 2
        n_drugs = len(top_ats) - n_parts_old
        logger.info(f"number of drugs: {n_drugs}")

        top_bnd = np.load(chk_path + "/sm_drg_bnd.npy")
        top = md.Topology.from_dataframe(top_ats, top_bnd)

        in_traj = md.load(chk_path + "/final_system_state.pdb")

        pdb = app.pdbfile.PDBFile(chk_path + "/sm_drg_traj.pdb")
        top = pdb.getTopology()

        # logger.info(f"in_traj top: {in_traj.topology}")

        xml_path = [f for f in os.listdir(chk_path) if f.endswith(".xml")][0]

        with open(chk_path + "/" + xml_path, "r") as f:
            system_ser = f.read()
            system = XmlSerializer.deserialize(system_ser)

    # This block sets up particle interactions.
    logger.info("Setting bonded and non-bonded interactions...")

    # Adding a regular harmonic bond force
    hb = openmm.openmm.HarmonicBondForce()

    # This function defines the custom potentials used for the simulation
    ah, yu = sim.set_custom_potentials(yukawa_kappa, lj_eps)

    # This function  sets the parameters for the potentials of the AA of our main
    # chains.
    sim.set_AA_params(hb, ah, yu, n_chains, yukawa_eps, prot, N, residues)

    if not os.path.isfile(check_point) and not args.resume:
        # Adding the small drug particles to the CustomNonbondedForce used in the system.
        # n_drugs (number of small molecules) by 2 (bimolecular).
        if sm_mol:
            lambdas, sigma = sim.add_small_molec(
                sm_mol,
                residues,
                logger,
                lambd_override,
                sigma_override,
                n_drugs,
                yu,
                ah,
                hb,
                n_parts_old,
                comp_dist,
            )

        logger.debug(f"ah:, {ah.getNumParticles()}")
        logger.debug(f"yu:, {yu.getNumParticles()}")

        yu.setNonbondedMethod(openmm.openmm.CustomNonbondedForce.CutoffPeriodic)
        ah.setNonbondedMethod(openmm.openmm.CustomNonbondedForce.CutoffPeriodic)
        hb.setUsesPeriodicBoundaryConditions(True)

    # Generating system if a checkpoint file is not found.
    if not os.path.isfile(check_point) and not args.resume:
        system.addForce(hb)
        system.addForce(yu)
        system.addForce(ah)

        serialized_system = XmlSerializer.serialize(system)

        try:
            outfile = open(f"./{folder}/{name}_{temp}_{sm_mol[0]}_system.xml", "w")
        except TypeError:
            outfile = open(f"./{folder}/{name}_{temp}_NODRG_system.xml", "w")

        print("")
        logger.info("Generating '.xml' system file...")
        outfile.write(serialized_system)
        outfile.close()

    # Skippng system generation if a checkpoint file is found.
    else:
        logger.info("\nCheckpoint file found, skipping system generation.")
        logger.debug(f"System num parts: {system.getNumParticles()}\n")

        logger.debug(pdb.topology)

    integrator = openmm.openmm.LangevinIntegrator(
        temp * unit.kelvin, 0.01 / unit.picosecond, 0.005 * unit.picosecond
    )

    platform, platform_props = sim.select_platform(plat_cpu, plat_gpu, logger)

    sim_time, time_units, log_report_interval, dcd_save_interval = sim.select_timestep(
        sim_time
    )

    top = pdb.topology

    simulation = app.simulation.Simulation(
        top,
        system,
        integrator,
        platform,
        platformProperties=platform_props,
    )

    if os.path.isfile(check_point) and args.resume:
        logger.info("\nResuming simulation from checkpoint file...")
        simulation.loadCheckpoint(check_point)
        logger.info("Checkpoint loaded!")

        try:
            simulation.reporters.append(
                app.dcdreporter.DCDReporter(
                    f"{folder}/{name}_{temp}_{sm_mol[0]}_report_continue.dcd",
                    dcd_save_interval,
                    append=False,
                )
            )
        except TypeError:
            simulation.reporters.append(
                app.dcdreporter.DCDReporter(
                    f"{folder}/{name}_{temp}_NODRG_report_continue.dcd",
                    dcd_save_interval,
                )
            )
    else:
        # Saving simulation parameters into a 'parameter.dat' file to facilitate
        # data analysis and results parsing later
        ut.write_params(
            path=f"./{folder}/parameters.dat",
            args=args,
            sm_mol=sm_mol,
            drg_param=drg_param,
            sim_time=sim_time,
            time_units=time_units,
        )

        print("")
        logger.info("Starting simulation...")
        simulation.context.setPositions(pdb.positions)

        logger.info(
            "Initial potential energy:"
            f" {simulation.context.getState(getEnergy=True).getPotentialEnergy()}"
        )

        if not os.path.isfile(check_point) and not args.resume:
            # Initial minimization of the system
            simulation.minimizeEnergy()

            logger.info("Energy minimized.")
            logger.info(
                "Potential energy after minimization:"
                f" {simulation.context.getState(getEnergy=True).getPotentialEnergy()}"
            )

        try:
            simulation.reporters.append(
                app.dcdreporter.DCDReporter(
                    folder + f"/{name}_{temp}_{sm_mol[0]}_report.dcd",
                    sim_time,
                )
            )

        except TypeError:
            simulation.reporters.append(
                app.dcdreporter.DCDReporter(
                    folder + f"/{name}_{temp}_NODRG_report.dcd",
                    sim_time,
                )
            )

    # Checks if there is an existing .log file in order to select the .log file
    # writing mode
    if os.path.isfile(check_point):
        append_mode = True
        log_name = f"{folder}/{name}_{temp}_continue.log"
    else:
        log_name = f"{folder}/{name}_{temp}.log"
        append_mode = False

    # Generates log file with information
    simulation.reporters.append(
        app.statedatareporter.StateDataReporter(
            log_name,
            reportInterval=log_report_interval,
            potentialEnergy=True,
            temperature=True,
            step=True,
            speed=True,
            volume=True,
            elapsedTime=True,
            separator="\t",
            append=append_mode,
        )
    )

    # Save positions before starting the simulation.
    positions = simulation.context.getState(getPositions=True).getPositions()
    app.PDBFile.writeFile(
        simulation.topology,
        positions,
        open(f"{folder}/initial_minimization.pdb", "w"),
    )

    sim_time_ns = ut.timesteps_to_ns(sim_time)

    simulation.saveState(f"{folder}/res_mc_minimization_state.xml")

    return simulation, system, top, sim_time, time_units, sim_time_ns


def minimize_potential_en(args, real_path, logger, folder):
    """This function attempts to reduce the energy by performing short MD simulations
    and comparing the final energy obtained with the initial one. If the new energy
    is smaller, the final state of the simulation is kept as the system's current state,
    otherwise the trajectory is discarded and computed again.

    Args:
        args (argparse.Namespace): Namespace containing the user input arguments.
    """

    (
        simulation,
        system,
        top,
        sim_time,
        time_units,
        sim_time_ns,
    ) = setup_montecarlo(args, real_path, logger, folder)

    logger.info("Minimizing potential energy.")

    energy_list = []

    first_energ = simulation.context.getState(getEnergy=True).getPotentialEnergy()

    for run_ind in range(args.runs[0]):

        in_pot_energ = simulation.context.getState(getEnergy=True).getPotentialEnergy()

        energy_list.append(in_pot_energ._value)

        print("")
        logger.info(f".:: MC minimization step: {run_ind+1}/{args.runs[0]} ::.")

        logger.info(
            f"Running simulation for {sim_time} {time_units} ({sim_time_ns} ns)."
        )

        # Running the simulations for the given timesteps.
        simulation.step(sim_time)

        run_pot_energ = simulation.context.getState(getEnergy=True).getPotentialEnergy()

        if run_pot_energ < in_pot_energ:
            logger.info(f"Current energy: {in_pot_energ}")
            logger.info(f"Run energy: {run_pot_energ}")
            logger.info(f"Energy difference: {run_pot_energ-in_pot_energ}")
            logger.info("[!] Run minimized energy.")
            simulation.saveState(f"{folder}/res_mc_minimization_state.xml")
            in_pot_energ = run_pot_energ

        else:
            logger.info(f"Current energy: {in_pot_energ}")
            logger.info(f"Run energy: {run_pot_energ}")
            logger.info(f"Energy unchanged.")
            simulation.loadState(f"{folder}/res_mc_minimization_state.xml")

    print("")

    energ_array = np.array([range(args.runs[0]), energy_list]).T

    logger.info("Minimization done.")
    logger.info(f"Energy difference = {abs(in_pot_energ - first_energ)}")

    print("")

    res_filename = f"{folder}/res_mc_minimization_energies.log"

    logger.info(f"Energies stored in {res_filename}.")
    np.savetxt(
        res_filename,
        energ_array,
        header=f"{args.runs[0]} {sim_time}",
    )

    # Saving final system position.
    positions = simulation.context.getState(getPositions=True).getPositions()
    app.PDBFile.writeFile(
        simulation.topology,
        positions,
        open(f"{folder}/res_mc_finalsystem.pdb", "w"),
    )

    # Plotting energy variation
    plt.plot(energ_array[1:, 0], energ_array[1:, 1])
    plt.xlabel("Iteration")
    plt.ylabel("Potential Energy (kJ/mol)")
    plt.savefig(f"{folder}/res_mc_energ_plot.png", dpi=200)


def compute_dispersion(simulation, top):

    # Gathering all particle coordinates from the current simulation state into a
    # numpy array
    coord = (
        simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
        / unit.nanometer
    )

    # Length of chain 0. All chains should be equal except the last one,
    # which is the one containing the small molecules
    for chain in top.chains():
        idp_chain_count = len(list(top.chains())) - 1
        idp_chain_length = len(list(chain.atoms()))
        break

    # Splitting the coordinate array into sub arrays size equal to the atom count in
    # the chains. The last chain is omitted as it contains the small molecules.
    idp_chain_coords = np.vsplit(
        coord[: idp_chain_count * idp_chain_length, :], idp_chain_count
    )

    # Empty list for the distance averages.
    idp_chain_averages = []

    # Finding the average point of every chain and adding it to a list.
    for chain in idp_chain_coords:
        idp_chain_averages.append(np.average(chain, axis=0))

    # Finding the distance between all of the average points
    distance_allchains = sci.spatial.distance.pdist(np.array(idp_chain_averages))

    # Computing the average of all chain average distances.
    avg_distance_allchains = np.average(distance_allchains)

    return avg_distance_allchains


def minimize_dispersion(args, real_path, logger, folder):
    """This function attempts to reduce the particle dispersion by performing short MD simulations
    and comparing the final dispersion obtained with the initial one. If the new value
    is smaller, the final state of the simulation is kept as the system's current state,
    otherwise the trajectory is discarded and computed again.

    Args:
        args (argparse.Namespace): Namespace containing the user input arguments.
    """

    (
        simulation,
        system,
        topology,
        sim_time,
        time_units,
        sim_time_ns,
    ) = setup_montecarlo(args, real_path, logger, folder)

    logger.info("Minimizing particle dispersion.")

    energy_list = []

    # Storing initial dispersion.
    first_disp = compute_dispersion(simulation, topology)

    for run_ind in range(args.runs[0]):

        # Computing initial dispersion (in nm).
        in_disp = compute_dispersion(simulation, topology)

        # Adding initial dispersion to a list.
        energy_list.append(in_disp)

        print("")
        logger.info(f".:: MC minimization step: {run_ind+1}/{args.runs[0]} ::.")

        logger.info(
            f"Running simulation for {sim_time} {time_units} ({sim_time_ns} ns)."
        )

        # Running the simulations for the given timesteps.
        simulation.step(sim_time)

        # Computing current run's dispersion (in nm).
        run_disp = compute_dispersion(simulation, topology)

        # Checking dispersion distance difference.
        if run_disp < in_disp:
            logger.info(f"Current dispersion: {in_disp} nm")
            logger.info(f"Run dispersion: {run_disp} nm")
            logger.info(f"Dispersion difference: {run_disp-in_disp} nm")
            logger.info("[!] Run decreased dispersion.")
            simulation.saveState(f"{folder}/res_mc_minimization_state.xml")
            in_disp = run_disp

        else:
            logger.info(f"Current dispersion: {in_disp} nm")
            logger.info(f"Run dispersion: {run_disp} nm")
            logger.info(f"Dispersion unchanged.")
            simulation.loadState(f"{folder}/res_mc_minimization_state.xml")

    print("")

    energ_array = np.array([range(args.runs[0]), energy_list]).T

    logger.info("Minimization done.")
    logger.info(f"Total dispersion variation = {abs(in_disp - first_disp)} nm")

    print("")

    res_filename = f"{folder}/res_mc_minimization_dispersion.log"

    logger.info(f"Dispersion values stored in '{res_filename}'.")
    np.savetxt(
        res_filename,
        energ_array,
        header=f"{args.runs[0]} {sim_time}",
    )

    # Saving final system position.
    positions = simulation.context.getState(getPositions=True).getPositions()
    app.PDBFile.writeFile(
        simulation.topology,
        positions,
        open(f"{folder}/res_mc_finalsystem.pdb", "w"),
    )

    # Plotting energy variation
    plt.plot(energ_array[1:, 0], energ_array[1:, 1])
    plt.xlabel("Iteration")
    plt.ylabel("Particle Dispersion (nm)")
    plt.savefig(f"{folder}/res_mc_disp_plot.png", dpi=200)
