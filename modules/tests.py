import time
from ast import literal_eval

import openmm
import openmm.unit as unit
import pandas as pd
import simulate as sim
from openmm import XmlSerializer, app

import modules.small_molecule as smol
import modules.utils as ut
from modules.analyse import *


def minimize_montecarlo(args, real_path, logger, folder):
    """This function attempts to reduce the energy by performing short MD simulations
    and comparing the final energy obtained with the initial one. If the new energy
    is smaller, the final state of the simulation is kept as the system's current state,
    otherwise the trajectory is discarded and computed again.

    Args:
        args (argparse.Namespace): Namespace containing the user input arguments.
    """

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

    # This block sets up particle interactions.
    logger.info("Setting bonded and non-bonded interactions...")

    # Adding a regular harmonic bond force
    hb = openmm.openmm.HarmonicBondForce()

    # This function defines the custom potentials used for the simulation
    ah, yu = sim.set_custom_potentials(yukawa_kappa, lj_eps)

    # This function  sets the parameters for the potentials of the AA of our main
    # chains.
    sim.set_AA_params(hb, ah, yu, n_chains, yukawa_eps, prot, N, residues)

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

    # Saving simulation parameters into a 'parameter.dat' file to facilitate
    # data analysis and results parsing later
    ut.write_params(
        path=f"./{folder}/parameters.dat",
        name=name,
        temp=temp,
        sm_mol=sm_mol,
        drg_param=drg_param,
        sim_time=sim_time,
        time_units=time_units,
        sigma=sigma,
        mass=mass_override,
        extension=args.extend_thermostat,
    )

    print("")
    logger.info("Starting simulation...")
    simulation.context.setPositions(pdb.positions)

    logger.info(
        "Initial potential energy:"
        f" {simulation.context.getState(getEnergy=True).getPotentialEnergy()}"
    )

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
                dcd_save_interval,
            )
        )

    except TypeError:
        simulation.reporters.append(
            app.dcdreporter.DCDReporter(
                folder + f"/{name}_{temp}_NODRG_report.dcd",
                dcd_save_interval,
            )
        )

    # Generates log file with information
    simulation.reporters.append(
        app.statedatareporter.StateDataReporter(
            f"{folder}/{name}_{temp}.log",
            reportInterval=log_report_interval,
            potentialEnergy=True,
            temperature=True,
            step=True,
            speed=True,
            volume=True,
            elapsedTime=True,
            separator="\t",
            append=False,
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

    simulation.saveCheckpoint(f"{folder}/res_minimization_checkpoint.chk")

    energy_list = []

    for run_ind in range(n_mc_runs[0]):

        print("")
        logger.info(f"MC minimization step: {run_ind}")

        in_pot_energ = simulation.context.getState(getEnergy=True).getPotentialEnergy()
        

        logger.info(
            f"Running simulation for {sim_time} {time_units} ({sim_time_ns} ns)."
        )

        # Running the simulations for the given timesteps.
        simulation.step(sim_time)

        run_pot_energ = simulation.context.getState(getEnergy=True).getPotentialEnergy()


        energy_list.append(run_pot_energ._value)

        if run_pot_energ < in_pot_energ:
            logger.info("Run minimized energy.")
            logger.info(f'Energy difference: {in_pot_energ-run_pot_energ}')
            simulation.saveCheckpoint(f"{folder}/res_minimization_checkpoint.chk")
        else:
            logger.info("Energy unchanged.")
            simulation.loadCheckpoint(f"{folder}/res_minimization_checkpoint.chk")

    print("")
    logger.info(f"MC minimization done.")
    logger.info(f"Energies stored in '{folder}/res_minimization_energies.log'.")
    np.savetxt(f"{folder}/res_minimization_energies.log", np.array([range(n_mc_runs[0]), energy_list]).T)
    
    # Saving final system position.
    positions = simulation.context.getState(getPositions=True).getPositions()
    app.PDBFile.writeFile(
        simulation.topology,
        positions,
        open(f"{folder}/res_finalsystem.pdb", "w"),
    )