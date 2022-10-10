import argparse
import os
import time
from ast import literal_eval

import openmm
import openmm.unit as unit
from openmm import XmlSerializer, app

import modules.small_molecule as smol
import modules.extensions as ext
import modules.utils as ut
from modules.analyse import *


def simulate(args):
    # arg_dict = args.__dict__

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
    # print('sim_time: ', sim_time)
    # print('sm_mol: ', sm_mol)
    # print('temp: ', temp)
    # print('prot: ', prot)
    # print('name: ', name)
    # print('verbosity: ', verbosity)
    # print('plat_cpu: ', plat_cpu)
    # print('plat_gpu: ', plat_gpu)
    # print('check_collision: ', check_collision)
    # print('lambd_override: ', lambd_override)
    # print('sigma_override: ', sigma_override)
    # print('mass_override: ', mass_override)
    # quit()

    folder = name + f"/{temp}/"

    try:
        check_point = folder + f"{name}_{temp}_{sm_mol[0]}_restart.chk"
    except TypeError:
        check_point = folder + f"{name}_{temp}_NODRG_restart.chk"

    residues = residues.set_index("one")

    # Generates the parameters for the LJ interaction using values in the
    # .csv file by employing the genParamsLJ function.
    lj_eps, fasta, types, MWs = genParamsLJ(residues, name, prot)

    # Generates the parameters for the Debye-Huckel long range interaction,
    # normally called Yukawa,  which includes an exponential term for ion
    # computed using values in the .csv file by employing the genParamsLJ
    # function.
    # This returns the epsilon (of 1/4*pi*eps) and the kappa for the
    # exponential term.
    yukawa_eps, yukawa_kappa = genParamsDH(residues, name, prot, temp)

    N = len(fasta)

    # setting slab parameters
    L = 15.0
    marg = 2
    if N > 400:
        L = 25.0
        Lz = 300.0
        marg = 8
        Nsteps = int(2e7)
    elif N > 200:
        L = 17.0
        Lz = 300.0
        marg = 4
        Nsteps = int(6e7)
    else:
        Lz = 10 * L
        Nsteps = int(6e7)

    system = openmm.System()

    # Set box vectors
    # The following lines define the size of each vector of the slab
    # that makes up our system.
    a = unit.Quantity(np.zeros([3]), unit.nanometers)
    a[0] = L * unit.nanometers
    b = unit.Quantity(np.zeros([3]), unit.nanometers)
    b[1] = L * unit.nanometers
    c = unit.Quantity(np.zeros([3]), unit.nanometers)
    c[2] = Lz * unit.nanometers

    # The following line sets the vectors for our system.
    system.setDefaultPeriodicBoxVectors(a, b, c)

    # Initial configuration
    xy = np.empty(0)

    # Generates a 1x2 vector with two x and y coordinates.
    xy = np.append(xy, np.random.rand(2) * (L - marg) - (L - marg) / 2).reshape((-1, 2))

    # Generates 100 new coordinates for the xy array, which must fulfill
    # certain conditions:
    for x, y in np.random.rand(1000, 2) * (L - marg) - (L - marg) / 2:
        x1 = x - L if x > 0 else x + L
        y1 = y - L if y > 0 else y + L
        if np.all(np.linalg.norm(xy - [x, y], axis=1) > 0.7):
            if np.all(np.linalg.norm(xy - [x1, y], axis=1) > 0.7):
                if np.all(np.linalg.norm(xy - [x, y1], axis=1) > 0.7):
                    xy = np.append(xy, [x, y]).reshape((-1, 2))
        if xy.shape[0] == 100:
            break

    n_chains = xy.shape[0]

    top = md.Topology()
    pos = []
    for x, y in xy:
        chain = top.add_chain()

        # Adds the positions at x, y and z.
        # Uses 0.38 because of the 3.8A separation between C in the backbone
        pos.append([[x, y, Lz / 2 + (i - N / 2.0) * 0.38] for i in range(N)])

        # Adds each residue to the chain, giving it an atom (C alpha atom)
        for resname in fasta:
            residue = top.add_residue(resname, chain)
            top.add_atom(resname, element=md.element.carbon, residue=residue)
        # Adds a bond between each atom and the next one
        for i in range(chain.n_atoms - 1):
            top.add_bond(chain.atom(i), chain.atom(i + 1))

    # Storing the topology into a trajectory with one frame
    in_traj = md.Trajectory(
        np.array(pos).reshape(n_chains * N, 3),
        top,
        0,
        [L, L, Lz],
        [90, 90, 90],
    )

    # Attempting to create directories in which to save the topology
    try:
        os.mkdir(f"./{name}")
    except FileExistsError:
        pass

    try:
        os.mkdir(f"./{name}/{int(temp)}/")
    except FileExistsError:
        pass

    # Saving a .pdb file with the current configuration. This file has 100
    # protein strands spanning a long stretch of the z axis and completely
    # straight.
    logger.info(f"Storing files in {os.getcwd()}/{name}/{int(temp)}/")
    in_traj.save_pdb(f"./{name}/{int(temp)}/top.pdb")
    pdb = app.pdbfile.PDBFile(f"./{name}/{int(temp)}/top.pdb")

    # Adding finally the particles to the system, with their charge and a
    # term for compensating for the terminal residues (2 for a terminal NH2
    # and 16 for the deprotonated OH of the ending carbonyl group)

    # TODO: Read mass from residues.csv in sm_molule.py

    for _ in range(n_chains):
        system.addParticle((residues.loc[prot.fasta[0]].MW + 2) * unit.amu)
        for a in prot.fasta[1:-1]:
            system.addParticle(residues.loc[a].MW * unit.amu)
        system.addParticle((residues.loc[prot.fasta[-1]].MW + 16) * unit.amu)

    logger.debug(in_traj.xyz.shape)
    logger.debug(top)
    n_parts_old = system.getNumParticles()

    # TODO: This should be defined in a separate input file and read as a dict
    # or be given as a parameter in the argument parser.

    if not os.path.isfile(check_point):
        if sm_mol:
            # TODO: Use these variables in the sm_mol function to get
            # parameters. Maybe save the small molecule params in a csv file?
            smol_name = sm_mol[0]
            smol_conc = float(sm_mol[1])
            comp_dist = float(sm_mol[2])

            drug_comp = smol_name.split("-")

            logger.info("\nAdding small molecules to the system...")
            # My function to add small particles to the system
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

            pdb = app.pdbfile.PDBFile(folder + "sm_drg_traj.pdb")
            logger.debug(f"Number of particles: {system.getNumParticles()}")
            logger.debug(f"Number of drugs: {n_drugs}")

        else:
            logger.info("No small molecule given. Proceeding with only protein.")
            drg_param = "None"
            sigma = "None"

    else:
        logger.info("\nReading small molecules from stored files...")
        top_ats = pd.read_csv(folder + "sm_drg_ats.csv")

        # TODO: Add a way of getting the comp_dist used in a simulation
        # when resuming

        # This was before: n_drugs = (len(top_ats) - n_parts_old) // 2
        n_drugs = len(top_ats) - n_parts_old
        logger.info(f"number of drugs: {n_drugs}")

        top_bnd = np.load(folder + "sm_drg_bnd.npy")
        top = md.Topology.from_dataframe(top_ats, top_bnd)

        in_traj = md.load(folder + "final_system_state.pdb")

        pdb = app.pdbfile.PDBFile(folder + "sm_drg_traj.pdb")
        top = pdb.getTopology()

        # logger.info(f"in_traj top: {in_traj.topology}")

        with open(
            f"./{name}/{int(temp)}/{name}_{temp}_{sm_mol[0]}_system.xml", "r"
        ) as f:
            system_ser = f.read()
            system = XmlSerializer.deserialize(system_ser)

        # logger.info(system.getForces())

    logger.info("\nSetting bonded and non-bonded interactions...")
    # Adding a regular harmonic bond force
    hb = openmm.openmm.HarmonicBondForce()

    # Adding our custom energy expression for a LJ type interaction to use in
    # the FF. It is non bonded as LJ is not bonded. Ashbaugh-Hatch (ah)
    # functional form.

    # Original lambda interaction definition: l=0.5*(l1+l2)
    energy_expression = "select(step(r-2^(1/6)*s),4*eps*l*((s/r)^12-(s/r)^6),4*eps*((s/r)^12-(s/r)^6)+eps*(1-l))"
    ah = openmm.openmm.CustomNonbondedForce(
        energy_expression + "; s=0.5*(s1+s2); l=0.5*(l1+l2)"
    )

    # Adding our custom energy expression for a Yukawa type interaction to use
    # in the FF. It is also non bonded.
    yu = openmm.openmm.CustomNonbondedForce(
        "q*(exp(-kappa*r)/r - exp(-kappa*4)/4); q=q1*q2"
    )

    # Adding our parameters to the FF
    yu.addGlobalParameter("kappa", yukawa_kappa / unit.nanometer)
    yu.addPerParticleParameter("q")

    ah.addGlobalParameter("eps", lj_eps * unit.kilojoules_per_mole)
    ah.addPerParticleParameter("s")
    ah.addPerParticleParameter("l")

    # This loop sets the parameters for the potentials of the AA of our main
    # chains.
    for j in range(n_chains):
        begin = j * N
        end = j * N + N

        for a, e in zip(prot.fasta, yukawa_eps):
            # print("e: ", e)
            yu.addParticle([e * unit.nanometer * unit.kilojoules_per_mole])
            ah.addParticle(
                [
                    residues.loc[a].sigmas * unit.nanometer,
                    residues.loc[a].lambdas * unit.dimensionless,
                ]
            )

        for i in range(begin, end - 1):
            hb.addBond(
                i,
                i + 1,
                0.38 * unit.nanometer,
                8033.28 * unit.kilojoules_per_mole / (unit.nanometer**2),
            )

            # Important, this makes pair not affect eachother with non bonded
            # potentials
            yu.addExclusion(i, i + 1)
            ah.addExclusion(i, i + 1)

    # Adding the small drug particles to the CustomNonbondedForce used in the system.
    # n_drugs (number of small molecules) by 2 (bimolecular).
    if sm_mol:
        for drg_ind, type_drg in enumerate(sm_mol[0].split("-")):

            # Finding the row of the AA corresponding to our type of molecule.
            res_row = residues.loc[residues["three"] == f"{type_drg}"]

            # Checking if non default lambda was given.
            if lambd_override or lambd_override == 0:
                lambdas = lambd_override[drg_ind]
                logger.info(
                    f"For {type_drg}-{drg_ind} using non-default lambda: {lambdas}"
                )
            else:
                lambdas = res_row["lambdas"][0]
                logger.info(f"For {type_drg}-{drg_ind} using default lambda: {lambdas}")

            # Checking if non default sigma was given.
            if sigma_override or sigma_override == 0:
                sigma = sigma_override[drg_ind]
                logger.info(
                    f"For {type_drg}-{drg_ind} using non-default sigma: {sigma}"
                )
            else:
                sigma = res_row["sigmas"][0]
                logger.info(f"For {type_drg}-{drg_ind} using default sigma: {sigma}")

            for i in range(n_drugs):

                # Yukawa Epsilon of the small molecules
                yu.addParticle(
                    [res_row["q"][0] * unit.nanometer * unit.kilojoules_per_mole]
                )

                ah.addParticle(
                    [
                        sigma * unit.nanometer,
                        lambdas * unit.dimensionless,
                    ]
                )

        # Adding bonds between the small molecules.
        sm_chain_length = len(sm_mol[0].split("-"))

        if sm_chain_length == 1:
            for i in range(
                n_parts_old,
                n_parts_old + (n_drugs * len(sm_mol[0].split("-"))) - 1,
            ):

                # Important, this makes pairs do not affect eachother with non bonded
                # potentials
                yu.addExclusion(i, i + 1)
                ah.addExclusion(i, i + 1)

        else:
            # Counter to keep track of the total number of added particles.
            sm_cnt = n_parts_old

            # print('sm_cnt: ', sm_cnt)
            # print('n_drugs: ', n_drugs)
            # quit()

            for i in range(n_drugs):
                # print("\nChain", i)
                for j in range(sm_chain_length - 1):
                    # print(f"bond: {sm_cnt}-{sm_cnt+1}")
                    hb.addBond(
                        sm_cnt,
                        sm_cnt + 1,
                        comp_dist * unit.nanometer,
                        8033.28 * unit.kilojoules_per_mole / (unit.nanometer**2),
                    )

                    # Important, this makes pairs do not affect eachother with non bonded
                    # potentials.
                    yu.addExclusion(sm_cnt, sm_cnt + 1)
                    ah.addExclusion(sm_cnt, sm_cnt + 1)

                    # Increasing the atom count by one for every atom in the chain.
                    sm_cnt += 1

                # Adding one to the atom count to avoid connecting two chains by their first and and first atoms.
                sm_cnt += 1

    logger.debug(f"ah:, {ah.getNumParticles()}")
    logger.debug(f"yu:, {yu.getNumParticles()}")

    yu.setNonbondedMethod(openmm.openmm.CustomNonbondedForce.CutoffPeriodic)
    ah.setNonbondedMethod(openmm.openmm.CustomNonbondedForce.CutoffPeriodic)
    hb.setUsesPeriodicBoundaryConditions(True)

    # Setting the cutoff. Using the value of 4nm, but it can be modified.
    # Another line above.
    yu.setCutoffDistance(4 * unit.nanometer)
    ah.setCutoffDistance(4 * unit.nanometer)

    if not os.path.isfile(check_point):
        system.addForce(hb)
        system.addForce(yu)
        system.addForce(ah)

        serialized_system = XmlSerializer.serialize(system)

        try:
            outfile = open(f"./{folder}{name}_{temp}_{sm_mol[0]}_system.xml", "w")
        except TypeError:
            outfile = open(f"./{folder}{name}_{temp}_NODRG_system.xml", "w")

        logger.info("Generating '.xml' system file...")
        outfile.write(serialized_system)
        outfile.close()

    else:
        logger.info("\nCheckpoint file found, skipping system generation.")
        logger.debug(f"System num parts: {system.getNumParticles()}\n")

        logger.debug(pdb.topology)

    # Defining a Langevin Integrator with the following parameters
    # temperature – the temperature of the heat bath (in Kelvin)
    # frictionCoeff – the frict. coeff. whcih couples the system to the heat bath
    # stepSize – the step size with which to integrate the system (in picoseconds)

    integrator = openmm.openmm.LangevinIntegrator(
        temp * unit.kelvin, 0.01 / unit.picosecond, 0.005 * unit.picosecond
    )

    if plat_cpu and not plat_gpu:
        # Using the CPU for the calculations
        platform = openmm.Platform.getPlatformByName("CPU")
        platform_props = None
        logger.info("\nUsing CPU for the calculations.")
        logger.info(f"This platform's speed score: {platform.getSpeed()}")
    else:
        # Uses CUDA as the platform for the GPU calculations.
        gpu_names = ut.get_gpus()
        logger.info(
            f"\nUsing CUDA for the calculations. GPU(s) available:\n{gpu_names}"
        )
        if not plat_gpu:
            logger.info(f"\nNo GPU index given, using index 0 as fallback.")
            platform = openmm.Platform.getPlatformByName("CUDA")
            platform_props = dict(CudaPrecision="mixed", DeviceIndex="0")
        else:
            plat_gpu_str = str(plat_gpu).replace("[", "").replace("]", "")
            logger.info(f"\nUsing GPU(s): {plat_gpu_str}.")
            platform = openmm.Platform.getPlatformByName("CUDA")
            platform_props = dict(CudaPrecision="mixed", DeviceIndex=plat_gpu_str)
            logger.info(f"This platform's speed score: {platform.getSpeed()}")

    # This conditional block checks the simulation time format given,
    # allowing to simulate a certain number of seconds or a number of
    # timesteps.

    log_report_interval = 1000

    if sim_time[0] and not sim_time[1]:
        sim_time = sim_time[0]
        time_units = "seconds"

    elif sim_time[1] and not sim_time[0]:
        sim_time = sim_time[1]
        time_units = "iterations"

    if time_units == "iterations":
        if sim_time < 100000:
            dcd_save_interval = 1000
        else:
            dcd_save_interval = 50000

    elif time_units == "seconds":
        if sim_time < 1000:
            log_report_interval = 10
            dcd_save_interval = 100
        else:
            dcd_save_interval = 50000

    if not os.path.isfile(check_point):
        top = pdb.topology

    simulation = app.simulation.Simulation(
        top,
        system,
        integrator,
        platform,
        platformProperties=platform_props,
    )

    if os.path.isfile(check_point):
        logger.info("\nResuming simulation from checkpoint file...")
        simulation.loadCheckpoint(check_point)
        logger.info("Checkpoint loaded!")

    if os.path.isfile(check_point):
        try:
            simulation.reporters.append(
                app.dcdreporter.DCDReporter(
                    name + f"/{temp}/{name}_{temp}_{sm_mol[0]}_report.dcd",
                    dcd_save_interval,
                    append=True,
                )
            )
        except TypeError:
            simulation.reporters.append(
                app.dcdreporter.DCDReporter(
                    name + f"/{temp}/{name}_{temp}_NODRG_report.dcd",
                    dcd_save_interval,
                )
            )
    else:
        # Saving simulation parameters into a 'parameter.dat' file to facilitate
        # data analysis and results parsing later
        ut.write_params(
            path=f"./{folder}parameters.dat",
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

        logger.info("\nStarting simulation...\n")
        simulation.context.setPositions(pdb.positions)

        logger.info(
            "Initial potential energy:"
            f" {simulation.context.getState(getEnergy=True).getPotentialEnergy()}"
        )

        # simulation.reporters.append(PDBReporter(‘output.pdb’, 1))

        simulation.minimizeEnergy()

        logger.info("\nEnergy minimized.")
        logger.info(
            "Potential energy after minimization:"
            f" {simulation.context.getState(getEnergy=True).getPotentialEnergy()}"
        )

        try:
            simulation.reporters.append(
                app.dcdreporter.DCDReporter(
                    name + f"/{temp}/{name}_{temp}_{sm_mol[0]}_report.dcd",
                    dcd_save_interval,
                )
            )

        except TypeError:
            simulation.reporters.append(
                app.dcdreporter.DCDReporter(
                    name + f"/{temp}/{name}_{temp}_NODRG_report.dcd",
                    dcd_save_interval,
                )
            )

    # Checks if there is an existing .log file in order to select the .log file
    # writing mode
    if os.path.isfile(check_point):
        append_mode = True
    else:
        append_mode = False

    # Generates log file with information
    simulation.reporters.append(
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
            append=append_mode,
        )
    )

    # Save positions before starting the simulation.
    positions = simulation.context.getState(getPositions=True).getPositions()
    app.PDBFile.writeFile(
        simulation.topology,
        positions,
        open(f"{folder}minimized_system.pdb", "w"),
    )

    logger.info(f"\nRunning simulation for {sim_time} {time_units}.")
    chk_reporter_flag = None

    # The checkpoint save time scales with the simulation length. The interval
    # can be adjusted.
    # Initial runtime was 20h.
    if time_units == "seconds":
        simulation.runForClockTime(
            sim_time * unit.second,
            checkpointFile=check_point,
            checkpointInterval=sim_time * 0.06 * unit.second,
        )

    # Alternative method of running the simulation
    elif time_units == "iterations":
        # Adding a checkpoint reporter manually, as the step method does not
        # have an option to automatically add the reporter.
        simulation.reporters.append(
            app.CheckpointReporter(
                file=check_point,
                reportInterval=sim_time * 0.05,
            )
        )
        # Adding a flag when the checkpoint file reporter is added to the simulation
        # to avoid adding it more than once.
        chk_reporter_flag = True

        # Running the simulations for the given timesteps.
        simulation.step(sim_time)

    # Saves final checkpoint file
    simulation.saveCheckpoint(check_point)

    # Save final system position.
    positions = simulation.context.getState(getPositions=True).getPositions()
    app.PDBFile.writeFile(
        simulation.topology,
        positions,
        open(f"{folder}final_system_state.pdb", "w"),
    )

    # Extending the simulation with a thermostat.
    if args.extend_thermostat:
        ext.extend_thermostat(
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
        )


if __name__ == "__main__":
    args = ut.arg_parse()
    # print('args: ', args)

    # Custom logger for easier debugging, using the python logging module.
    logger, verbosity = ut.custom_logger(args)

    # Getting the main script path
    real_path = os.path.split(os.path.realpath(__file__))[0]

    # Checking the current script version.
    ut.check_version(real_path, logger)

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

    logger.info(f"\nWorking with protein {args.name[0]} at {args.temp[0]} K.")

    t0 = time.time()

    # Adding additional arguments to the argument parser.
    vars(args)["proteins"] = proteins
    vars(args)["residues"] = residues
    vars(args)["logger"] = logger
    vars(args)["verbosity"] = verbosity

    # Running main simulation function
    simulate(args)

    logger.info(f"Simulation Done. Total time: {time.time()-t0:.1f} s.")

    # Attempting to send a push notification to a pushbullet account to notify the
    # end of the simulation. Needs a pushbullet API key which has to be stored in the file
    # ./modules/.pbtoken
    try:
        with open(f"{real_path}/modules/.pbtoken", "r") as f:
            pb_token = f.readline().strip()

        notification_body = f"Total time: \n{time.time()-t0:.1f} s\n\nParameters:"
        for arg in vars(args):
            if arg != "residues" and arg != "proteins":
                notification_body += f"\n{arg}: {args.__dict__[arg]}".replace(
                    "[", ""
                ).replace("]", "")

        ut.send_notif(
            title="Simulation Complete",
            body=notification_body,
            pb_token=pb_token,
        )

    except FileNotFoundError:
        logger.info("'.pbtoken' not found. Ommiting push notification.")
