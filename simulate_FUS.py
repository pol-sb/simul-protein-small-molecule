import openmm
import openmm.unit as unit
from openmm import app
from openmm import XmlSerializer
from analyse import *
import time
import small_molecule as smol
import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("--name", nargs="?", const="", type=str)
parser.add_argument("--temp", nargs="?", const="", type=int)
args = parser.parse_args()


def simulate(residues, name, prot, temp):

    folder = name + "/{:d}/".format(temp)
    check_point = folder + "restart.chk"

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
    margin = 2
    if N > 400:
        L = 25.0
        Lz = 300.0
        margin = 8
        Nsteps = int(2e7)
    elif N > 200:
        L = 17.0
        Lz = 300.0
        margin = 4
        Nsteps = int(6e7)
    else:
        Lz = 10 * L
        Nsteps = int(6e7)

    system = openmm.System()

    # Set box vectors
    # The following block of code defines the size of each vector of the slab
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
    # TODO: Protein chains are randomly distributed here?
    xy = np.empty(0)

    # Generates a 1x2 vector with two x and y coordinates.
    xy = np.append(
        xy, np.random.rand(2) * (L - margin) - (L - margin) / 2
    ).reshape((-1, 2))

    # Generates 100 new coordinates for the xy array, which must fulfill
    # certain conditions:
    #
    #
    for x, y in np.random.rand(1000, 2) * (L - margin) - (L - margin) / 2:
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

    # Saving a .pdb file with the current configuration. This file has 100
    # protein strands spanning a long stretch of the z axis and completely
    # straight
    in_traj = md.Trajectory(
        np.array(pos).reshape(n_chains * N, 3),
        top,
        0,
        [L, L, Lz],
        [90, 90, 90],
    )

    in_traj.save_pdb(name + "/{:d}/top.pdb".format(temp))
    pdb = app.pdbfile.PDBFile(name + "/{:d}/top.pdb".format(temp))

    # Adding finally the particles to the system, with their charge and a
    # term for compensating for the terminal residues (2 for a terminal NH2
    # and 16 for the deprotonated OH of the ending carbonyl group)
    for _ in range(n_chains):
        system.addParticle((residues.loc[prot.fasta[0]].MW + 2) * unit.amu)
        for a in prot.fasta[1:-1]:
            system.addParticle(residues.loc[a].MW * unit.amu)
        system.addParticle((residues.loc[prot.fasta[-1]].MW + 16) * unit.amu)

    print(in_traj.xyz.shape)
    print(top)

    if not os.path.isfile(check_point):
        print("\nAdding small molecules to the system...")
        # My function to add small particles to the system
        in_traj, top, system, n_drugs = smol.add_drugs(
            system=system,
            in_traj=in_traj,
            in_top=top,
            conc=0.005,
            mass=1,
            dist_threshold=2,
            drug_components=2,
            directory=folder,
        )

    else:
        print("\nReading small molecules from stored files...")
        top_ats = pd.read_csv(folder + "sm_drg_ats.csv")
        top_bnd = np.load(folder + "sm_drg_bnd.npy")
        top = md.Topology.from_dataframe(top_ats, top_bnd)

        in_traj = md.load(folder + "sm_drg_traj.pdb")

        # print(in_traj.xyz.shape)
        # print(top)

    # TODO: This is giving me a problem, the number of particles in the pdb
    # is 16636 (larger than 16468 which is init system + drugs) and results in
    # an error.
    pdb = app.pdbfile.PDBFile(folder + "sm_drg_traj.pdb")
    print("len pdb: ", pdb.getTopology().getNumAtoms())

    #######
    # TODO: Add function or block of code to add coarse-grained chemical
    # compounds here.
    #
    # Initially they must be monovalent compounds (just one particle) without
    # charge, but this should be able to be changed. We will do bivalent
    # compounds with charge eventually.
    #
    # How many particles will be set with the concentration MUST be able
    # to be changed to make several tests.
    #######

    # Adding a regular harmonic bond force
    hb = openmm.openmm.HarmonicBondForce()

    # Adding our custom energy expression for a LJ type interaction to use in
    # the FF. It is non bonded as LJ is not bonded. Ashbaugh-Hatch (ah)
    # functional form.
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

    print("Number of particles:", system.getNumParticles())
    print("Number of drugs:", n_drugs)

    for j in range(n_chains):
        begin = j * N
        end = j * N + N

        for a, e in zip(prot.fasta, yukawa_eps):
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

            # Important
            yu.addExclusion(i, i + 1)
            ah.addExclusion(i, i + 1)

    # Adding the small drug particles to the CustomNonbondedForce used in
    # # the system. n_drugs (number of small molecules) by 2 (bimolecular).
    for i in range(n_drugs * 2):

        # Yukawa Epsilon of the small molecules
        yu.addParticle([1 * unit.nanometer * unit.kilojoules_per_mole])
        ah.addParticle(
            [
                # Sigma of the small molecules.
                1 * unit.nanometer,
                # Lambda of the small molecules.
                1 * unit.dimensionless,
            ]
        )

    # TODO: I will have to add bonds to my particles!! AAAA

    print("ah:", ah.getNumParticles())
    print("yu:", yu.getNumParticles())

    yu.setNonbondedMethod(openmm.openmm.CustomNonbondedForce.CutoffPeriodic)
    ah.setNonbondedMethod(openmm.openmm.CustomNonbondedForce.CutoffPeriodic)
    hb.setUsesPeriodicBoundaryConditions(True)
    # tambe amunt
    yu.setCutoffDistance(4 * unit.nanometer)
    ah.setCutoffDistance(4 * unit.nanometer)

    system.addForce(hb)
    system.addForce(yu)
    system.addForce(ah)

    serialized_system = XmlSerializer.serialize(system)
    outfile = open("system.xml", "w")
    outfile.write(serialized_system)
    outfile.close()

    integrator = openmm.openmm.LangevinIntegrator(
        temp * unit.kelvin, 0.01 / unit.picosecond, 0.005 * unit.picosecond
    )  # 322

    # Uses CUDA as the platform for the GPU calculations.
    platform = openmm.Platform.getPlatformByName("CPU")

    simulation = app.simulation.Simulation(
        pdb.topology,
        system,
        integrator,
        platform,
    )

    if os.path.isfile(check_point):
        print("\nResuming simulation from checkpoint file\n")
        simulation.loadCheckpoint(check_point)
        simulation.reporters.append(
            app.dcdreporter.DCDReporter(
                name + "/{:d}/{:s}.dcd".format(temp, name),
                int(10),
                append=True,
            )
        )
    else:
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()
        simulation.reporters.append(
            app.dcdreporter.DCDReporter(
                name + "/{:d}/{:s}.dcd".format(temp, name), int(10)
            )
        )

    # Generates log file with information
    simulation.reporters.append(
        app.statedatareporter.StateDataReporter(
            "{:s}_{:d}.log".format(name, temp),
            10,
            potentialEnergy=True,
            temperature=True,
            step=True,
            speed=True,
            elapsedTime=True,
            separator="\t",
        )
    )

    # TODO: Defines total runtime (20h) and checkpoint save interval (5h)?
    simulation.runForClockTime(
        2 * unit.minute,
        checkpointFile=check_point,
        checkpointInterval=30 * unit.second,
    )

    # Saves checkpoint file
    simulation.saveCheckpoint(check_point)

    # TODO: What is this.
    # genDCD(residues, name, prot, temp, n_chains)
    positions = simulation.context.getState(getPositions=True).getPositions()
    app.PDBFile.writeFile(
        simulation.topology,
        positions,
        open(f"{folder}final_system_state.pdb", "w"),
    )


residues = pd.read_csv("residues.csv").set_index("three", drop=False)
proteins = pd.read_pickle("proteins.pkl")
print(args.name, args.temp)
t0 = time.time()
simulate(residues, args.name, proteins.loc[args.name], args.temp)
print("Simulation Done. Total time: {:.3f}".format(time.time() - t0))
