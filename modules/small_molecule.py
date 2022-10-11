import logging
import pprint as pp
import time

import mdtraj as md
import numpy as np
import openmm
import openmm.app as app
import openmm.unit as unit
import pandas as pd

# TODO: Add Charge

logger = logging.getLogger("small molecule")


def add_drugs(
    system: openmm.openmm.System,
    in_top: md.Topology,
    in_traj: md.Trajectory,
    conc: float,
    dist_threshold: float,
    drug_components: list,
    directory: str,
    comp_dist,
    verbosity,
    residues,
    col_chk_flag: bool,
    lambd_override: float,
    mass_override: float,
    sigma_override: float,
):
    # Get verbosity level
    # logger.setLevel(verbosity[0])

    # This list will contain a list of small molecule related parameters
    # for result reporting.
    drg_param = []

    # Decide how many particles to place by using the concentration.
    sim_box_data = system.getDefaultPeriodicBoxVectors()
    sim_box_unit = sim_box_data[0].unit
    len_x = sim_box_data[0][0].value_in_unit(sim_box_unit)
    len_y = sim_box_data[1][1].value_in_unit(sim_box_unit)
    len_z = sim_box_data[2][2].value_in_unit(sim_box_unit)

    # In nanometers
    sim_box_vol = len_x * len_y * len_z

    # Converting the volume to liter
    sim_box_vol_L = sim_box_vol * 1e-24
    logger.debug(f"Simulation box volume: {sim_box_vol_L} L")
    logger.debug(f"Small molecule concentration: {conc}")

    lambda_list = []

    # TODO: Update this for several components.
    if lambd_override or lambd_override == 0:
        drg_param.append(lambd_override)
    else:
        for drg_typ in drug_components:
            lam = residues.loc[residues["three"] == f"{drg_typ}"]["lambdas"][0]
            lambda_list.append(lam)
        drg_param.append(lambda_list)

    # Computing the total number of particles needed to fulfill
    # the given concentration.
    # mass / vol = conc -> mass = conc*vol
    # conc is in mmol/L
    # sim_box_vol_L in L
    # therefore, tot_mmol is mmol
    tot_mmol = conc * sim_box_vol_L
    logger.debug(f"Total mmols: {tot_mmol} mmol")

    # Converting the total mmol to mol, and then to particles using the
    # Avogadro's constant.
    num_part = int((tot_mmol * 1e-3) * 6.02214076e23)
    logger.info(f"Generating {num_part} particles.")
    drg_param.append(num_part)

    sigmas = []
    if sigma_override or sigma_override == 0:
        sigmas.append(sigma_override[0])
    else:
        for drg_typ in drug_components:
            sig = residues.loc[residues["three"] == f"{drg_typ}"]["sigmas"][0]
            sigmas.append(sig)

    # Computing the total volume occupied by all of the small molecules.
    # BUG: Volume is multiplied every time by num part, which is not correct if there
    # is more than one sigma given
    volume = 0
    for sigma in sigmas:
        volume += ((4 / 3) * np.pi * ((sigma / 2) ** 3)) * num_part
    logger.info(f"Volume* occupied by the small molecules: {round(volume, 2)}")

    # Adding the small molecules to the system. Each AA has its own molecular
    # weight.
    if mass_override:
        for drg_ind, drg_typ in enumerate(drug_components):
            logger.info(
                f"For {drg_typ}-{drg_ind} using non-default mass:"
                f" {mass_override[drg_ind]}"
            )
            for part in range(num_part):
                system.addParticle(mass_override[drg_ind] * unit.amu)
    else:
        for drg_ind, drg_typ in enumerate(drug_components):
            res_row = residues.loc[residues["three"] == f"{drg_typ}"]
            logger.info(
                f"For {drg_typ}-{drg_ind} using default mass: {res_row['MW'][0]}"
            )
            for part in range(num_part):
                system.addParticle(res_row["MW"][0] * unit.amu)

    drugs = CGdrug(
        n_drugs=num_part,
        n_components=drug_components,
        system=system,
        in_top=in_top,
        in_traj=in_traj,
        dist_threshold=25,
        distance=comp_dist,
        verbosity=verbosity,
        col_chk_flag=col_chk_flag,
    )

    # TODO: The add_charge function is unused.
    # drugs.add_charge(-1, 0)

    traj = drugs.get_trajectory()
    top = drugs.get_topology()

    # Saving the Trajectory and Topology so it can be loaded later,
    # for example, if the calculation is stopped and then resumed.
    logger.debug(
        f"Saving small molecule trajectories and topologies in: '{directory}'"
    )
    traj.save_pdb(directory + "sm_drg_traj.pdb")
    top_df = top.to_dataframe()
    top_ats = top_df[0]
    top_bnd = top_df[1]
    top_ats.to_csv(directory + "sm_drg_ats.csv")
    np.save(directory + "sm_drg_bnd.npy", top_bnd)

    return traj, top, system, num_part, drg_param


class CGdrug:
    def __init__(
        self,
        n_drugs,
        n_components,
        system,
        in_traj,
        in_top,
        dist_threshold,
        distance,
        verbosity,
        col_chk_flag,
    ):

        self.verbosity = verbosity
        self.n_drugs = n_drugs
        self.n_components = n_components
        self.tot_atoms = n_drugs * len(n_components)
        self.prot_traj = in_traj
        self.prot_top = in_top

        centers = self._gen_centers(n_drugs, system)
        self._add_components(
            centers, n_components, distance, dist_threshold, col_chk_flag
        )
        self._create_topology()
        self._create_trajectory()

    def _gen_centers(self, n_drugs, sys):
        sim_box_data = sys.getDefaultPeriodicBoxVectors()
        sim_box_unit = sim_box_data[0].unit

        self.Lx = sim_box_data[0][0].value_in_unit(sim_box_unit)
        self.Ly = sim_box_data[1][1].value_in_unit(sim_box_unit)
        self.Lz = sim_box_data[2][2].value_in_unit(sim_box_unit)

        x_coords = (self.Lx / 2 - (-self.Lx / 2)) * np.random.ranf(int(n_drugs)) + (
            -self.Lx / 2
        )
        y_coords = (self.Ly / 2 - (-self.Ly / 2)) * np.random.ranf(int(n_drugs)) + (
            -self.Ly / 2
        )
        z_coords = (self.Lz - (0.0)) * np.random.ranf(int(n_drugs)) + (0.0)
        centers = np.vstack((x_coords, y_coords, z_coords)).T

        return centers

    def collision_check(self, dist_threshold: float, enabled: bool):

        if enabled:
            logger.info(
                "Small drug particles generated. Starting protein collision check..."
            )

            dist_threshold = enabled[0]

            t1 = time.time()
            changes = 0
            for aa_ind, aa_cord in enumerate(self.prot_traj.xyz[0]):

                aa_cord = np.array(aa_cord)

                if self.verbosity[0] == logging.INFO:
                    print(
                        f"Checking collisions with residue: {aa_ind+1} out of"
                        f" {len(self.prot_traj.xyz[0][:,0])} - Changes: {changes}",
                        end="\r",
                    )

                drg_pos_arr = []
                for drg in self.description.items():
                    drg_pos_arr.append(drg[1]["coordinates"])
                drg_pos_arr = np.array(drg_pos_arr)

                dist_arr = np.linalg.norm(drg_pos_arr - aa_cord, axis=1)

                for chg_drg_ind in np.where(dist_arr < dist_threshold)[0]:
                    drg_pos_arr[chg_drg_ind][1] = (
                        self.Ly / 2 - (-self.Ly / 2)
                    ) * np.random.ranf(1) + (-self.Ly / 2)
                    changes += 1

            t2 = time.time()
            logger.info(f"Collision check done. Elapsed time: {t2-t1:.2f} s.")
        else:
            logger.warning("Collision check disabled.")

    def _add_components(self, centers, n_comp, distance, dist_threshold, col_chk_flag):

        # If two components are detected, add components to the drug description
        # objects.
        if len(n_comp) == 2:

            # Printing debug information
            logger.debug(f"Detected 2 components: {n_comp}")
            logger.debug(f"Distance between components: {distance}")

            # Creating description property which will be used inside the class
            self.description = {}

            # Working with each small molecule center
            for ind, center in enumerate(centers):

                # List that will contain the coordinates for a small molecule
                drug_coor = []

                # Computing the positions of each of the two elements in the small
                # molecule and adding them to a list.
                for i in [-1, 1]:
                    comp_coor = [
                        center[0] + ((distance / 2) * i),
                        center[1],
                        center[2],
                    ]
                    comp_coor = np.array(comp_coor)
                    drug_coor.append(comp_coor)

                self.description[ind] = {"coordinates": drug_coor}

            # Creating a new object with all of the small molecule coordinates in
            # an array
            coord_arr = []
            for drug in self.description.values():
                for coord in drug["coordinates"]:
                    coord_arr.append(coord)
            self.coord_arr_drg = np.array(coord_arr)

        # Case for small molecules with just one component
        elif len(n_comp) == 1:

            # Printing debug information.
            logger.debug(f"Detected 1 component: {n_comp}")

            self.description = {}

            # Declaring the position of the single component small molecule.
            for ind, center in enumerate(centers):
                comp_coor = [
                    center[0],
                    center[1],
                    center[2],
                ]

                self.description[ind] = {"coordinates": comp_coor}

            coord_arr = []
            for drug in self.description.values():
                coord_arr.append(drug["coordinates"])
            self.coord_arr_drg = np.array(coord_arr)

        # If a larger than 2 odd number of small molecule residues are given, use this
        # method of assignment.
        elif len(n_comp) > 2 and len(n_comp) % 2 != 0:

            logger.debug(f"User defined small molecules have {len(n_comp)} components.")
            logger.debug(f"Distance between components: {distance}")

            self.description = {}

            # Working with each small molecule center
            for ind, center in enumerate(centers):
                drug_coor = []

                # Here, the odd numbered chain is split into two segments in order to be
                # able to assign the distances differently for each segment.
                # The first segment goes from the initial bead to the center bead, and
                # the second segment goes from the middle bead until the last one.

                # This loop assigns the positions of the beads in the first segment of
                # the small molecule chain.
                for i in range(-(len(n_comp) // 2), 1):
                    comp_coor = [
                        center[0],
                        center[1],
                        center[2] + ((distance) * i),
                    ]
                    comp_coor = np.array(comp_coor)
                    drug_coor.append(comp_coor)

                # This loop assigns the positions of the beads in the second segment of
                # the small molecule chain.
                for i in range(1, (len(n_comp) // 2) + 1):
                    comp_coor = [
                        center[0],
                        center[1],
                        center[2] + ((distance) * i),
                    ]

                    comp_coor = np.array(comp_coor)
                    drug_coor.append(comp_coor)

                self.description[ind] = {"coordinates": drug_coor}

            coord_arr = []
            for drug in self.description.values():
                for coord in drug["coordinates"]:
                    coord_arr.append(coord)
            self.coord_arr_drg = np.array(coord_arr)

        # If a larger than 2 even number of small molecule residues are given, use this
        # method of assignment.
        elif len(n_comp) > 2 and len(n_comp) % 2 == 0:

            logger.debug(f"User defined small molecules have {len(n_comp)} components.")
            logger.debug(f"Distance between components: {distance}")

            self.description = {}
            for ind, center in enumerate(centers):
                drug_coor = []

                # Searching for the bead next to the chain center
                near_middl_bead = center - np.array([0, 0, distance / 2])
                init_bead = near_middl_bead

                # Searching for the initial bead of the chain using the next to middle
                # bead.
                for i in range((len(n_comp) - 1) // 2):
                    init_bead -= np.array([0, 0, distance])

                drug_coor.append(init_bead)

                # Starting from the initial bead, a bond distance is added to every
                # bead's position in the z axis in order to build the chain.
                for i in range(len(n_comp) - 1):
                    comp_coor = [
                        init_bead[0],
                        init_bead[1],
                        init_bead[2] + distance,
                    ]
                    init_bead = comp_coor
                    comp_coor = np.array(comp_coor)
                    drug_coor.append(comp_coor)

                self.description[ind] = {"coordinates": drug_coor}

            coord_arr = []
            for drug in self.description.values():
                for coord in drug["coordinates"]:
                    coord_arr.append(coord)
            self.coord_arr_drg = np.array(coord_arr)

        # After the small molecules are added, the collision check is performed if
        # the correct flag is given by the user at launch.
        self.collision_check(dist_threshold, enabled=col_chk_flag)

    def _create_topology(self):

        self.top = self.prot_top

        chain = self.top.add_chain()
        resname = "DRG"

        for drug in range(self.tot_atoms):
            residue = self.top.add_residue(resname, chain)
            self.top.add_atom(resname, element=md.element.carbon, residue=residue)

        cntr = 0
        if len(self.n_components) == 2:
            for i in range(self.n_drugs):
                self.top.add_bond(
                    chain.atom(i + cntr), chain.atom((i + 1) + cntr), order=1
                )
                cntr += 1
        else:
            logger.critical("[!] Topology not correctly implemented.")
            cntr = 0

            for drug_ind in range(self.n_drugs):
                for drug_comp in range(len(self.n_components)):
                    if cntr < (self.n_drugs*len(self.n_components))-1:

                        self.top.add_bond(
                            chain.atom((drug_ind * len(self.n_components)) + drug_comp),
                            chain.atom(
                                (drug_ind * len(self.n_components)) + (drug_comp + 1)
                            ),
                            order=1,
                        )

                        cntr += 1
                    #  cntr = 0
            # quit()

    def _create_trajectory(self):

        total_coords = np.vstack((self.prot_traj.xyz[0], self.coord_arr_drg))

        self.trajectory = md.Trajectory(
            total_coords,
            self.top,
            unitcell_lengths=[self.Lx, self.Ly, self.Lz],
            unitcell_angles=[90.0, 90.0, 90.0],
        )

    def add_charge(self, charge, comp_index):
        for drug in self.description.values():

            if "charge" in drug:
                drug["charge"][comp_index] = charge
            else:
                drug["charge"] = np.zeros([len(drug["coordinates"])])
                drug["charge"][comp_index] = charge

    def get_topology(self):
        return self.top

    def get_trajectory(self):
        return self.trajectory


if __name__ == "__main__":

    print(
        "This file is intented to be used as a module file and should only"
        " called from other scripts, not run by itself."
    )
