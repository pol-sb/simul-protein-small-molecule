import logging
import pprint as pp
import time

import mdtraj as md
import numpy as np
import openmm
import openmm.app as app
import openmm.unit as unit
import pandas as pd

# TODO: Make it so that the added particles dont collide with the residues after
# being placed.

# TODO: Make the system bivalent

# TODO: Use avogadro number for particle collocation!!

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
):
    # Get verbosity level
    # logger.setLevel(verbosity[0])

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

    # comp_molarmass = 0
    # for drg_typ in drug_components:
    #     mmass = residues.loc[residues["three"] == f"{drg_typ}"]["MW"][0]
    #     comp_molarmass += mmass

    # print("comp_molarmass: ", comp_molarmass, "g/mol")

    # TODO: Check this
    # Computing the total number of particles needed to fulfill
    # the given concentration.
    # mass / vol = conc -> mass = conc*vol
    # conc is in mmol/L
    # sim_box_vol_L in L
    # therefore, tot_mmol is mmol
    tot_mmol = conc * sim_box_vol_L
    print(f"Total mmols: {tot_mmol} mmol")

    # Converting the total mmol to mol, and then to particles using the
    # Avogadro's constant.
    num_part = int((tot_mmol * 1e-3) * 6.02214076e23)
    logger.debug(f"Number of particles:{num_part} particles.")

    # Adding the small molecules to the system. Each AA has its own molecular
    # weight.
    for drg_typ in drug_components:
        res_row = residues.loc[residues["three"] == f"{drg_typ}"]
        for part in range(num_part):
            system.addParticle(res_row["MW"][0] * unit.amu)

    drugs = CGdrug(
        n_drugs=num_part,
        n_components=drug_components,
        system=system,
        in_top=in_top,
        in_traj=in_traj,
        dist_threshold=1,
        distance=comp_dist,
        verbosity=verbosity,
    )

    # TODO: The add_charge function is unused.
    # drugs.add_charge(-1, 0)

    traj = drugs.get_trajectory()
    top = drugs.get_topology()

    # Saving the Trajectory and Topology so it can be loaded later,
    # for example, if the calculation is stopped and then resumed.
    logger.debug(
        f"\nSaving small molecule trajectories and topologies in:\n'{directory}'"
    )
    traj.save_pdb(directory + "sm_drg_traj.pdb")
    top_df = top.to_dataframe()
    top_ats = top_df[0]
    top_bnd = top_df[1]
    top_ats.to_csv(directory + "sm_drg_ats.csv")
    np.save(directory + "sm_drg_bnd.npy", top_bnd)

    return traj, top, system, num_part


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
    ):

        self.verbosity = verbosity
        self.n_drugs = n_drugs
        self.n_components = n_components
        self.tot_atoms = n_drugs * len(n_components)
        self.prot_traj = in_traj
        self.prot_top = in_top

        centers = self._gen_centers(n_drugs, system)
        self._add_components(centers, n_components, distance, dist_threshold)
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

    def collision_check(self, dist_threshold):

        if self.verbosity[0] == logging.DEBUG:
            verb = True
        else:
            verb = False

        logger.info(
            "Small drug particles generated. Starting protein collision check...\n"
        )

        t1 = time.time()
        changes = 0
        for aa_ind, aa_cord in enumerate(self.prot_traj.xyz[0]):

            if verb:
                print(
                    f"Checking collisions with residue: {aa_ind+1} out of"
                    f" {len(self.prot_traj.xyz[0][:,0])} - Changes: {changes}",
                    end="\r",
                )

            for drug_ind, drug_coord in enumerate(self.description.items()):

                for part in drug_coord[1]["coordinates"]:
                    dist = np.linalg.norm(part - aa_cord)

                    while dist <= dist_threshold:
                        new_y_pos = (self.Ly / 2 - (-self.Ly / 2)) * np.random.ranf(
                            1
                        ) + (-self.Ly / 2)

                        # This won't work if we add more than two particles...
                        drug_coord[1]["coordinates"][0][1] = new_y_pos
                        drug_coord[1]["coordinates"][1][1] = new_y_pos

                        dist = np.linalg.norm(part - aa_cord)
                        changes += 1

        t2 = time.time()
        logger.info(f"\nCollision check done. Elapsed time: {t2-t1:.2f} s.")

    def _add_components(
        self,
        centers,
        n_comp,
        distance,
        dist_threshold,
    ):

        if len(n_comp) == 2:
            logger.debug(f"Detected 2 components: {n_comp}")
            logger.debug(f"Distance between components: {distance}")
            self.description = {}
            for ind, center in enumerate(centers):
                drug_coor = []
                for i in [-1, 1]:
                    comp_coor = [
                        center[0] + ((distance / 2) * i),
                        center[1],
                        center[2],
                    ]
                    comp_coor = np.array(comp_coor)
                    drug_coor.append(comp_coor)

                self.description[ind] = {"coordinates": drug_coor}

            coord_arr = []
            for drug in self.description.values():
                for coord in drug["coordinates"]:
                    coord_arr.append(coord)
            self.coord_arr_drg = np.array(coord_arr)

        elif len(n_comp) == 1:
            logger.debug(f"Detected 1 component: {n_comp}")

            self.description = {}
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

        else:
            logger.critical(f"{len(n_comp)} components not supported.")

        # logger.critical("Remember to delete this and reenable the collision check")
        self.collision_check(dist_threshold)

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

    def _create_trajectory(self):

        # coord_arr = []

        # for drug in self.description.values():
        #     for coord in drug["coordinates"]:
        #         coord_arr.append(coord)

        # coord_arr = np.array(coord_arr)

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
