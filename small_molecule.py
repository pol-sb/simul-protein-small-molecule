import numpy as np
import openmm
import pandas as pd
import pprint as pp
import openmm.unit as unit
import openmm.app as app
import mdtraj as md
import time


# TODO: Make it so that the added particles dont collide with the residues after
# being placed.

# TODO: Make the system bivalent

# TODO: Use avogadro number for particle collocation!!

# TODO: Add Charge


def add_drugs(
    system: openmm.openmm.System,
    in_top: md.Topology,
    in_traj: md.Trajectory,
    conc: float,
    mass: float,
    dist_threshold: float,
    drug_components: int,
    directory: str,
    components_distance=0.5,
):
    # Decide how many particles to place by using the concentration.
    sim_box_data = system.getDefaultPeriodicBoxVectors()
    sim_box_unit = sim_box_data[0].unit
    len_x = sim_box_data[0][0].value_in_unit(sim_box_unit)
    len_y = sim_box_data[1][1].value_in_unit(sim_box_unit)
    len_z = sim_box_data[2][2].value_in_unit(sim_box_unit)

    sim_box_vol = len_x * len_y * len_z

    # print("mass: ", mass)
    # print("conc: ", conc)

    # TODO: Check this
    # Computing the total number of particles needed to fulfill
    # the given concentration.
    # mass / vol = conc -> mass = conc*vol
    tot_mass = conc * sim_box_vol
    num_part = int(tot_mass / mass)


    # Generating random coordinates array with the small particle random
    # coordinates
    # if drug_components == 1:

    # coor_arr = gen_monovalent(len_x, len_y, len_z, num_part)

    drugs = CGdrug(
        n_drugs=num_part,
        n_components=2,
        system=system,
        in_top=in_top,
        in_traj=in_traj,
        dist_threshold=1,
    )
    drugs.add_charge(-1, 0)

    traj = drugs.get_trajectory()
    top = drugs.get_topology()

    traj.save_pdb(directory+"sm_drg_traj.pdb")
    top_df = top.to_dataframe()
    top_df.data_to_csv(directory+"sm_drg_top.csv")

    return traj, top


# def gen_monovalent(L_x, L_y, L_z, num_part):
#     x_coords = (L_x / 2 - (-L_x / 2)) * np.random.ranf(int(num_part)) + (
#         -L_x / 2
#     )
#     y_coords = (L_y / 2 - (-L_y / 2)) * np.random.ranf(int(num_part)) + (
#         -L_y / 2
#     )
#     z_coords = (L_z - 0) * np.random.ranf(int(num_part)) + 0
#     coor_arr = np.vstack((x_coords, y_coords, z_coords)).T
#     return coor_arr


# def add_polyvalent(coor_arr, num_comp, distance):

#     x_coords = coor_arr[:, 0]
#     y_coords = coor_arr[:, 1]
#     z_coords = coor_arr[:, 2]

#     # TODO: This adds the new second molecule of the bivalent compound at the
#     # end of the coord array. Would this work for molecules with more than two
#     # components?

#     for comp in range(1, num_comp):
#         new_xcoord = x_coords + distance
#         new_ycoord = y_coords
#         new_zcoord = z_coords

#     new_coor_arr = np.vstack((new_xcoord, new_ycoord, new_zcoord)).T

#     final_arr = np.vstack((coor_arr, new_coor_arr))

#     return final_arr


class CGdrug:
    def __init__(
        self,
        n_drugs,
        n_components,
        system,
        in_traj,
        in_top,
        dist_threshold,
        distance=0.5,
    ):

        self.n_drugs = n_drugs
        self.n_components = n_components
        self.tot_atoms = n_drugs * n_components
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

        x_coords = (self.Lx / 2 - (-self.Lx / 2)) * np.random.ranf(
            int(n_drugs)
        ) + (-self.Lx / 2)
        y_coords = (self.Ly / 2 - (-self.Ly / 2)) * np.random.ranf(
            int(n_drugs)
        ) + (-self.Ly / 2)
        z_coords = (self.Lz - (0.0)) * np.random.ranf(int(n_drugs)) + (0.0)
        centers = np.vstack((x_coords, y_coords, z_coords)).T

        return centers

    def _add_components(self, centers, n_comp, distance, dist_threshold):

        # For now n_comp HAS TO BE 2 or this function will EXPLODE!

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
        coord_arr = np.array(coord_arr)

        print(
            "Small drug particles generated. Starting protein collision"
            " check.\n"
        )

        t1 = time.time()
        changes = 0
        for aa_ind, aa_cord in enumerate(self.prot_traj.xyz[0]):

            print(
                f"Checking collisions with residue: {aa_ind+1} out of"
                f" {len(self.prot_traj.xyz[0][:,0])} - Changes: {changes}",
                end="\r",
            )

            for drug_ind, drug_coord in enumerate(self.description.items()):

                for part in drug_coord[1]["coordinates"]:
                    dist = np.linalg.norm(part - aa_cord)

                    while dist <= dist_threshold:
                        new_y_pos = (
                            self.Ly / 2 - (-self.Ly / 2)
                        ) * np.random.ranf(1) + (-self.Ly / 2)

                        # This won't work if we add more than two particles...
                        drug_coord[1]["coordinates"][0][1] = new_y_pos
                        drug_coord[1]["coordinates"][1][1] = new_y_pos

                        dist = np.linalg.norm(part - aa_cord)
                        changes += 1

        t2 = time.time()
        print(f"\nCollision check done. Elapsed time: {t2-t1:.2f} s.")

    def _create_topology(self):

        self.top = self.prot_top

        chain = self.top.add_chain()
        resname = "DRG"

        for drug in range(self.tot_atoms):
            residue = self.top.add_residue(resname, chain)
            self.top.add_atom(
                resname, element=md.element.carbon, residue=residue
            )

        cntr = 0
        for i in range(self.n_drugs):
            self.top.add_bond(
                chain.atom(i + cntr), chain.atom((i + 1) + cntr), order=1
            )
            cntr += 1

    def _create_trajectory(self):

        coord_arr = []

        for drug in self.description.values():
            for coord in drug["coordinates"]:
                coord_arr.append(coord)

        coord_arr = np.array(coord_arr)
        total_coords = np.vstack((self.prot_traj.xyz[0], coord_arr))

        self.trajectory = md.Trajectory(
            total_coords,
            self.top,
            unitcell_lengths=[self.Lx, self.Ly, self.Lz],
            unitcell_angles=[90.0, 90.0, 90.0],
        )

        self.trajectory.save_pdb("traj_test2.pdb")

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
