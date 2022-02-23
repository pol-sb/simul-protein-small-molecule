from asyncio.proactor_events import _ProactorDuplexPipeTransport
from fcntl import LOCK_MAND
from matplotlib.pyplot import axis
import numpy as np
import openmm
import openmm.unit as unit
import openmm.app as app
import mdtraj as md
import os
import time


# TODO: Add Charge

def add_drugs(
    system: openmm.openmm.System,
    traj: md.Trajectory,
    conc: float,
    mass: float,
    dist_threshold: float,
    monovalent: bool,
):
    # Decide how many particles to place by using the concentration.
    sim_box_data = system.getDefaultPeriodicBoxVectors()
    sim_box_unit = sim_box_data[0].unit
    len_x = sim_box_data[0][0].value_in_unit(sim_box_unit)
    len_y = sim_box_data[1][1].value_in_unit(sim_box_unit)
    len_z = sim_box_data[2][2].value_in_unit(sim_box_unit)

    sim_box_vol = len_x * len_y * len_z

    # TODO: What happens with the units of mass and concentration?
    # Our length comes in nm, I think should do a conversion.
    print("mass: ", mass)
    print("conc: ", conc)

    # Computing the total number of particles needed to fulfill
    # the given concentration.
    # mass / vol = conc -> mass = conc*vol
    tot_mass = conc * sim_box_vol
    num_part = int(tot_mass / mass)
    print("tot_mass: ", tot_mass)
    print("num_part: ", num_part)


    # Generating random coordinates array with the small particle random
    # coordinates
    if monovalent:
        coor_arr = gen_monovalent(len_x, len_y, len_z, num_part)
    

    print(f"\n{num_part} small molec. generated. Computing collisions...")

    # TODO: Now I should check if the particles collide or are close to any of
    # the particles in our given trajectory
    # Maybe check the distance and if it is lower than a given threshold compute
    # it again.
    # np.linalg.norm(b-a) -> distance

    #######
    # TODO: Optimize this block
    #
    t1 = time.time()
    changes = 0
    for aa_ind, aa_cord in enumerate(traj.xyz[0]):
        print(
            f"Checking collisions with residue: {aa_ind+1} out of"
            f" {len(traj.xyz[0][:,0])} - Changes: {changes}",
            end="\r",
        )
        for drug_ind, drug_coord in enumerate(coor_arr):
            dist = np.linalg.norm(drug_coord - aa_cord)
            while dist <= dist_threshold:
                coor_arr[drug_ind] = (len_x - 0) * np.random.ranf(3) + 0
                dist = np.linalg.norm(drug_coord - aa_cord)
                changes += 1

    t2 = time.time()
    print(f"\nCollision check done. Elapsed time: {t2-t1:.2f} s.")
    #
    #######

    for part in range(num_part):
        system.addParticle(mass)

    

    # It seems that I need a topology for the small molecules
    # in order to be able to merge both trajectories.

    drug_top = md.Topology()
    drug_chain = drug_top.add_chain()
    drug_top.add_residue('TEST', drug_chain)
    resid = drug_top.residue(0)
    for m_ind, molec in enumerate(coor_arr):
        drug_top.add_atom(f"A{m_ind}", None, resid)

    drug_traj = md.Trajectory(coor_arr, drug_top)
    traj = traj.stack(drug_traj)

    # TODO: This returns a slab with the small molecules on it,
    # but won't run the simulation.

    traj.save('test_geom.xyz')

    return traj


def gen_monovalent(L_x, L_y, L_z, num_part):
    x_coords = (L_x - 0) * np.random.ranf(int(num_part)) + 0
    y_coords = (L_y - 0) * np.random.ranf(int(num_part)) + 0
    z_coords = (L_z - 0) * np.random.ranf(int(num_part)) + 0
    coor_arr = np.vstack((x_coords, y_coords, z_coords)).T
    return coor_arr

if __name__ == "__main__":

    print(
        "This file is intented to be used as a module file and should only"
        " called from other scripts, not run by itself."
    )
