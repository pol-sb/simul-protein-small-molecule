import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np
import os
import sys


def get_trajectory():
    sim_path = os.getcwd()
    sim_name = [f for f in os.listdir(sim_path) if f.endswith("report.dcd")]

    if len(sim_name) == 0:
        print(
            "ERROR -  No trajectory files found. Please make sure you are in the"
            " correct folder."
        )
        sys.exit(1)

    print("\n  [*] Loading trajectory", end="")

    t = md.load(
        f"{sim_path}/{sim_name[0]}",
        top=f"{sim_path}/final_system_state.pdb",
    )[0:]
    print(" - DONE")

    return t, sim_name


def create_histogram(z, dens):

    z_lim = 75
    # z_lim = max([abs(np.min(z)), abs(np.max(z))])*1.25

    print("\n  [*] Creating histogram", end="")
    for zi in z:
        # print('zi: ', max(zi), min(zi))
        dens.append(np.histogram(zi, bins=np.arange(-z_lim, z_lim))[0])
    print(" - DONE")
    return dens


def prepare_profile(traj, sim_name):

    z = traj.xyz[:, :, 2]
    print('zmax: ', z.max())
    print('zmin: ', z.min())
    quit()
    dens = []
    dens_drg = []

    # Check if the last chain is the DRG chain
    # assert np.all(traj.top.select("chainid 100") == traj.top.select("resname DRG"))

    drg = traj.top.select("resname DRG")
    if len(drg) == 0:
        print("\n  [*] Centering trajectory", end="")

        z -= z.mean(1, keepdims=True)
        print(" - DONE")

        S = traj.unitcell_lengths[0, 0] ** 2

        dens = create_histogram(z, dens)

        plt.xlabel("z")
        plt.ylabel("time")
        plt.title(sim_name[0][:-11] + "_prot")
        plt.imshow(dens, aspect=0.1, cmap="plasma")
        plt.tight_layout()
        plt.savefig(f"{sim_name[0][:-11]}_NODRG_densprof.png", dpi=200)

    else:

        lmbda_val = input("      Plase input the lambda value: ")

        print("\n  [*] Centering trajectory", end="")

        z_prot = z[:, 0 : drg[0]]
        z_drg = z[:, drg]
        z_prot_mean = z_prot.mean(1, keepdims=True)
        z -= z_prot_mean
        print(" - DONE")

        S = traj.unitcell_lengths[0, 0] ** 2

        dens = create_histogram(z_prot, dens)

        dens_drg = create_histogram(z_drg, dens_drg)

        # for zi in z_chains:
        #    print(zi.shape)
        #    break

        dens = np.array(dens)

        dens_ext_bottom = dens.shape[0]

        fig, ax = plt.subplots(ncols=2)

        plt_prot = ax[0]
        plt_drg = ax[1]

        plt_prot.set_xlabel("z")
        plt_prot.set_ylabel("time")
        plt_prot.set_title(sim_name[0][:-11] + "_prot")
        a = plt_prot.imshow(
            dens,
            aspect=(dens.shape[1] / dens.shape[0]) * 0.9,
            cmap="plasma",
            # vmax=abs(dens).max(),
            # vmin=-abs(dens).max(),
            extent=[-75, 75, dens_ext_bottom, 0],
        )

        plt_drg.set_xlabel("z")
        plt_drg.set_ylabel("time")
        plt_drg.set_title(
            sim_name[0][:-11] + "_drg\n" + r"$\lambda =$" + f"{lmbda_val}"
        )

        plt_drg.imshow(
            dens_drg,
            aspect=(dens.shape[1] / dens.shape[0]) * 0.9,
            cmap="plasma",
            extent=[-75, 75, dens_ext_bottom, 0],
        )

        plt.tight_layout()
        plt.savefig(
            f'{sim_name[0][:-11]}_lmbda-{lmbda_val.replace(".","")}_densprof.png',
            dpi=200,
        )


if __name__ == "__main__":
    traj, s_name = get_trajectory()
    prepare_profile(traj, s_name)
