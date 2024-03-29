import os
import sys

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np

import utils as ut

ASP_RAT_MOD = 1.8
0
Z_LIM = 75
EXTENT = 75
EQUIL_TIME = 2000

# plt.rc('text', usetex=True)
plt.rc("font", family="serif")
plt.rcParams.update({"font.size": 16})
# plt.tick_params(labelsize=14)
plt.rcParams["axes.linewidth"] = 0.5


def get_trajectory(stride):
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
        stride=stride,
    )[0:]
    print(" - DONE")

    return t, sim_name


def create_histogram(z, dens, z_lim):

    # z_lim = max([abs(np.min(z)), abs(np.max(z))])*1.25

    print("\n  [*] Creating histogram", end="")
    for zi in z:
        dens.append(np.histogram(zi, bins=np.arange(-z_lim, z_lim))[0])
    print(" - DONE")
    return dens


def read_lambda():

    param_dict = ut.read_parameters()
    lambd = param_dict["DRG_LAMB"].replace("[", "").replace("]", "")

    return lambd


def prepare_profile(traj, sim_name, Z_LIM, EXTENT):

    z = traj.xyz[:, :, 2]
    box_z_length = traj.unitcell_lengths[0][2]
    dens = []
    dens_drg = []

    # Check if the last chain is the DRG chain
    # assert np.all(traj.top.select("chainid 100") == traj.top.select("resname DRG"))

    drg = traj.top.select("resname DRG")

    # Systems with no small molecule
    if len(drg) == 0:

        print("\n  [*] Centering trajectory", end="")

        z -= np.median(z, 1, keepdims=True)
        z = (z + box_z_length / 2) % box_z_length - box_z_length / 2

        print(" - DONE")

        # S = traj.unitcell_lengths[0, 0] ** 2

        plot_flag = ""

        while plot_flag != "ok":

            dens = []

            dens = create_histogram(z, dens, z_lim=Z_LIM)
            dens = np.array(dens)

            dens_filename = f"avgdens_frame_prot_{sim_name[0][:-11]}_NODRG.out"

            dens_ext_bottom = dens.shape[0]

            # Saving density average for each frame into a file so it can
            # be used later for plotting.
            dens_avg = np.mean(dens[EQUIL_TIME:, :], axis=0)
            dens_avg_dat = np.array([range(dens[0].shape[0]), dens_avg])
            np.savetxt(dens_filename, dens_avg_dat.T)

            plt.xlabel("z", fontsize=25)
            plt.ylabel("Time", fontsize=25)
            plt.title(sim_name[0][:-11] + "_prot")
            plt.imshow(
                dens,
                aspect=(dens.shape[1] / dens.shape[0]) * ASP_RAT_MOD,
                cmap="viridis",
                extent=[-EXTENT, EXTENT, dens_ext_bottom, 0],
            )
            fig.colorbar(im1, shrink=0.5, pad=0.025)
            plt.tight_layout()
            plt.ion()
            plt.show()
            plot_flag = input("\n  [?] Input 'ok' if the plot is ready: ").lower()

            if plot_flag != "ok":
                Z_LIM = float(
                    input(f"\n  [?] Please input 'z_lim' value (previous {Z_LIM}): ")
                )

                EXTENT = float(
                    input(
                        f"\n  [?] Please input the 'extent' value (previous {EXTENT}): "
                    )
                )

        plt.savefig(f"{sim_name[0][:-11]}_NODRG_densprof.png", dpi=400)

    # For systems with small molecule
    else:

        lmbda_val = read_lambda()

        print("\n  [*] Centering trajectory", end="")

        z_prot = z[:, 0 : drg[0]]
        z_drg = z[:, drg]

        # Computing the median which will be used to center our plots. The
        # median is used as it is more robust against outliers than the mean.
        z_prot_mean = np.median(z_prot, 1, keepdims=True)
        z -= z_prot_mean
        z = (z + box_z_length / 2) % box_z_length - box_z_length / 2

        z_prot = z[:, 0 : drg[0]]
        z_drg = z[:, drg]

        print(" - DONE")

        # S = traj.unitcell_lengths[0, 0] ** 2

        plot_flag = ""

        while plot_flag != "ok":

            dens = []
            dens_drg = []

            dens = create_histogram(z_prot, dens, z_lim=Z_LIM)
            dens_drg = create_histogram(z_drg, dens_drg, z_lim=Z_LIM)

            dens = np.array(dens)
            dens_drg = np.array(dens_drg)

            dens_filename = (
                f"avgdens_frame_prot_{sim_name[0][:-11]}_lmbda-"
                f'{lmbda_val.replace(".","")}.out'
            )

            dens_drg_filename = (
                f"avgdens_frame_drg_{sim_name[0][:-11]}_lmbda-"
                f'{lmbda_val.replace(".","")}.out'
            )

            dens_ext_bottom = dens.shape[0]

            # Saving the protein density average for each frame into a file so it
            # can be used later for plotting.
            dens_avg = np.mean(dens[EQUIL_TIME:, :], axis=0)
            dens_avg_dat = np.array([range(dens[0].shape[0]), dens_avg])
            np.savetxt(dens_filename, dens_avg_dat.T)

            # Saving the drug density average for each frame into a file so it can
            # be used later for plotting.
            dens_drg_avg = np.mean(dens_drg[EQUIL_TIME:, :], axis=0)
            dens_drg_avg_dat = np.array([range(dens_drg[0].shape[0]), dens_drg_avg])
            np.savetxt(dens_drg_filename, dens_drg_avg_dat.T)

            # ax = plt.subplots(122)
            # ax2 = plt.subplots(122)
            # # fig2, ax2 = plt.subplots(2,1)
            fig, (ax1, ax2), = plt.subplots(nrows=1, ncols=2)

            plt_prot = ax1
            plt_drg = ax2

            # print(axd)
            # quit()

            plt_prot.set_xlabel("z", fontsize=25)
            plt_prot.set_ylabel("Time", fontsize=25)
            plt_prot.tick_params(axis='x', labelsize=20)
            plt_prot.tick_params(axis='y', labelsize=20)
            plt_prot.set_title(sim_name[0][:-11] + "_prot")
            im1 = plt_prot.imshow(
                dens,
                aspect=(dens.shape[1] / dens.shape[0]) * ASP_RAT_MOD,
                cmap="viridis",
                # vmax=abs(dens).max(),
                # vmin=-abs(dens).max(),
                extent=[-EXTENT, EXTENT, dens_ext_bottom, 0],
            )

            plt_drg.set_xlabel("z", fontsize=25)
            plt_drg.set_ylabel("Time", fontsize=25)
            plt_drg.tick_params(axis='x', labelsize=20)
            plt_drg.tick_params(axis='y', labelsize=20)
            plt_drg.set_title(
                sim_name[0][:-11] + "_drg\n" + r"$\lambda =$" + f"{lmbda_val}"
            )

            im2 = plt_drg.imshow(
                dens_drg,
                aspect=(dens.shape[1] / dens.shape[0]) * ASP_RAT_MOD,
                cmap="viridis",
                extent=[-EXTENT, EXTENT, dens_ext_bottom, 0],
            )
            fig.colorbar(im1, ax=plt_prot, shrink=0.45, pad=0.135)
            fig.colorbar(im2, ax=plt_drg, shrink=0.45, pad=0.135)
            fig.tight_layout(pad=0.4, w_pad=0.2, h_pad=0.2)
            plt.ion()
            plt.show()
            plot_flag = input("\n  [?] Input 'ok' if the plot is ready: ").lower()

            if plot_flag != "ok":
                Z_LIM = float(
                    input(f"\n  [?] Please input 'z_lim' value (previous {Z_LIM}): ")
                )

                EXTENT = float(
                    input(
                        f"\n  [?] Please input the 'extent' value (previous {EXTENT}): "
                    )
                )
            else:
                plt.savefig(
                    f'{sim_name[0][:-11]}_lmbda-{lmbda_val.replace(".","")}_densprof.png',
                    dpi=1000,
                )
                plt.savefig(
                    f'{sim_name[0][:-11]}_lmbda-{lmbda_val.replace(".","")}_densprof.svg',
                    dpi=1000,
                )


if __name__ == "__main__":
    traj, s_name = get_trajectory(stride=None)
    prepare_profile(traj, s_name, Z_LIM, EXTENT)
