import itertools
import os
import pprint

import matplotlib.pyplot as plt
import numpy as np

import utils as ut

import traceback

# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
plt.rcParams.update({"font.size": 14})
plt.rcParams["axes.linewidth"] = 0.6
# plt.style.use('grayscale')


def choose_cmap_color(value, prop):

    TEMPS = [
        # 290,
        #300,
        # 310,
        323,
        # 330,
        #340,
    ]
    LAMBDAS = [
        0.0,
        0.175,
        0.35,
        0.525,
        0.7,
    ]
    CONC = [
        0.2,
        2.0,
        #10,
        20.0,
    ]

    SIGMA = [0.45, 0.97, 2.099]
    # print("conc:", CONC)
    # print("type conc:", type(CONC))

    if prop in ["temperature", "T", "temp", "t"]:
        col_arr = plt.cycler("color", plt.cm.inferno(np.linspace(0, 1, len(TEMPS))))
        col_val = np.array(col_arr)[TEMPS.index(value)]

    elif prop in ["lambda", "l", "lam", "lamb"]:
        col_arr = plt.cycler("color", plt.cm.plasma(np.linspace(0, 1, len(LAMBDAS))))
        col_val = np.array(col_arr)[LAMBDAS.index(value)]

    elif prop in ["conc", "c", "conc"]:
        col_arr = plt.cycler("color", plt.cm.winter(np.linspace(0, 1, len(CONC))))
        col_val = np.array(col_arr)[CONC.index(value)]

    elif prop in ["both"]:
        combin = list(itertools.product(TEMPS, LAMBDAS))
        # print("combin: ", combin)
        col_arr = plt.cycler("color", plt.cm.magma(np.linspace(0, 1, len(combin))))
        col_val = np.array(col_arr)[combin.index(value)]

    elif prop in ["conc_temp"]:
        combin = list(itertools.product(TEMPS, CONC))
        col_arr = plt.cycler("color", plt.cm.cool(np.linspace(0, 1, len(combin))))
        col_val = np.array(col_arr)[combin.index(value)]

    elif prop in ["conc_lambda"]:
        combin = list(itertools.product(LAMBDAS, CONC))
        #print("combin: ", combin)
        col_arr = plt.cycler("color", plt.cm.spring(np.linspace(0, 1, len(combin))))
        col_val = np.array(col_arr)[combin.index(value)]

    elif prop in ["conc_sigma", "csig"]:
        combin = list(itertools.product(SIGMA, CONC))
        #print("combin: ", combin)
        col_arr = plt.cycler("color", plt.cm.viridis(np.linspace(0, 1, len(combin))))
        col_val = np.array(col_arr)[combin.index(value)]

    elif prop in ["all"]:
        combin = list(itertools.product(TEMPS, LAMBDAS, CONC))
        col_arr = plt.cycler("color", plt.cm.prism(np.linspace(0, 1, len(combin))))
        col_val = np.array(col_arr)[combin.index(value)]

    return col_val["color"]


def plot_densities(file_paths, color: str, conc_units: str):

    fig, (ax1, ax2) = plt.subplots(ncols=2)
    fig.set_size_inches(12, 6)

    for folder in file_paths:

        params = ut.read_parameters(folder)

        print("params: ", params)
        # quit()

        T = params["TEMP_(K)"]
        d = params["DRG_NAME"]
        p = params["PROT_NAME"]
        lam = params["DRG_LAMB"].replace("[", "").replace("]", "")
        cc = params["DRG_CONC_(mM)"]

        if "DRG_SIGMA" in params.keys():
            sigma = float(params["DRG_SIGMA"])
        else:
            sigma = float(input("Please input sigma for the current simulation:"))

        drg_volume = ((4 / 3) * np.pi * ((sigma / 2) ** 3)) * int(params["DRG_NUMB"])

        # In nanometers
        sim_box_vol = 15 * 15 * 150

        # Converting the simbox volume to liters
        sim_box_vol_L = sim_box_vol * 1e-24
        drg_L = drg_volume / sim_box_vol_L

        f_path = [f for f in os.listdir(folder) if f.endswith(".out") and "prot" in f]

        arr = np.loadtxt(folder + "/" + f_path[0])

        if color in ["lambda", "lambdas", "lam"]:
            c_R, c_G, c_B, c_A = choose_cmap_color(lam, "lambda")
        elif color in ["temp", "temperature", "T"]:
            c_R, c_G, c_B, c_A = choose_cmap_color(int(T), "temp")
        elif color in ["both"]:
            c_R, c_G, c_B, c_A = choose_cmap_color((T, lam), "both")
        elif color in ["conc", "c"]:
            c_R, c_G, c_B, c_A = choose_cmap_color(float(cc), "conc")
        elif color in ["conc_temp", "ct"]:
            c_R, c_G, c_B, c_A = choose_cmap_color((int(T), float(cc)), "conc_temp")
        elif color in ["conc_sigma", "csig"]:
            c_R, c_G, c_B, c_A = choose_cmap_color((float(sigma), float(cc)), "csig")
        elif color in ["conc_lambda", "cl"]:
            c_R, c_G, c_B, c_A = choose_cmap_color(
                (float(lam), float(cc)), "conc_lambda"
            )
        elif color == "all":
            c_R, c_G, c_B, c_A = choose_cmap_color((T, lam, float(cc)), "all")

        if d != "NODRG":
            if "conc" in color:
                if conc_units == "mM":
                    label = f"{T}K λ={lam}\n{cc}mM" + r" $\sigma$" + f"={sigma}"
                else:
                    drg_L_str = f"{drg_L:.2e}".replace("+", "")
                    label = (
                        f"{T}K λ={lam}\n{drg_L_str} "
                        + r"$nm^3/L$"
                        + r" $\sigma$"
                        + f"={sigma}"
                    )
            else:
                label = f"{T}K\nλ={lam}" + r" $\sigma$" + f"={sigma}"
        else:
            label = f"{T}K" + r" $\sigma$" + f"={sigma}"

        print(arr.shape)

        ax1.plot(
            arr[:, 0],
            arr[:, 1],
            label=label,
            color=(c_R, c_G, c_B),
        )

        if d != "NODRG":
            drg_path = f_path[0].replace("prot", "drg")
            arr_drg = np.loadtxt(folder + "/" + drg_path)

            if conc_units == "mM":
                label_drg = f"{d} {T}K\nλ={lam} {cc}mM" + r" $\sigma$" + f"={sigma}"
            else:
                drg_L_str = f"{drg_L:.2e}".replace("+", "")
                label_drg = (
                    f"{d} {T}K λ={lam}\n{drg_L_str} "
                    + r"$nm^3/L$"
                    + r" $\sigma$"
                    + f"={sigma}"
                )

            ax2.plot(
                arr_drg[:, 0],
                arr_drg[:, 1],
                label=label_drg,
                color=(c_R, c_G, c_B),
            )

        else:
            print("Ploteando solo prot.")
            d = "NODRG"
            try:
                f_path = [
                    f for f in os.listdir(folder) if f.endswith(".out") and "prot" in f
                ]
                arr = np.loadtxt(folder + "/" + f_path[0])

                ax1.plot(
                    arr[:, 0],
                    arr[:, 1],
                    label=f"{T}K",
                    linestyle="dotted",
                )

            except KeyError:
                print("prot ommitting pr, d, c, T, lam:", pr, d, c, T, lam)
                # print(traceback.format_exc())
                pass

    ax1.legend(loc="upper right", fontsize=10.5, handlelength=1, labelspacing=0.7)
    ax2.legend(
        loc="center right", fontsize=10.5, ncol=1, handlelength=1, labelspacing=0.7
    )
    ax1.set_xlabel("z")
    ax1.set_ylabel("Average Density")
    ax2.set_xlabel("z")
    ax2.set_ylabel("Average Density")

    plt.tight_layout()

    temps_name = (
        str(T).replace("[", "").replace("]", "").replace(",", "").replace(" ", "_")
    )

    lambdas_name = (
        str(lam)
        .replace("[", "")
        .replace("]", "")
        .replace(",", "")
        .replace(" ", "_")
        .replace("''", "NODRG")
    )

    plt.savefig(f"{p}_{lam}_{temps_name}_plot.png", dpi=500)
    plt.show()


if __name__ == "__main__":

    sim_fold = "/home/tetsuo420/Documentos/MASTER/MMCAMFQB/TFM/simulations/"

    file_paths = [
        # f"{sim_fold}sim_modifier_3/Q5-8_20_20mM_sigma02/340/",
        f"{sim_fold}sim_modifier_3/Q5-8_20_20mM_sigma02/323/",
        f"{sim_fold}sim_modifier_3/Q5-8_20_20mM_sigma2/323/",
        # f"{sim_fold}sim_modifier_sigma_mass_2/Q5-8_20_02mM/300/",
        f"{sim_fold}sim_modifier_sigma_mass_2/Q5-8_20_02mM/323/",
        f"{sim_fold}sim_GLY_lambda/Q5-8_20_GLY_lmbda_0/323/",
    ]

    plot_densities(
        file_paths,
        color="conc_sigma",
        conc_units="mM",
    )
