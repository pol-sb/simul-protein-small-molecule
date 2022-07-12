import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_pickle("simulations_df.pkl")

def choose_cmap_color(value, prop, PROTEINS, DRUGS, CONC, TEMPS, LAMBDAS):

    # TEMPS = [290, 300, 310, 315, 320, 323, 325, 330]
    # LAMBDAS = ["None", 0.0, 0.175, 0.35, 0.525, 0.7]
    # CONC = [0.2, 2.0, 10, 20.0]
    print("conc:", CONC)

    if prop in ["temperature", "T", "temp", "t"]:
        # origial -> magma
        col_arr = plt.cycler("color", plt.cm.viridis(np.linspace(0, 1, len(TEMPS))))
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

    elif prop in ["conc_lambda", "lambda_conc"]:
        combin = list(itertools.product(LAMBDAS, CONC))
        print("combin: ", combin)
        col_arr = plt.cycler("color", plt.cm.spring(np.linspace(0, 1, len(combin))))
        col_val = np.array(col_arr)[combin.index(value)]

    elif prop in ["all", "conc_lambda_temp"]:
        combin = list(itertools.product(TEMPS, LAMBDAS, CONC))
        col_arr = plt.cycler("color", plt.cm.prism(np.linspace(0, 1, len(combin))))
        col_val = np.array(col_arr)[combin.index(value)]

    return col_val["color"]



def plot_densities(
    sim_dict,
    proteins: list,
    drugs: list,
    conc: list,
    temps: list,
    lambdas: list,
    color: str,
    conc_units: str,
):

    if drugs[0] == "NODRG":
        for pr, T in itertools.product(proteins, temps):
            try:
                print("\nsim_dict: ", sim_dict[pr]["NODRG"]["0"][str(T)]["None"])
                f_path = sim_dict[pr]["NODRG"]["0"][str(T)]["None"]["path"]
                print("f_path: ", f_path)
                arr = np.loadtxt(f_path)

                c_R, c_G, c_B, c_A = choose_cmap_color(
                    T, "temp", proteins, drugs, conc, temps, lambdas
                )

                plt.plot(
                    arr[:, 0], arr[:, 1], label=f"{T}K", color=(c_R, c_G, c_B), lw=3
                )

            except KeyError:
                print(pr, T, "skipped")
                pass

        plt.legend(framealpha=0, fontsize=19)
        plt.xlabel("z", fontsize=25)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        # plt.yticks([0,350,700], [0,350,700])
        # plt.xticks([0,75,150], [0,75,150])
        plt.yticks([], [])
        plt.xticks([], [])
        plt.ylabel("Avg. Density", fontsize=25)

        # plt.title(f"{proteins[0]}, NO DRG")

        temps_name = (
            str(temps)
            .replace("[", "")
            .replace("]", "")
            .replace(",", "")
            .replace(" ", "_")
        )
        plt.tight_layout()
        plt.savefig(
            f"{proteins[0]}_NODRG_{temps_name}_plot.png", dpi=500, transparent=True
        )
        plt.show()

    else:

        fig, (ax1, ax2) = plt.subplots(ncols=2)
        fig.set_size_inches(12, 6)

        pprint.pprint(sim_dict)
        for pr, d, c, T, lam in itertools.product(
            proteins, drugs, conc, temps, lambdas
        ):
            print("\ntrying:", pr, d, c, T, lam)
            # print("lam: ", lam, "type:", type(lam))
            # print(lam == "")

            if lam not in ["", "None"] and c != 0:
                print("entra primer if")

                if lam == 0:
                    # print("lam", lam, "es 0")
                    lam = 0.0

                try:
                    # Debugging prints
                    # print("Dict_1")
                    # print(sim_dict[pr])
                    # print("Dict_2")
                    # print('d: ', d)
                    # print(sim_dict[pr][d])
                    # print("Dict_3")
                    # print(sim_dict[pr][d][str(c)])
                    print("Dict_4")
                    print(sim_dict[pr][d][str(c)][str(T)])
                    print("Fin Dicts")

                    if lam == 0.0:
                        lam_mod = str(lam)
                    else:
                        lam_mod = f"{lam:.3f}"
                        if lam_mod == "0.700":
                            lam_mod = "0.7"
                        if lam_mod == "0.900":
                            lam_mod = "0.9"
                        # if lam_mod == '0.350':
                        #     lam_mod = '0.35'

                    if lam == "":
                        lam_mod = "None"

                    print("lam_mod: ", lam_mod)
                    # print(sim_dict[pr][d][str(c)][str(T)])

                    f_path = sim_dict[pr][d][str(c)][str(T)][lam_mod]["path"]
                    # print("f_path: ", f_path)
                    fold_n = "/".join(f_path.split("/")[:-1]) + "/"

                    params = ut.read_parameters(fold_n)
                    cc = params["DRG_CONC_(mM)"]

                    if "DRG_SIGMA" in params.keys():
                        sigma = float(params["DRG_SIGMA"])
                    else:
                        sigma = float(
                            input("Please input sigma for the current simulation:")
                        )

                    drg_volume = ((4 / 3) * np.pi * ((sigma / 2) ** 3)) * int(
                        params["DRG_NUMB"]
                    )

                    # In nanometers
                    sim_box_vol = 15 * 15 * 150

                    # Converting the simbox volume to liters
                    sim_box_vol_L = sim_box_vol * 1e-24

                    drg_L = drg_volume / sim_box_vol_L

                    arr = np.loadtxt(f_path)

                    if color in ["lambda", "lambdas", "lam"]:
                        c_R, c_G, c_B, c_A = choose_cmap_color(
                            lam, "lambda", proteins, drugs, conc, temps, lambdas
                        )
                    elif color in ["temp", "temperature", "T"]:
                        c_R, c_G, c_B, c_A = choose_cmap_color(
                            T, "temp", proteins, drugs, conc, temps, lambdas
                        )
                    elif color in ["temp_lambda", "both", "lambda_temp"]:
                        c_R, c_G, c_B, c_A = choose_cmap_color(
                            (T, lam), "both", proteins, drugs, conc, temps, lambdas
                        )
                    elif color in ["conc", "c"]:
                        c_R, c_G, c_B, c_A = choose_cmap_color(
                            float(cc), "conc", proteins, drugs, conc, temps, lambdas
                        )
                    elif color in ["conc_temp", "ct"]:
                        c_R, c_G, c_B, c_A = choose_cmap_color(
                            (T, float(cc)),
                            "conc_temp",
                            proteins,
                            drugs,
                            conc,
                            temps,
                            lambdas,
                        )
                    elif color in ["conc_lambda", "cl"]:
                        c_R, c_G, c_B, c_A = choose_cmap_color(
                            (lam, float(cc)),
                            "conc_lambda",
                            proteins,
                            drugs,
                            conc,
                            temps,
                            lambdas,
                        )
                    elif color in ["all", "conc_lambda_temp", "conc_temp_lambda"]:
                        c_R, c_G, c_B, c_A = choose_cmap_color(
                            (T, lam, float(cc)),
                            "all",
                            proteins,
                            drugs,
                            conc,
                            temps,
                            lambdas,
                        )

                    if d != "NODRG":
                        if "conc" in color:
                            if conc_units == "mM":
                                label = f"{T}K λ={lam} {cc}mM"
                            else:
                                drg_L_str = f"{drg_L:.2e}".replace("+", "")
                                label = f"{T}K λ={lam}\n{drg_L_str} " + r"$nm^3/L$"
                        else:
                            label = f"{T}K\nλ={lam}"
                    else:
                        label = f"{T}K"

                    print(arr.shape)

                    ax1.plot(
                        arr[:, 0],
                        arr[:, 1],
                        label=label,
                        color=(c_R, c_G, c_B),
                    )

                    if d != "NODRG":
                        drg_path = f_path.replace("prot", "drg")
                        arr_drg = np.loadtxt(drg_path)

                        if conc_units == "mM":
                            label_drg = f"{d} {T}K λ={lam} {cc}mM"
                        else:
                            drg_L_str = f"{drg_L:.2e}".replace("+", "")
                            label_drg = f"{d} {T}K λ={lam}\n{drg_L_str} " + r"$nm^3/L$"

                        ax2.plot(
                            arr_drg[:, 0],
                            arr_drg[:, 1],
                            label=label_drg,
                            color=(c_R, c_G, c_B),
                        )

                except KeyError:
                    print(f"Skipping {pr} {d} {c} mM, {T} K, lam={lam}")
                    # print(traceback.format_exc())
                    pass
            else:
                print("Ploteando solo prot.")
                d = "NODRG"
                try:
                    f_path = sim_dict[pr][d][str(c)][str(T)][str(lam)]["path"]
                    arr = np.loadtxt(f_path)
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
            str(temps)
            .replace("[", "")
            .replace("]", "")
            .replace(",", "")
            .replace(" ", "_")
        )

        lambdas_name = (
            str(lambdas)
            .replace("[", "")
            .replace("]", "")
            .replace(",", "")
            .replace(" ", "_")
            .replace("''", "NODRG")
        )

        plt.savefig(f"{proteins[0]}_{lambdas_name}_{temps_name}_plot.png", dpi=500)
        plt.show()
