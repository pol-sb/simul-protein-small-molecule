import numpy as np
import pandas as pd
import os
import utils as ut


def prepare_dict():

    paths = []
    for tup in os.walk("."):
        for fname in tup[2]:
            if (".out" in fname) and ("_drg_" not in fname):
                paths.append(tup[0] + "/" + fname)

    sim_dict = {}

    # this part fills up the dictionary with all of the paths, it is probably
    # the least elegant solution in the world, but works for now.
    for path in paths:

        fold_n = "/".join(path.split("/")[:-1]) + "/"
        params = ut.read_parameters(fold_n)
        # print("params: ", params)
        prot_name = params["PROT_NAME"]

        try:
            sim_dict[prot_name]
        except KeyError:
            sim_dict[prot_name] = {}

        drg_name = params["DRG_NAME"]
        if drg_name == "None":
            drg_name = "NODRG"

        try:
            sim_dict[prot_name][drg_name]
        except KeyError:
            sim_dict[prot_name][drg_name] = {}

        drg_conc = params["DRG_CONC_(mM)"]

        try:
            sim_dict[prot_name][drg_name][drg_conc]
        except KeyError:
            sim_dict[prot_name][drg_name][drg_conc] = {}

        lambd = params["DRG_LAMB"].replace("[", "").replace("]", "")

        try:
            lambd = float(lambd)
        except ValueError:
            pass

        try:
            sim_dict[prot_name][drg_name][drg_conc][lambd]
        except KeyError:
            sim_dict[prot_name][drg_name][drg_conc][lambd] = {}

        # print('path: ', path)
        sigma = params["DRG_SIGMA"]
        try:
            sigma = float(sigma)
        except ValueError:
            pass

        try:
            sim_dict[prot_name][drg_name][drg_conc][lambd][sigma]
        except KeyError:
            sim_dict[prot_name][drg_name][drg_conc][lambd][sigma] = {}

        temp = params["TEMP_(K)"]

        try:
            sim_dict[prot_name][drg_name][drg_conc][lambd][sigma][int(temp)]
        except KeyError:
            sim_dict[prot_name][drg_name][drg_conc][lambd][sigma][
                int(temp)
            ] = np.loadtxt(path)

    return sim_dict


np.set_printoptions(threshold=2)
data_dict = prepare_dict()
# print('data_dict: ', data_dict["Q5-8_20"]["GLY"])
# print('data_dict: ', data_dict)
data_df = pd.DataFrame(
    [
        (k1, k2, k3, k4, k5, k6, v)
        for k1, k23456v in data_dict.items()
        for k2, k3456v in k23456v.items()
        for k3, k456v in k3456v.items()
        for k4, k56v in k456v.items()
        for k5, k6v in k56v.items()
        for k6, v in k6v.items()
    ],
    columns=["protein", "small_molec", "conc", "lambd", "sigma", "temp", "average"],
)

save = input("Save the database? (y/n) ")

if save.lower() == "y":
    data_df.to_pickle("simulations_df.pkl")
    data_df.to_csv("simulations_df.csv")

print(data_df)