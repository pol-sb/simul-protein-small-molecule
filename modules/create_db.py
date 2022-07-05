import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import utils as ut
import pprint as pp


def prepare_dict_list():
    idp_paths = []
    drg_paths = []

    # Gathering the '.out' file paths
    for tup in os.walk("."):
        for fname in tup[2]:
            if (".out" in fname) and ("_drg_" not in fname):
                idp_paths.append(tup[0] + "/" + fname)
            elif (".out" in fname) and ("_drg_" in fname):
                drg_paths.append(tup[0] + "/" + fname)


    # Filling the dict_list with every record available
    dict_list = []
    for path, path2 in zip(idp_paths, drg_paths):

        fold_n = "/".join(path.split("/")[:-1]) + "/"
        params = ut.read_parameters(fold_n)
        sim_dict = {
            "protein": params["PROT_NAME"],
            "small_molec": params["DRG_NAME"],
            "conc": float(params["DRG_CONC_(mM)"]),
            "lambda": float(params["DRG_LAMB"].translate(str.maketrans("", "", "[]"))),
            "sigma": float(params["DRG_SIGMA"]),
            "temp": float(params["TEMP_(K)"]),
            "idp_average": np.loadtxt(path),
            "drg_average": np.loadtxt(path2),
        }

        dict_list.append(sim_dict)

    return dict_list


# Limiting the np.array decimal length
np.set_printoptions(threshold=2)

# Calling the function to fill the dictionary
rec_list = prepare_dict_list()


data_df = pd.json_normalize(rec_list)


# Defining new lists to be used as a new column for the plateau average.
plat_idp_list = []
plat_drg_list = []

# Adding the plateau averages to the dataframe.
for avg in data_df["idp_average"]:
    plat_idp_list.append(np.average(avg[:, 1][65:84]))

for avg in data_df["drg_average"]:
    plat_drg_list.append(np.average(avg[:, 1][65:84]))

data_df["idp_plat_avg"] = plat_idp_list
data_df["drg_plat_avg"] = plat_drg_list

print(data_df)

save = input("Save the database? (y/n) ")

if save.lower() == "y":
    data_df.to_pickle("simulations_df.pkl")
    data_df.to_csv("simulations_df.csv")
