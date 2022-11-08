import argparse
import logging
import os
import subprocess as sb
import sys
import time

import numpy as np
import pandas as pd
import requests


LOGO = (
    " _     _                 _                 _ \n"
    "(_) __| |_ __        ___(_)_ __ ___  _   _| |\n"
    "| |/ _` | '_ \ _____/ __| | '_ ` _ \| | | | |\n"
    "| | (_| | |_) |_____\__ \ | | | | | | |_| | |\n"
    "|_|\__,_| .__/      |___/_|_| |_| |_|\__,_|_|\n"
    "        |_|                                  \n"
)


def check_version(path):
    """This function checks the program's github repository for updates.

    Args:
        path (string): main path of the repository containing the program.
    """
    # Saving original path
    orig_path = os.getcwd()

    # Checking the current repository commit
    try:
        os.chdir(path)

        sb.call(["git", "fetch"], stdout=sb.DEVNULL, stderr=sb.STDOUT)

        local_version = (
            sb.check_output(["git", "describe", "--abbrev=40", "--always", "--long"])
            .decode("utf-8")
            .strip()
        )
        remote_version = (
            sb.check_output(
                ["git", "describe", "--abbrev=40", "--always", "--long", "origin/main"]
            )
            .decode("utf-8")
            .strip()
        )

        print(LOGO)

        if local_version == remote_version:
            print(f"Version {local_version}.")
            print(
                "\x1b[1;32;49m"
                + "[OK]"
                + "\x1b[0m"
                + " Current script version is up-to-date!"
            )
            sys.exit(69)

        else:
            print(f"Version {local_version}.")
            print("\x1b[1;33;49m" + "[!]" + "\x1b[0m" + " Updates are available!")

    # If the update check fails this message is printed.
    except Exception:
        print(LOGO)
        print("\x1b[2;31;49m" + "[!!!] Could not check version!" + "\x1b[0m")

    # Returning to the original path.
    os.chdir(orig_path)


def find_checkpoint():

    chk_path = ""
    chk_file = ""

    for root, dirs, files in os.walk("."):
        chk_file = [f for f in files if f.endswith("restart.chk")]
        if len(chk_file) > 0:
            chk_path = root
            chk_file = chk_file[0]
            break

    # print("chk_path: ", chk_path)
    # print("chk_file: ", chk_file)
    # quit()
    return chk_path, chk_file


def add_simulargs_to_subparser(subparser):
    """This function adds simulation-related arguments to the subparser

    Args:
        subparser (argparse.ArgumentParser): 'simulation' subparser that will contain the simulation arguments.

    Returns:
        argparse.ArgumentParser: 'simulation' subparser containing simulation arguments.
    """
    subp1 = subparser.add_parser(
        "simulate",
        aliases=["sim"],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="This command allows to run IDP simulations including optional small molecules.",
        description=(
            LOGO
            + "This command allows to run IDP simulations including optional small molecules.\n"
            "See below for some usage examples."
        ),
        epilog=(
            "\n\nUsage:\n\n - An example"
            " command to run a FUS protein simulation at 323K with a 20mM concentration"
            " of a GLY-like small molecule, using 0 for its interaction strength, and "
            "defining a size (sigma) and mass (mass) for the small molecule. The"
            " simulation will run for 10000000 steps and using the gpu with ID 0:"
            "\n\n\tidp-simul simulate --name FUS --temp 323 --small_molec GLY 20"
            " 0 --lmb 0 --sigma 0.45 -m 57.5 --nsteps 10000000 --cc 10 --gpu 0"
            "\n\n - An example command to run a Q5-8_20 simulation at 315K using a "
            "small molecule chain with 5 residues, while defining the interaction "
            "strength of each residue. The simulation will run for 43200 seconds and"
            " using the cpu:"
            "\n\n\tidp-simul simulate --name Q5-8_20 --temp 315 --small_molec GLY-GLY-GLY-ARG-"
            "GLY 1 0.5 --lmb 0.350 0.900 0.1 0.1 0.1 --time 43200 --cc 10 --cpu\n"
            "When a simulation is started, a folder named after the protein name will be created"
        ),
    )

    g1 = subp1.add_argument_group(title="Simulation parameters")
    g1.add_argument(
        "--name",
        "-N",
        nargs=1,
        type=str,
        metavar="NAME",
        help="Name of the protein sequence to be simulated.",
        required=True,
    )

    g1.add_argument(
        "--temp",
        "-T",
        nargs=1,
        type=int,
        metavar="TEMP",
        help="Temperature (in K) of the system.",
        required=True,
    )

    g1.add_argument(
        "--resume",
        "--continue" "--res",
        action="store_true",
        help="If this argument is passed, the program will attempt to resume a simulation.",
    )

    g3 = subp1.add_argument_group("Small molecule parameters")

    g3.add_argument(
        "--small_molec",
        "-s",
        nargs=3,
        default=None,
        const=None,
        type=str,
        metavar=("RES", "CONC", "DISTANCE"),
        help=(
            "Residue Names (3 letter name, if more than one, joined by a hyphen),"
            " concentration (in mM) and distance between particles (in A) of the small"
            " molecules to be added. \nFor example: "
            "--small_molec ARG-LYS-GLY-GLY-GLY 20 0.5"
        ),
    )

    g3.add_argument(
        "--lambd",
        "--lam",
        "--lmb",
        nargs="+",
        type=float,
        help=(
            "List of float lambda values to use for the small molecules, given "
            "in the same order as the small molecules."
            "\nFor example: --lmbd 0.350 0.900 0.1 0.6 0.1"
        ),
    )

    g3.add_argument(
        "--sigma",
        "--sig",
        nargs="+",
        type=float,
        help=(
            "List of float sigma values to use for the small molecules, given "
            "in the same order as the small molecules."
        ),
    )

    g3.add_argument(
        "--mass",
        "-m",
        nargs="+",
        type=float,
        help=(
            "List of float molar mass values to use for the small molecules, given in"
            " the same order as the small molecules."
        ),
    )

    g3.add_argument(
        "--check_collision",
        "--cc",
        type=float,
        nargs=1,
        metavar="DISTANCE",
        help=(
            "If present, enables collision check after adding the small"
            " molecules into the system. Omit this option to disable the check."
            "Takes a distance value (in A) between the small molecules and the protein"
            " to check for collision."
        ),
    )

    g4 = subp1.add_argument_group("Output options")

    g4.add_argument(
        "-v",
        "--verbose",
        help="Increase output verbosity",
        action="store_true",
    )

    g4.add_argument(
        "-q",
        "--quiet",
        help="Decrease output verbosity",
        action="store_true",
    )

    g4.add_argument(
        "--notif",
        help="Send PB notification",
        action="store_true",
    )

    g2 = subp1.add_argument_group("Simulation time selection")
    g2_1 = g2.add_mutually_exclusive_group(required=True)

    g2_1.add_argument(
        "--nsteps",
        "--nstep",
        "--tstep",
        "--step",
        "--steps",
        "-n",
        nargs="?",
        metavar="NSTEPS",
        const=0,
        type=int,
        help="Number of timesteps to run the simulation.",
    )

    g2_1.add_argument(
        "--time",
        "--tsec",
        "-t",
        nargs="?",
        metavar="NSECONDS",
        const=0,
        type=int,
        help="Number of seconds to run the simulation.",
    )

    g_sim = subp1.add_argument_group("Simulation configuration")
    g_sim_2 = g_sim.add_mutually_exclusive_group(required=True)

    g_sim_2.add_argument(
        "--cpu",
        help="Use the CPU as platform for the openmm.simulation.",
        action="store_true",
    )
    g_sim_2.add_argument(
        "--gpu",
        nargs="+",
        metavar="GPUID",
        type=int,
        help=(
            "Use the GPU as platform for the openmm.simulation. "
            "Several GPUs can be used at once by passing their indexes separated by a space."
        ),
    )

    ex_sim = subp1.add_argument_group("Simulation post-treatment")

    ex_sim.add_argument(
        "--extend_thermostat",
        "--ethermo",
        "--et",
        type=float,
        metavar=("TEMP", "NSTEPS"),
        nargs=2,
        help=(
            "If present, after finishing the main dynamics, modify the thermostat"
            " temperature and run the simulation for a given number of steps. "
            "--extend_thermostat takes two arguments: the first is the new "
            "temperature and the second is the number of steps to run the "
            "simulation."
        ),
    )

    return subparser


def add_simulargs_to_parser(parser):
    """This function adds simulation-related arguments to the subparser

    Args:
        subparser (argparse.ArgumentParser): 'simulation' subparser that will contain the simulation arguments.

    Returns:
        argparse.ArgumentParser: 'simulation' subparser containing simulation arguments.
    """

    g1 = parser.add_argument_group(title="Simulation parameters")
    g1.add_argument(
        "--name",
        "-N",
        nargs=1,
        type=str,
        metavar="NAME",
        help="Name of the protein sequence to be simulated.",
        required=True,
    )

    g1.add_argument(
        "--temp",
        "-T",
        nargs=1,
        type=int,
        metavar="TEMP",
        help="Temperature (in K) of the system.",
        required=True,
    )

    g1.add_argument(
        "--resume",
        "--continue",
        "--res",
        action="store_true",
        help="If this argument is passed, the program will attempt to resume a simulation.",
    )

    g3 = parser.add_argument_group("Small molecule parameters")

    g3.add_argument(
        "--small_molec",
        "-s",
        nargs=3,
        default=None,
        const=None,
        type=str,
        metavar=("RES", "CONC", "DISTANCE"),
        help=(
            "Residue Names (3 letter name, if more than one, joined by a hyphen),"
            " concentration (in mM) and distance between particles (in A) of the small"
            " molecules to be added. \nFor example: "
            "--small_molec ARG-LYS-GLY-GLY-GLY 20 0.5"
        ),
    )

    g3.add_argument(
        "--lambd",
        "--lam",
        "--lmb",
        nargs="+",
        type=float,
        help=(
            "List of float lambda values to use for the small molecules, given "
            "in the same order as the small molecules."
            "\nFor example: --lmbd 0.350 0.900 0.1 0.6 0.1"
        ),
    )

    g3.add_argument(
        "--sigma",
        "--sig",
        nargs="+",
        type=float,
        help=(
            "List of float sigma values to use for the small molecules, given "
            "in the same order as the small molecules."
        ),
    )

    g3.add_argument(
        "--mass",
        "-m",
        nargs="+",
        type=float,
        help=(
            "List of float molar mass values to use for the small molecules, given in"
            " the same order as the small molecules."
        ),
    )

    g3.add_argument(
        "--check_collision",
        "--cc",
        type=float,
        nargs=1,
        metavar="DISTANCE",
        help=(
            "If present, enables collision check after adding the small"
            " molecules into the system. Omit this option to disable the check."
            "Takes a distance value (in A) between the small molecules and the protein"
            " to check for collision."
        ),
    )

    g4 = parser.add_argument_group("Output options")

    g4.add_argument(
        "-v",
        "--verbose",
        help="Increase output verbosity",
        action="store_true",
    )

    g4.add_argument(
        "-q",
        "--quiet",
        help="Decrease output verbosity",
        action="store_true",
    )

    g4.add_argument(
        "--notif",
        "--notify",
        help="Send push notification",
        action="store_true",
    )

    g2 = parser.add_argument_group("Simulation time selection")
    g2_1 = g2.add_mutually_exclusive_group(required=True)

    g2_1.add_argument(
        "--nsteps",
        "--nstep",
        "--tstep",
        "--step",
        "--steps",
        "-n",
        nargs="?",
        metavar="NSTEPS",
        const=0,
        type=int,
        help="Number of timesteps to run the simulation.",
    )

    g2_1.add_argument(
        "--time",
        "--tsec",
        "-t",
        nargs="?",
        metavar="NSECONDS",
        const=0,
        type=int,
        help="Number of seconds to run the simulation.",
    )

    g_sim = parser.add_argument_group("Simulation configuration")
    g_sim_2 = g_sim.add_mutually_exclusive_group(required=True)

    g_sim_2.add_argument(
        "--cpu",
        help="Use the CPU as platform for the openmm.simulation.",
        action="store_true",
    )
    g_sim_2.add_argument(
        "--gpu",
        nargs="+",
        metavar="GPUID",
        type=int,
        help=(
            "Use the GPU as platform for the openmm.simulation. "
            "Several GPUs can be used at once by passing their indexes separated by a space."
        ),
    )

    ex_sim = parser.add_argument_group("Simulation post-treatment")

    ex_sim.add_argument(
        "--extend_thermostat",
        "--ethermo",
        "--et",
        type=float,
        metavar=("TEMP", "NSTEPS"),
        nargs=2,
        help=(
            "If present, after finishing the main dynamics, modify the thermostat"
            " temperature and run the simulation for a given number of steps. "
            "--extend_thermostat takes two arguments: the first is the new "
            "temperature and the second is the number of steps to run the "
            "simulation."
        ),
    )

    return parser


def arg_parse():
    parser = argparse.ArgumentParser(
        prog="idp-simul",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            LOGO
            + "\nOpenMM script to run intrinsically disordered protein (IDP) molecular\n"
            "dynamics simulations with the capability of including user-defined\n"
            "small molecules in the system."
        ),
    )

    sim_par = parser.add_subparsers(
        title="Operation Modes",
        metavar="subparser",
        help="Subcommands that are used to select the operation mode of the program.",
        description="Subcommands that are used to select the operation mode of the program.",
        required=True,
        dest="subparser_name",
    )

    subp1 = add_simulargs_to_subparser(sim_par)

    subp2 = sim_par.add_parser(
        "check_version",
        help="This command allows to check the program's version.",
        description="This command allows to check the program's version.",
    )

    subp3 = sim_par.add_parser(
        "database",
        help="Commands related to working with the IDP simulations database.",
        description="Commands related to working with the IDP simulations database.",
    )

    subp4 = sim_par.add_parser(
        "tests",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="Run several tests in order to check several IDP properties.",
        description="Run several tests in order to check several IDP properties.\n"
        "Available tests shown below.",
    )

    subp4_tests = subp4.add_subparsers(
        title="Test types",
        description="Subcommands that are used to select which tests are run.",
        required=True,
        dest="test_name",
    )

    subp_montecarlo = subp4_tests.add_parser(
        "minimize",
        # parents=[sim_par],
        aliases=["dynmin", "montecarlo", "dynamic-minimize"],
        help="This command allows to iteratively modify a reaction coordinate of an IDP system.",
        description="Run short molecular dynamics simulations in a certain IDP system containing"
        " small molecules and accept the trajectory according to a reaction coordinate modification. "
        "Available parameters shown below:",
    )

    # Creating group for the minimization parameter
    minim_group = subp_montecarlo.add_argument_group(title="Minimization parameters")

    # Mutually exclusive group so only one of the parameters is allowed at once.
    minim_excl_group = minim_group.add_mutually_exclusive_group(required=True)

    minim_excl_group.add_argument(
        "--potential",
        "--potenerg",
        "--Uenerg",
        action="store_true",
        help="If this argument is passed, the program will minimize the potential energy.",
    )

    minim_excl_group.add_argument(
        "--dispersion",
        "--disp",
        "--particle_dispersion",
        "--particle_disp",
        "--p_disp",
        action="store_true",
        help="If this argument is passed, the program will minimize the particle dispersion.",
    )

    subp5 = add_simulargs_to_parser(subp_montecarlo)

    g0 = subp_montecarlo.add_argument_group(title="Minimization parameters")
    g0.add_argument(
        "--runs",
        "--mcrun",
        "--ntimes",
        nargs=1,
        type=int,
        metavar="NAME",
        help="Number of minimization runs to be performed.",
        required=True,
    )

    subp_contnodrg = subp4_tests.add_parser(
        "cont-nodrg",
        aliases=["nodrg", "contnodrg", "continue-nodrg"],
        help="This command allows to resume a simulation removing the small molecules added previously.",
        description="Continue molecular dynamics simulations for IDP systems using previously generated"
        " checkpoints, removing any small molecules present. ",
    )

    args = parser.parse_args()

    return args


def create_dirs(args):
    """Generate the folders required for a given type of simulation.

    Args:
        args (argparse.Namespace): argparse.Namespace containing all of the launch arguments for the program.
    """

    folder_date = time.strftime("%Y%m%d_%H%M%S")

    if args.subparser_name == "simulate":
        prefix = ""
    elif args.subparser_name == "tests":
        if args.test_name == "montecarlo":
            prefix = "mc_minimize_"
        else:
            prefix = "test_"

    try:
        folder_path = f"./{prefix}{args.name[0]}_{folder_date}"
        os.mkdir(folder_path)
    except FileExistsError:
        pass

    try:
        folder_path = f"./{prefix}{args.name[0]}_{folder_date}/{int(args.temp[0])}"
        os.mkdir(folder_path)
    except FileExistsError:
        pass

    return folder_path, prefix


def custom_logger(args):
    # Attempting to create directories in which to save the topology
    folder_path, prefix = create_dirs(args)

    # fmt_str = "\n%(asctime)s - %(levelname)s:\n\t %(message)s"
    # logname = f"./{args.name[0]}/{int(args.temp[0])}/idp-simul_logger.log"
    logname = f"{folder_path}/idp-simul_logger.log"

    # Format strings for the debug logger
    fmt_debug = "%(asctime)s | %(levelname)s | %(message)s"
    datefmt_debug = "%d-%m-%Y %H:%M:%S"

    # Format strings for the regular logger
    fmt_info = "%(message)s"

    # Format strings for the error logger
    fmt_error = "\n%(asctime)s - \N{ESC}[31m%(levelname)s:\u001b[0m\n\t %(message)s"

    # Logger configuration
    if args.verbose and not args.quiet:
        verb = logging.DEBUG
        logging.basicConfig(
            filename=logname,
            filemode="a",
            level=verb,
            format=fmt_debug,
            datefmt=datefmt_debug,
        )

        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.DEBUG)
        out_fmt_debug = logging.Formatter(fmt_debug)
        stdout_handler.setFormatter(out_fmt_debug)

        fmt_string = fmt_debug

    elif args.quiet and not args.verbose:
        verb = logging.DEBUG
        logging.basicConfig(
            filename=logname,
            filemode="a",
            level=verb,
            format=fmt_debug,
            datefmt=datefmt_debug,
        )

        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.ERROR)
        out_fmt_error = logging.Formatter(fmt_error)
        stdout_handler.setFormatter(out_fmt_error)

        fmt_string = fmt_error

    else:
        verb = logging.INFO
        logging.basicConfig(
            filename=logname,
            filemode="a",
            level=verb,
            format=fmt_debug,
            datefmt=datefmt_debug,
        )

        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.DEBUG)
        out_fmt_info = logging.Formatter(fmt_info)
        stdout_handler.setFormatter(out_fmt_info)

        fmt_string = fmt_info

    # create logger
    logger = logging.getLogger("idp-simul logger")
    logger.addHandler(stdout_handler)
    verbosity = [verb, fmt_string]

    return logger, verbosity, folder_path


def read_parameters(path=""):
    # Opening the parameter file
    with open(f"{path}parameters.dat", "r") as f:
        params_raw = f.readlines()

    # Initializing parameter dictionary
    param_dict = {}

    # Reading the parameters and storing them into a dictionary iteratively
    for param in params_raw[3:]:
        param_clean = param.strip().split("\t")
        param_dict[param_clean[0]] = param_clean[1]

    param_dict["START_DATE"] = params_raw[1].replace("#", "").split(" - ")[0].strip()
    param_dict["START_TIME"] = params_raw[1].replace("#", "").split(" - ")[1].strip()

    return param_dict


def create_hash_path(path=""):
    """
    Creates a unique hash for each simulation, using the file path.
    This will allow to differentiate between simulations with the same parameters but
    executed at different times/locations.

    Parameters
    ----------
    path : str, optional
        Path of a folder containing a 'parameters.dat' file, by default empty.

    Returns
    -------
    str
        Unique hash value for a given simulation.
    """

    # Getting parameters
    p_dic = read_parameters(path)

    # Preparing a tuple from the dictionary values. A tuple is needed because it is
    # inmutable, needed for the hash function to work.
    hash_tupl = tuple(p_dic.values())

    # Creating a hash from the tuple and converting it to string.
    hash_str = str(hash(hash_tupl))

    return hash_str


def create_hash(params):
    """
    Creates a unique hash for each simulation. This will allow to differentiate between
    simulations with the same parameters but executed at different times/locations.

    Parameters
    ----------
    params : list
        List containing all of the simulation parameters.

    Returns
    -------
    str
        Unique hash value for a given simulation.
    """

    # Preparing a tuple from the dictionary values. A tuple is needed because it is
    # inmutable, needed for the hash function to work.
    tupled_list = []

    for param in params:
        if type(param) == list:
            for p in param:
                if type(p) == list:
                    tupled_list.append(tuple(p))
                else:
                    tupled_list.append(p)
        else:
            tupled_list.append(param)

    hash_tupl = tuple(tupled_list)

    # Creating a hash from the tuple and converting it to string.
    hash_str = str(hash(hash_tupl))

    return hash_str


def write_params(
    path: str,
    args,
    sm_mol,
    drg_param,
    sim_time,
    time_units,
):

    param_list = [
        args.__dict__.values(),
        sm_mol,
        drg_param,
        sim_time,
        time_units,
        time.strftime("%d.%m.%Y - %H:%M:%S"),
    ]

    hash_str = create_hash(param_list)

    with open(path, "w+") as f:
        f.write("# Simulation parameters\n")
        f.write(f"# {time.strftime('%d.%m.%Y - %H:%M:%S')}\n\n")

        for key, value in args.__dict__.items():
            if key not in ["proteins", "residues", "verbosity"]:
                key = key.upper()
                value = str(value).replace("'", "").replace("[", "").replace("]", "")

                f.write(f"{key}\t{value}\n")

        f.write(f"HASH\t{hash_str}\n")


def send_notif(title, body, topic_name, folder_path):
    """
    Send push notifications to a smartphone using the open source ntfy.sh service.
    The ntfy app should be installed in the target smartphone, and a topic must be
    created in the app. This topic's name must be added to the 'modules/.ntfy_topic'
    file, which will be used to send the notification via a POST request.
    See: https://ntfy.sh/ and https://ntfy.sh/docs/subscribe/phone/ for more information.

    Args:
        title (str): Title of the notification
        body (str): Main body of the notification
        pb_token (str): Name of the ntfy topic. Used for the request.
    """

    files = [f for f in os.listdir(folder_path) if f.endswith(".png")]
    if len(files) > 0:
        img_name = files[0]
        image_path = f"{folder_path}/{img_name}"
    else:
        image_path = ""

    requests.post(
        f"https://ntfy.sh/{topic_name}",
        data=body,
        headers={
            "Title": title,
            "Tags": "white_check_mark,atom_symbol",
        },
    )

    if image_path:
        requests.post(
            f"https://ntfy.sh/{topic_name}",
            data=open(image_path, "rb"),
            headers={
                "Title": title,
                "Tags": "white_check_mark,atom_symbol",
                "Filename": img_name,
            },
        )


def timesteps_to_ns(timesteps):

    # timestep to picosecond
    tstep_res = timesteps * 0.005

    # timestep to nanosecond
    tstep_res /= 1000

    return tstep_res


def get_gpus():
    gpu_list = sb.check_output(["nvidia-smi", "-L"]).decode("utf-8").strip()
    return gpu_list


def prepare_dict_list():
    idp_paths = []
    drg_paths = []

    # Gathering the '.out' file paths
    for tup in os.walk("."):
        drg_cnt = 0
        idp_cnt = 0

        for fname in tup[2]:

            if (".out" in fname) and ("_drg_" not in fname):
                idp_paths.append(tup[0] + "/" + fname)
                idp_cnt += 1

            elif (".out" in fname) and ("_drg_" in fname):
                # print(tup[0] + "/" + fname)
                drg_paths.append(tup[0] + "/" + fname)

                drg_cnt += 1

        if idp_cnt > 0 and drg_cnt == 0:
            drg_paths.append("None")

    # Filling the dict_list with every record available
    dict_list = []
    for path, path2 in zip(idp_paths, drg_paths):

        fold_n = "/".join(path.split("/")[:-1]) + "/"
        params = read_parameters(fold_n)
        # hash_str = create_hash(fold_n)

        if params["DRG_LAMB"] in [None, "None"]:
            lamb = None
        else:
            lamb = float(params["DRG_LAMB"].translate(str.maketrans("", "", "[]")))

        if params["DRG_NAME"] == "NODRG":
            drg_name = None
        else:
            drg_name = params["DRG_NAME"]

        if params["TEMP_(K)"] in [None, "None"]:
            temp = None
        else:
            temp = float(params["TEMP_(K)"])

        if params["DRG_SIGMA"] in [None, "None"]:
            sigma = None
        else:
            sigma = float(params["DRG_SIGMA"])

        if params["DRG_CONC_(mM)"] in [None, "None"]:
            conc = None
        else:
            conc = float(params["DRG_CONC_(mM)"])

        try:
            drg_avg = np.loadtxt(path2)
        except FileNotFoundError:
            drg_avg = "None"
        except IOError:
            drg_avg = "None"

        try:
            if params["EXTENSION"] in [None, "None"]:
                extension = None
            else:
                extension = params["EXTENSION"]
        except KeyError:
            extension = None

        try:
            if params["TIME_UNIT"] in [None, "None"]:
                time_unit = None
            else:
                time_unit = params["TIME_UNIT"]
        except KeyError:
            time_unit = None

        try:
            if params["SIM_TIME"] in [None, "None"]:
                sim_time = None
            else:
                sim_time = params["SIM_TIME"]
        except KeyError:
            sim_time = None

        try:
            if params["HASH"] in [None, "None"]:
                hash_str = None
            else:
                hash_str = params["HASH"]
        except KeyError:
            hash_str = None

        sim_dict = {
            "protein": params["PROT_NAME"],
            "small_molec": drg_name,
            "conc": conc,
            "lambd": lamb,
            "sigma": sigma,
            "temp": temp,
            "idp_average": np.loadtxt(path),
            "drg_average": drg_avg,
            "hash": hash_str,
            "extension": extension,
            "sim_time": sim_time,
            "time_unit": time_unit,
        }

        dict_list.append(sim_dict)

    return dict_list


def generate_db():

    # Calling the function to fill the dictionary
    rec_list = prepare_dict_list()

    data_df = pd.json_normalize(rec_list)

    # Defining new lists to be used as a new column for the plateau average.
    plat_idp_list = []
    plat_drg_list = []
    dilute_idp_list = []
    dilute_drg_list = []

    # Adding the plateau averages to the dataframe.
    for avg in data_df["idp_average"]:
        plat_idp_list.append(np.average(avg[:, 1][65:84]))

        dilute_region = np.hstack([avg[:, 1][:65], avg[:, 1][84:]])
        dilute_idp_list.append(np.average(dilute_region))

    for avg in data_df["drg_average"]:

        if type(avg) == np.ndarray:
            plat_drg_list.append(np.average(avg[:, 1][65:84]))
            dilute_drg_region = np.hstack([avg[:, 1][:65], avg[:, 1][84:]])
            dilute_drg_list.append(np.average(dilute_drg_region))

        else:
            plat_drg_list.append(None)
            dilute_drg_list.append(None)

    data_df["idp_plat_avg"] = plat_idp_list
    data_df["drg_plat_avg"] = plat_drg_list

    data_df["idp_dilu_avg"] = dilute_idp_list
    data_df["drg_dilu_avg"] = dilute_drg_list

    return data_df


def save_db(data_df):

    save = input("Save the database? (y/n) ")

    if save.lower() == "y":
        data_df.to_pickle("simulations_df.pkl")
        data_df.to_csv("simulations_df.csv")


if __name__ == "__main__":
    print(
        "This file is intented to be used as a module file and should only"
        " called from other scripts, not run by itself."
    )
