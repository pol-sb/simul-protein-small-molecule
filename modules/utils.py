import logging
import argparse
import pprint as pp
import requests


def arg_parse():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "OpenMM script to run intrinsically disordered protein MD simulations"
            " including small molecules in the system.\nAn example"
            " command:\n\n\tsimulate.py --name Q5-8_20 --temp 323 --small_molec GLY 20"
            " 0 --lmb 0 --sigma 0.45 -m 57.5 --time 43200 --cc"
            " 10\n\n────────────────────────── Arguments ──────────────────────────"
        ),
    )

    g1 = parser.add_argument_group("Simulation parameters")

    g1.add_argument(
        "--name",
        "-N",
        nargs=1,
        type=str,
        help="Name of the protein sequence to be simulated.",
        required=True,
    )

    g1.add_argument(
        "--temp",
        "-T",
        nargs=1,
        type=int,
        help="Temperature (in K) of the system.",
        required=True,
    )

    g3 = parser.add_argument_group("Small molecule parameters")

    g3.add_argument(
        "--small_molec",
        "-s",
        nargs=3,
        default=None,
        const=None,
        type=str,
        help=(
            "Names (3 letter name), concentration (in mM) and distance between "
            "particles (in A) "
            "of the small molecules to be added. \nFor example: "
            "ARG-LYS 0.005 0.5 "
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
        help=(
            "If present, enables collision check after adding the small"
            " molecules into the system. Omit this option to disable the check."
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

    g2 = parser.add_argument_group("Simulation time selection")
    g2_1 = g2.add_mutually_exclusive_group(required=True)

    g2_1.add_argument(
        "--nsteps",
        "--steps",
        "-n",
        nargs="?",
        const=0,
        type=int,
        help="Number of timesteps to run the simulation.",
    )

    g2_1.add_argument(
        "--time",
        "--tsec",
        "-t",
        nargs="?",
        const=0,
        type=int,
        help="Number of seconds to run the simulation.",
    )

    g_sim = parser.add_argument_group("Simulation configuration")

    g_sim.add_argument(
        "--cpu",
        help="Use only the CPU as platform for the openmm.simulation.",
        action="store_true",
    )

    args = parser.parse_args()

    return args


def custom_logger(args):

    # Format strings for the logger
    fmt_str = "\n%(asctime)s - %(levelname)s:\n\t %(message)s"
    datefmt = "[%d-%m-%Y %H:%M:%S]"

    # Logger configuration
    if args.verbose and not args.quiet:
        verb = logging.DEBUG
        logging.basicConfig(
            level=verb,
            format=fmt_str,
            datefmt=datefmt,
        )
    elif args.quiet and not args.verbose:
        fmt_str = "\n%(asctime)s - \N{ESC}[31m%(levelname)s:\u001b[0m\n\t %(message)s"
        verb = logging.ERROR
        logging.basicConfig(
            level=verb,
            format=fmt_str,
            datefmt=datefmt,
        )
    else:
        verb = logging.INFO
        fmt_str = "%(message)s"
        logging.basicConfig(
            level=verb,
            format=fmt_str,
            datefmt=datefmt,
        )

    verbosity = [verb, fmt_str]
    # create logger
    logger = logging.getLogger("main simulation")
    verbosity = [verb, fmt_str]

    return logger, verbosity


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


def send_notif(title, body, pb_token):
    """Send notifications to pushbullet with the given API key."""

    url = 'https://api.pushbullet.com/v2/pushes'
    headers = {f'Access-Token': pb_token, 'Content-Type': 'application/json'}
    data = {"type": "note", "title": title, "body": body}
    req = requests.post(url, auth=(pb_token, ''), data=data)


if __name__ == "__main__":

    print(
        "This file is intented to be used as a module file and should only"
        " called from other scripts, not run by itself."
    )
