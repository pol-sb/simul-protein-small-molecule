import logging
from argparse import ArgumentParser


def arg_parse():
    parser = ArgumentParser()
    parser.add_argument(
        "--name",
        "-n",
        nargs=1,
        type=str,
        help="Name of the protein sequence to be simulated.",
    )
    parser.add_argument(
        "--temp",
        "-T",
        nargs=1,
        type=int,
        help="Temperature (in K) of the system.",
    )
    parser.add_argument(
        "--small_molec",
        "-sm",
        nargs=3,
        type=str,
        help=(
            "Name (3 letter name), concentration (in mM) and distance (in A)"
            " of the small molecules to be added. \nFor example:"
            " ARG-LYS 0.005 0.5"
        ),
    )

    parser.add_argument(
        "-v",
        "--verbose",
        help="Increase output verbosity",
        action="store_true",
    )

    parser.add_argument(
        "-q",
        "--quiet",
        help="Decrease output verbosity",
        action="store_true",
    )

    parser.add_argument(
        "--cpu",
        help="Use only the CPU as platform for the openmm.simulation.",
        action="store_true",
    )

    parser.add_argument(
        "--time",
        "-tsec",
        "-t",
        nargs="?",
        const=0,
        type=int,
        help="Number of seconds to run the simulation.",
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


if __name__ == "__main__":

    print(
        "This file is intented to be used as a module file and should only"
        " called from other scripts, not run by itself."
    )
