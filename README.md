# Protein Simulation

## Table of Contents

- [About](#about)
- [Getting Started](#getting_started)
- [Usage](#usage)
- [Contributing](../CONTRIBUTING.md)

## About <a name = "about"></a>

Scripts for the simulation of proteins with the capacity of introducing small molecules to the system and study their effects on protein separation.

## Getting Started <a name = "getting_started"></a>

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See [deployment](#deployment) for notes on how to deploy the project on a live system.

### Prerequisites

What things you need to install the software and how to install them.

```
Give examples
```

### Installing

A step by step series of examples that tell you how to get a development env running.

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo.

## Usage <a name = "usage"></a>

To run simulations, execute the `simulate.py` script with your python environment. The script has an argument parser which requires two arguments:
- `name`: name of the protein (e.g. `FUS`).
- `temp`: temperature in Kelvin.

An example of a command to run a simulation:

```
python3 simulate.py --name Q5-8_20 --temp 310 --small-molec GLY 20 0 --time 43200
```
