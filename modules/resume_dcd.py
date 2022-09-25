import openmm
import openmm.unit as unit
from openmm import XmlSerializer, app

import modules.small_molecule as smol
import extensions as ext
import utils as ut
import mdtraj as md
import os

pdb = app.pdbfile.PDBFile("sm_drg_traj.pdb")
top = pdb.getTopology()

chk_files = [fil for fil in os.listdir() if fil.endswith(".chk")]
xml_files = [fil for fil in os.listdir() if fil.endswith(".xml")]
dcd_files = [fil for fil in os.listdir() if fil.endswith(".dcd")]

params = ut.read_parameters()
temp = params["TEMP"]

# Change this
platform_str = "CPU"

if platform_str == "CPU":
    platform = openmm.Platform.getPlatformByName("CPU")
    platform_props = {}
else:
    platform = openmm.Platform.getPlatformByName("CUDA")
    platform_props = dict(CudaPrecision="mixed", DeviceIndex="0, 1")

integrator = openmm.openmm.LangevinIntegrator(
    temp * unit.kelvin, 0.01 / unit.picosecond, 0.005 * unit.picosecond
)

with open(xml_files[0], "r") as f:
    system_ser = f.read()
    system = XmlSerializer.deserialize(system_ser)

traj = md.load_dcd(dcd_files[0], top="sm_drg_traj.pdb")
posc = traj.xyz[-1]

simulation = app.simulation.Simulation(
    top,
    system,
    integrator,
    platform,
    platformProperties=platform_props,
)

simulation.context.setPositions(posc)

simulation.reporters.append(
    app.dcdreporter.DCDReporter("resumed_report.dcd",
        10,
        append=False,
    )
)

simulation.reporters.append(
    app.statedatareporter.StateDataReporter(
        f"resumed_simulation.log",
        reportInterval=10,
        potentialEnergy=True,
        temperature=True,
        step=True,
        speed=True,
        volume=True,
        elapsedTime=True,
        separator="\t",
        append=False,
    )
)

simulation.step(1000)


    