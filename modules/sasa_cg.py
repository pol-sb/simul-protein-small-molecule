##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Robert McGibbon
# Contributors:
#
# MDTraj is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 2.1
# of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with MDTraj. If not, see <http://www.gnu.org/licenses/>.
##############################################################################
#
# Modification by Ramon Crehuet to work with coarse-grained models.
# Atom radii are replaced by residue radii
#
##############################################################################

##############################################################################
# Imports
##############################################################################

from __future__ import print_function, division
from copy import deepcopy
import numpy as np
from mdtraj.utils import ensure_type
from mdtraj.geometry import _geometry

__all__ = ['shrake_rupley']

# These van der waals radii are taken from
# https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
# which references 
# A. Bondi (1964). "van der Waals Volumes and Radii". J. Phys. Chem. 68: 441. doi:10.1021/j100785a001 and doi:10.1021/jp8111556. 
# M. Mantina et al. (2009). "Consistent van der Waals Radii for the Whole Main Group". J. Phys. Chem. A. 113 (19): 5806--12. doi:10.1021/jp8111556
# Where no van der Waals value is known, a default of 2 angstroms is used.
# However, because certain atoms in biophysical simulations have a high 
# chance of being completely ionized, we have decided to give the 
# following atoms their ionic radii their ionic radii:
# +2: Be, Mg, Ca, Ba
# +1: Li, Na, K, Cs
# -1: Cl
# These ionic radii are were taken from:
# Shannon, R. D. Revised effective ionic radii and systematic studies of interatomic distances in halides and chalcogenides. Acta Crystallographica Section A 32, 751--767 (1976). doi:10.1107/S0567739476001551
# For most atoms, adding electrons usually doesn't change the radius much 
# (<10%), while removing them changes it substantially (>50%). Further, 
# when atoms like N, S, and P, are positive, they are bound to atoms in such 
# a way that would "hide" their radii anyway. We have therefore chosen to just 
# use their vdW radii.


##############################################################################
# Functions
##############################################################################


def shrake_rupley(traj, change_radii, probe_radius=0.14, n_sphere_points=960, get_mapping=False):
    """Compute the solvent accessible surface area of each atom or residue in each simulation frame.

    Parameters
    ----------
    traj : Trajectory
        An mtraj trajectory.
    probe_radius : float, optional
        The radius of the probe, in nm.
    n_sphere_points : int, optional
        The number of points representing the surface of each atom, higher
        values leads to more accuracy.
    change_radii : dict, optional
        A partial or complete dict containing the radii to change from the 
        defaults. Should take the form {"Symbol" : radii_in_nm }, e.g. 
        {"Cl" : 0.175 } to change the radii of Chlorine to 175 pm.
    get_mapping : bool, optional
        Instead of returning only the areas, also return the indices of the
        atoms or the residue-to-atom mapping. If True, will return a tuple
        that contains the areas and the mapping (np.array, shape=(n_atoms)).

    Returns
    -------
    areas : np.array, shape=(n_frames, n_features)
        The accessible surface area of each atom or residue in every frame.
        If mode == 'atom', the second dimension will index the atoms in
        the trajectory, whereas if mode == 'residue', the second
        dimension will index the residues.

    Notes
    -----
    This code implements the Shrake and Rupley algorithm, with the Golden
    Section Spiral algorithm to generate the sphere points. The basic idea
    is to great a mesh of points representing the surface of each atom
    (at a distance of the van der waals radius plus the probe
    radius from the nuclei), and then count the number of such mesh points
    that are on the molecular surface -- i.e. not within the radius of another
    atom. Assuming that the points are evenly distributed, the number of points
    is directly proportional to the accessible surface area (its just 4*pi*r^2
    time the fraction of the points that are accessible).

    There are a number of different ways to generate the points on the sphere --
    possibly the best way would be to do a little "molecular dyanmics" : put the
    points on the sphere, and then run MD where all the points repel one another
    and wait for them to get to an energy minimum. But that sounds expensive.

    This code uses the golden section spiral algorithm
    (picture at http://xsisupport.com/2012/02/25/evenly-distributing-points-on-a-sphere-with-the-golden-sectionspiral/)
    where you make this spiral that traces out the unit sphere and then put points
    down equidistant along the spiral. It's cheap, but not perfect.

    The gromacs utility g_sas uses a slightly different algorithm for generating
    points on the sphere, which is based on an icosahedral tesselation.
    roughly, the icosahedral tesselation works something like this
    http://www.ziyan.info/2008/11/sphere-tessellation-using-icosahedron.html

    References
    ----------
    .. [1] Shrake, A; Rupley, JA. (1973) J Mol Biol 79 (2): 351--71.
    """

    xyz = ensure_type(traj.xyz, dtype=np.float32, ndim=3, name='traj.xyz', shape=(None, None, 3), warn_on_cast=False)
    dim1 = xyz.shape[1]
    atom_mapping = np.arange(dim1, dtype=np.int32)

    out = np.zeros((xyz.shape[0], dim1), dtype=np.float32)
    resi_radii = [change_radii[resi.name] for resi in traj.topology.residues()]
    radii = np.array(resi_radii, np.float32) + probe_radius

    _geometry._sasa(xyz, radii, int(n_sphere_points), atom_mapping, out)

    if get_mapping == True:
        return out, atom_mapping
    else:
        return out

