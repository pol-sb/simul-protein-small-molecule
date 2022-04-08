import os

import MDAnalysis
import mdtraj as md
import numpy as np
import pandas as pd
from MDAnalysis import transformations


def initProteins():
    proteins = pd.DataFrame(
        index=[
            "A1star200",
            "A1star300",
            "A1star500",
            "ht40",
            "FUS",
            "M4D",
            "A1",
            "aSyn",
            "LAF1",
            "LAF1_80",
            "LAF1D2130",
            "LAF1shuf",
        ],
        columns=["eps_factor", "pH", "ionic", "fasta"],
    )
    fasta_ht40 = """MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKESPLQTPTEDGSEEP
GSETSDAKSTPTAEDVTAPLVDEGAPGKQAAAQPHTEIPEGTTAEEAGIGDTPSLEDEAAGHVTQARMVSKSKDGTGSDDKK
AKGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAPKTPPSSGEPPKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREP
KKVAVVRTPPKSPSSAKSRLQTAPVPMPDLKNVKSKIGSTENLKHQPGGGKVQIINKKLDLSNVQSKCGSKDNIKHVPGGGS
VQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTD
HGAEIVYKSPVVSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL""".replace(
        "\n", ""
    )
    fasta_FUS = """MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYG
QSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQS
SSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNS""".replace(
        "\n", ""
    )
    fasta_M4D = """GSMASASSSQRGRSGSGNFGGGRGGGFGGNGNFGRGGNFSGRGGFGGSRGGGGYGGSGGGYNGFGNSGSNFGGGGSYNGFGNYNNQ
SSNFGPMKGGNFGGRSSGPYGGGGQYFAKPRNQGGYGGSSSSSSYGSGRRF""".replace(
        "\n", ""
    )
    fasta_A1 = """GSMASASSSQRGRSGSGNFGGGRGGGFGGNDNFGRGGNFSGRGGFGGSRGGGGYGGSGDGYNGFGNDGSNFGGGGSYNDFGNYNNQ
SSNFGPMKGGNFGGRSSGGSGGGGQYFAKPRNQGGYGGSSSSSSYGSGRRF""".replace(
        "\n", ""
    )
    fasta_A1star = """GSMASASSSQRGRSGSGNFGGGRGGGFGGNDNFGRGGNFSGRGGFGGSRGGGGYGGSGDGYNGFGNDGSNFGGGGNYNNQ
SSNFGPMKGGNFGGRSSGPYGGGGQYFAKPRNQGGYGGSSSSSSYGSGRRF""".replace(
        "\n", ""
    )
    fasta_aSyn = """MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTK
EQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA""".replace(
        "\n", ""
    )
    fasta_LAF1 = """MESNQSNNGGSGNAALNRGGRYVPPHLRGGDGGAAAAASAGGDDRRGGAGGGGYRRGGGNSGGGGGGG
YDRGYNDNRDDRDNRGGSGGYGRDRNYEDRGYNGGGGGGGNRGYNNNRGGGGGGYNRQDRGDGGSSNFSRGG
YNNRDEGSDNRGSGRSYNNDRRDNGGDGLEHHHHHH""".replace(
        "\n", ""
    )
    fasta_LAF1D2130 = """MESNQSNNGGSGNAALNRGGDGGAAAAASAGGDDRRGGAGGGGYRRGGGNSGGGGGGGYDRGYNDNRDDR
DNRGGSGGYGRDRNYEDRGYNGGGGGGGNRGYNNNRGGGGGGYNRQDRGDGGSSNFSRGGYNNRDEGSDN
RGSGRSYNNDRRDNGGDGLEHHHHHH""".replace(
        "\n", ""
    )
    fasta_LAF1shuf = """MNNSGDNDRGSGNYGLRNSFGDDGYGDNGNDEGNSGYRNRGLGGDRADEYGNSGGNGDNEAAPNASDRDD
AHYYDSDDYDDGGGGRGSGGAGGGGARGPGSNRAGRYGGGGRRGRGRGNGYNGNRSQRRRGGGRGRGNRG
YRVGNGNGQSGGRNSRGGGGGNGGANYGLEHHHHHH""".replace(
        "\n", ""
    )

    proteins.loc["LAF1"] = dict(
        eps_factor=0.2, pH=7.5, fasta=list(fasta_LAF1), ionic=0.15
    )
    proteins.loc["LAF1D2130"] = dict(
        eps_factor=0.2, pH=7.5, fasta=list(fasta_LAF1D2130), ionic=0.15
    )
    proteins.loc["LAF1shuf"] = dict(
        eps_factor=0.2, pH=7.5, fasta=list(fasta_LAF1shuf), ionic=0.15
    )
    proteins.loc["LAF1_80"] = dict(
        eps_factor=0.2, pH=7.4, fasta=list(fasta_LAF1), ionic=0.08
    )
    proteins.loc["FUS"] = dict(
        eps_factor=0.2, pH=5.5, fasta=list(fasta_FUS), ionic=0.15
    )
    proteins.loc["A1"] = dict(eps_factor=0.2, pH=7.0, fasta=list(fasta_A1), ionic=0.15)
    proteins.loc["A1star200"] = dict(
        eps_factor=0.2, pH=7.0, fasta=list(fasta_A1star), ionic=0.2
    )
    proteins.loc["A1star300"] = dict(
        eps_factor=0.2, pH=7.0, fasta=list(fasta_A1star), ionic=0.3
    )
    proteins.loc["A1star500"] = dict(
        eps_factor=0.2, pH=7.0, fasta=list(fasta_A1star), ionic=0.5
    )
    proteins.loc["M4D"] = dict(
        eps_factor=0.2, pH=7.0, fasta=list(fasta_M4D), ionic=0.15
    )
    proteins.loc["ht40"] = dict(
        eps_factor=0.2, pH=7.4, fasta=list(fasta_ht40), ionic=0.005
    )
    proteins.loc["aSyn"] = dict(
        eps_factor=0.2, pH=7.4, fasta=list(fasta_aSyn), ionic=0.005
    )
    return proteins


def genParamsLJ(df, name, prot):

    # I added a [0] in here.
    fasta = prot.fasta.copy()

    r = df.copy()
    # print('r: ', r)
    r.loc["X"] = r.loc[fasta[0]]
    r.loc["Z"] = r.loc[fasta[-1]]
    r.loc["X", "MW"] += 2
    r.loc["Z", "MW"] += 16
    fasta[0] = "X"
    fasta[-1] = "Z"
    types = list(np.unique(fasta))
    MWs = [r.loc[a, "MW"] for a in types]
    lj_eps = prot.eps_factor * 4.184
    return lj_eps, fasta, types, MWs


def genParamsDH(df, name, prot, temp):
    kT = 8.3145 * temp * 1e-3
    fasta = prot.fasta.copy()
    r = df.copy()

    # Set the charge on HIS based on the pH of the protein solution
    r.loc["H", "q"] = int(1.0 / (1 + 10 ** (prot.pH - 6)))

    r.loc["X"] = r.loc[fasta[0]]
    r.loc["Z"] = r.loc[fasta[-1]]
    fasta[0] = "X"
    fasta[-1] = "Z"
    r.loc["X", "q"] = r.loc[prot.fasta[0], "q"] + 1.0
    r.loc["Z", "q"] = r.loc[prot.fasta[-1], "q"] - 1.0
    # Calculate the prefactor for the Yukawa potential
    fepsw = (
        lambda T: 5321 / T
        + 233.76
        - 0.9297 * T
        + 0.1417 * 1e-2 * T * T
        - 0.8292 * 1e-6 * T**3
    )
    epsw = fepsw(temp)
    lB = 1.6021766**2 / (4 * np.pi * 8.854188 * epsw) * 6.022 * 1000 / kT
    yukawa_eps = [r.loc[a].q * np.sqrt(lB * kT) for a in fasta]
    # Calculate the inverse of the Debye length
    yukawa_kappa = np.sqrt(8 * np.pi * lB * prot.ionic * 6.022 / 10)
    return yukawa_eps, yukawa_kappa


def genDCD(residues, name, prot, temp, n_chains):
    top = md.Topology()
    for _ in range(n_chains):
        chain = top.add_chain()
        for resname in prot.fasta:
            residue = top.add_residue(residues.loc[resname, "three"], chain)
            top.add_atom(
                residues.loc[resname, "three"],
                element=md.element.carbon,
                residue=residue,
            )
        for i in range(chain.n_atoms - 1):
            top.add_bond(chain.atom(i), chain.atom(i + 1))

    t = md.load(name + "/{:d}".format(temp) + "/{:s}.dcd".format(name), top=top)
    t.xyz *= 10
    t.unitcell_lengths *= 10
    lz = t.unitcell_lengths[0, 2]
    edges = np.arange(-lz / 2.0, lz / 2.0, 1)
    dz = (edges[1] - edges[0]) / 2.0
    z = edges[:-1] + dz
    h = np.apply_along_axis(lambda a: np.histogram(a, bins=edges)[0], 1, t.xyz[:, :, 2])
    zmid = np.apply_along_axis(lambda a: z[a.argmax()], 1, h)
    indices = np.argmin(np.abs(t.xyz[:, :, 2] - zmid[:, np.newaxis]), axis=1)
    t[0].save_pdb(name + "/{:d}".format(temp) + "/top.pdb")
    t.save_dcd(name + "/{:d}".format(temp) + "/traj3.dcd")

    u = MDAnalysis.Universe(
        name + "/{:d}".format(temp) + "/top.pdb",
        name + "/{:d}".format(temp) + "/traj3.dcd",
    )
    ag = u.atoms
    with MDAnalysis.Writer(name + "/{:d}".format(temp) + "/traj2.dcd", ag.n_atoms) as W:
        for ts, ndx in zip(u.trajectory, indices):
            ts = transformations.unwrap(ag)
            ts = transformations.center_in_box(
                u.select_atoms("index {:d}".format(ndx)), center="geometry"
            )(ts)
            ts = transformations.wrap(ag)(ts)
            W.write(ag)

    t = md.load(
        name + "/{:d}".format(temp) + "/traj2.dcd",
        top=name + "/{:d}".format(temp) + "/top.pdb",
    )
    lz = t.unitcell_lengths[0, 2]
    edges = np.arange(0, lz, 1)
    dz = (edges[1] - edges[0]) / 2.0
    z = edges[:-1] + dz
    h = np.apply_along_axis(lambda a: np.histogram(a, bins=edges)[0], 1, t.xyz[:, :, 2])
    zmid = np.apply_along_axis(lambda a: z[a > np.quantile(a, 0.98)].mean(), 1, h)
    indices = np.argmin(np.abs(t.xyz[:, :, 2] - zmid[:, np.newaxis]), axis=1)

    u = MDAnalysis.Universe(
        name + "/{:d}".format(temp) + "/top.pdb",
        name + "/{:d}".format(temp) + "/traj2.dcd",
    )
    ag = u.atoms
    with MDAnalysis.Writer(name + "/{:d}".format(temp) + "/traj1.dcd", ag.n_atoms) as W:
        for ts, ndx in zip(u.trajectory, indices):
            ts = transformations.unwrap(ag)(ts)
            ts = transformations.center_in_box(
                u.select_atoms("index {:d}".format(ndx)), center="geometry"
            )(ts)
            ts = transformations.wrap(ag)(ts)
            W.write(ag)

    t = md.load(
        name + "/{:d}".format(temp) + "/traj1.dcd",
        top=name + "/{:d}".format(temp) + "/top.pdb",
    )
    lz = t.unitcell_lengths[0, 2]
    edges = np.arange(0, lz, 1)
    dz = (edges[1] - edges[0]) / 2.0
    z = edges[:-1] + dz
    h = np.apply_along_axis(lambda a: np.histogram(a, bins=edges)[0], 1, t.xyz[:, :, 2])
    h1 = np.mean(h[:1000], axis=0)
    maxoverlap = np.apply_along_axis(
        lambda a: np.correlate(h1, np.histogram(a, bins=edges)[0], "full").argmax()
        - h1.size
        + dz,
        1,
        t.xyz[:, :, 2],
    )

    u = MDAnalysis.Universe(
        name + "/{:d}".format(temp) + "/top.pdb",
        name + "/{:d}".format(temp) + "/traj1.dcd",
    )
    ag = u.atoms
    with MDAnalysis.Writer(name + "/{:d}".format(temp) + "/traj.dcd", ag.n_atoms) as W:
        for ts, mo in zip(u.trajectory, maxoverlap):
            ts = transformations.unwrap(ag)(ts)
            ts = transformations.translate([0, 0, mo * 10])(ts)
            ts = transformations.wrap(ag)(ts)
            W.write(ag)

    t = md.load(
        name + "/{:d}".format(temp) + "/traj.dcd",
        top=name + "/{:d}".format(temp) + "/top.pdb",
    )
    h = np.apply_along_axis(lambda a: np.histogram(a, bins=edges)[0], 1, t.xyz[:, :, 2])
    np.save("{:s}_{:d}.npy".format(name, temp), h, allow_pickle=False)
    os.remove(name + "/{:d}".format(temp) + "/traj1.dcd")
    os.remove(name + "/{:d}".format(temp) + "/traj2.dcd")
    os.remove(name + "/{:d}".format(temp) + "/traj3.dcd")
    t.xyz /= 10
    t.unitcell_lengths /= 10
    t[0].save_pdb(name + "/{:d}".format(temp) + "/top.pdb")
    t.save_dcd(name + "/{:d}".format(temp) + "/traj.dcd")
