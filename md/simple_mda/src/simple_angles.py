"""
Analysis of dihedral angles and Ramachandran plot data.
"""

import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Ramachandran, Dihedral


def simple_dihedrals(u, resid):
    """
    Compute Ramachandran dihedral angles (phi, psi) for a residue.

    Analyzes phi (N-CA-C-N) and psi (CA-C-N-CA) backbone dihedral angles
    across all frames in the trajectory for a single residue.

    Parameters
    ----------
    u : MDAnalysis.Universe
        The molecular dynamics universe/trajectory.
    resid : int
        Residue ID to analyze.

    Returns
    -------
    MDAnalysis.analysis.dihedrals.Ramachandran
        Ramachandran analysis object containing:
        - results.angles: phi and psi angles per frame
        - results.times: time values for each frame

    Notes
    -----
    Uses MDAnalysis built-in Ramachandran analysis, which automatically
    computes standard backbone dihedral angles.
    """

    # Select the residue of interest
    res = u.select_atoms(f"resid {resid}")

    # Run Ramachandran analysis on trajectory
    R = Ramachandran(res).run()

    return R

