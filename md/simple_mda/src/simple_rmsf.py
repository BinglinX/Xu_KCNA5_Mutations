"""
Root Mean Square Fluctuation (RMSF) analysis of molecular dynamics trajectories.
"""

import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np
from MDAnalysis.analysis import rms, align


def simple_rmsf(u, chain=None):
    """
    Calculate C-alpha Root Mean Square Fluctuation (RMSF) per residue.

    Computes RMSF by aligning the trajectory to its average structure and
    measuring the per-residue C-alpha displacement. RMSF quantifies residue
    flexibility/dynamism across the trajectory.

    Parameters
    ----------
    u : MDAnalysis.Universe
        The molecular dynamics universe/trajectory to analyze.
    chain : str, optional
        Chain identifier ('A', 'B', 'C', or 'D') to analyze single chain.
        If None, analyzes all chains with protein atoms.
        Maps to internal segids: A->seg_0_PROA, B->seg_1_PROB, etc.

    Returns
    -------
    resids : np.ndarray
        Residue IDs corresponding to each RMSF value. Shape: (n_residues,).
    R : MDAnalysis.analysis.rms.RMSF
        RMSF analysis object containing per-residue fluctuations.
        Access values via R.results.rmsf or R.results.rmsf_std.

    Examples
    --------
    >>> resids, rmsf_obj = simple_rmsf(u, chain='A')
    >>> rmsf_values = rmsf_obj.results.rmsf
    >>> flexible_residues = resids[rmsf_values > 2.0]  # > 2 Angstrom

    Notes
    -----
    - RMSF = sqrt( <(r_i(t) - <r_i>)^2> )
    - Higher RMSF indicates more flexible/dynamic residue
    - Typical flexible regions: loops, termini; rigid: secondary structure
    - Segids are hardcoded for 4-chain system; may need customization.
    """

    # Map chain letters to internal MDAnalysis segids
    segid_list = {
        "A": "seg_0_PROA",
        "B": "seg_1_PROB",
        "C": "seg_2_PROC",
        "D": "seg_3_PROD"
    }

    # Build selection string based on chain parameter
    if chain:
        filter = f"protein and name CA and segid {segid_list[chain]}"
    else:
        # If no chain specified, select all protein CA atoms
        filter = "protein and name CA"

    # Step 1: Compute average structure (reference for RMSF)
    average = align.AverageStructure(u, u, select=filter,
                                    ref_frame=0).run()
    ref = average.results.universe

    # Step 2: Align trajectory to average structure
    # (removes rotation/translation to focus on internal fluctuations)
    aligner = align.AlignTraj(u, ref,
                          select=filter,
                          in_memory=True).run()

    # Step 3: Calculate RMSF for each residue
    c_alphas = u.select_atoms(filter)
    R = rms.RMSF(c_alphas).run()

    return c_alphas.resids, R