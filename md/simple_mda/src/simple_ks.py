"""
Kolmogorov-Smirnov (KS) test analysis of residue fluctuations.

Compares residue fluctuation distributions between trajectories to identify
residues with significantly different dynamics.
"""

import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ks_2samp, kstwo
from MDAnalysis.analysis import rms, align


def residue_fluctuations(u):
    """
    Compute per-residue displacement from average structure across frames.

    Aligns trajectory to average structure and measures C-alpha displacement
    for each residue at each frame.

    Parameters
    ----------
    u : MDAnalysis.Universe
        The molecular dynamics universe/trajectory.

    Returns
    -------
    np.ndarray
        Displacements per frame and residue. Shape: (n_frames, n_residues).
        Each value is the Euclidean distance from average position (Angstrom).

    Notes
    -----
    - Hardcoded for 'seg_0_PROA' segid; may need customization for other systems.
    - Alignment removes rotation and translation before computing displacement.
    """

    # Define selection and compute average structure
    filter = "protein and name CA and segid seg_0_PROA"
    average = align.AverageStructure(
        u, u, select=filter, ref_frame=0
    ).run()
    ref = average.results.universe
    
    # Align trajectory to average structure
    align.AlignTraj(u, ref, select=filter, in_memory=True).run()

    # Get C-alpha atoms and store mean positions
    ca = u.select_atoms(filter)
    mean_pos = ca.positions.copy()

    fluctuations = []

    # Compute displacement for each frame
    for ts in u.trajectory:
        disp = np.linalg.norm(ca.positions - mean_pos, axis=1)
        fluctuations.append(disp)

    return np.array(fluctuations)  # shape: (n_frames, n_residues)


def split_halves(data):
    """
    Split trajectory data into first and second halves.

    Used for jackknife-style variance estimation.

    Parameters
    ----------
    data : np.ndarray
        Data array with frames as first dimension.

    Returns
    -------
    tuple of np.ndarray
        (first_half, second_half)
    """
    n = data.shape[0] // 2
    return data[:n], data[n:]


def ks_per_residue(sample_a, sample_b):
    """
    Compute Kolmogorov-Smirnov test statistic per residue.

    Tests whether two samples of residue fluctuations come from the same
    distribution.

    Parameters
    ----------
    sample_a : np.ndarray
        Fluctuation data, shape (n_frames, n_residues).
    sample_b : np.ndarray
        Fluctuation data, shape (n_frames, n_residues).

    Returns
    -------
    np.ndarray
        KS test statistic per residue. Shape: (n_residues,).
        Higher values indicate greater difference in distributions.
    """
    n_res = sample_a.shape[1]
    ks = np.zeros(n_res)
    pvals = np.zeros(n_res)

    # Apply KS test to each residue independently
    for i in range(n_res):
        ks[i], pvals[i] = ks_2samp(
            sample_a[:, i],
            sample_b[:, i]
        )

    return ks


def simple_ks(u1, u2):
    """
    Compare residue fluctuation distributions between two trajectories.

    Performs KS test between all combinations of trajectory halves to produce
    a robust comparison statistic and p-values.

    Parameters
    ----------
    u1 : MDAnalysis.Universe
        First trajectory (e.g., wild-type).
    u2 : MDAnalysis.Universe
        Second trajectory (e.g., mutant).

    Returns
    -------
    ks : np.ndarray
        Mean KS test statistic per residue. Shape: (n_residues,).
    p : np.ndarray
        P-values per residue (significance of difference).

    Notes
    -----
    - Compares all 4 combinations of halves: 1a vs 2a, 1b vs 2a, 1a vs 2b, 1b vs 2b.
    - Results averaged across combinations for robustness.
    - P-values computed using two-sample KS distribution.
    """

    # Compute fluctuations for each trajectory
    fluct_1, fluct_2 = residue_fluctuations(u1), residue_fluctuations(u2)

    # Split each into halves for jackknife estimation
    fluct_1a, fluct_1b = split_halves(fluct_1)
    fluct_2a, fluct_2b = split_halves(fluct_2)

    # Average KS statistic across all half-pair combinations
    ks = (
        ks_per_residue(fluct_1a, fluct_2a) +
        ks_per_residue(fluct_1b, fluct_2a) +
        ks_per_residue(fluct_1a, fluct_2b) +
        ks_per_residue(fluct_1b, fluct_2b)
    ) / 4.0

    # Effective sample size (frames per half)
    n_eff = fluct_1.shape[0]

    # Compute p-values from two-sample KS distribution
    p = kstwo.sf(ks, n_eff)

    return ks, p


def simple_ks_self(u1):
    """
    Compute KS test statistic by comparing trajectory halves (self-comparison).

    Estimates the noise/variability in fluctuation measurements within a
    single trajectory by comparing first and second halves.

    Parameters
    ----------
    u1 : MDAnalysis.Universe
        The trajectory to analyze.

    Returns
    -------
    ks : np.ndarray
        KS test statistic per residue. Shape: (n_residues,).
    p : np.ndarray
        P-values per residue.

    Notes
    -----
    - Useful as a control to assess measurement noise.
    - If ks values here are non-negligible, then differences between
      trajectories near this magnitude may not be statistically significant.
    """

    # Compute fluctuations and split into halves
    fluct_1 = residue_fluctuations(u1)
    fluct_1a, fluct_1b = split_halves(fluct_1)

    # Compare halves within same trajectory
    ks = ks_per_residue(fluct_1a, fluct_1b)

    # Effective sample size
    n_eff = fluct_1.shape[0]

    # Compute p-values
    p = kstwo.sf(ks, n_eff)

    return ks, p