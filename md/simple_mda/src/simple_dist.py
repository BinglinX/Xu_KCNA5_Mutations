"""
Utilities for computing inter-atomic distances in MD trajectories.
"""

import MDAnalysis as mda
from MDAnalysis.analysis import align
import numpy as np
import pickle


def simple_dist_any(u, res1, res2, atom1, atom2, end_frame=None):
    """
    Compute distance between two specific atoms across trajectory frames.

    Aligns the trajectory to C-alpha atoms before calculating distances to
    remove rotational/translational motion artifacts.

    Parameters
    ----------
    u : MDAnalysis.Universe
        The molecular dynamics universe/trajectory.
    res1 : int
        Residue ID for the first atom.
    res2 : int
        Residue ID for the second atom.
    atom1 : str
        Atom name in first residue (e.g., 'CA', 'CB', 'OG').
    atom2 : str
        Atom name in second residue.
    end_frame : int, optional
        If specified, only analyze up to this frame index.
        If None, analyze entire trajectory.

    Returns
    -------
    np.ndarray
        Array of distances (in Angstrom) per frame. Shape (n_frames,).

    Notes
    -----
    - Trajectory is aligned to protein C-alpha atoms before distance calculation.
    - Returns Euclidean distance without periodic boundary corrections.
    """

    # Align trajectory to C-alpha atoms (remove rotation/translation)
    align.AlignTraj(
        u, u,
        select="protein and name CA",
        in_memory=True
    ).run()

    # Select the two atoms of interest
    atom1_atoms = u.select_atoms(f"resid {res1} and name {atom1}")
    atom2_atoms = u.select_atoms(f"resid {res2} and name {atom2}")

    distances = []
    
    # Iterate over frames (optionally truncated at end_frame)
    trajectory = u.trajectory[:end_frame]
    for ts in trajectory:
        # Compute Euclidean distance between atom positions
        d = np.linalg.norm(atom1_atoms.positions - atom2_atoms.positions)
        distances.append(d)

    return np.array(distances)


def system_dist(systems, res1, res2, atom1, atom2, save=None, end_ns=None):
    """
    Compute inter-atomic distance for multiple trajectory systems/repeats.

    Useful for comparing the same distance measurement across different
    simulations (e.g., WT vs mutant, different simulation runs).

    Parameters
    ----------
    systems : list of MDAnalysis.Universe
        List of universe objects to analyze (e.g., different trajectories).
    res1 : int
        Residue ID for the first atom.
    res2 : int
        Residue ID for the second atom.
    atom1 : str
        Atom name in first residue (e.g., 'CA', 'CB').
    atom2 : str
        Atom name in second residue.
    save : str, optional
        If specified, save results to pickle file at this path.
    end_ns : float, optional
        If specified, truncate each trajectory to this nanosecond value.
        Conversion: frame_index = end_ns * 1000 / dt (assumes dt in ps).

    Returns
    -------
    list of np.ndarray
        List of distance arrays, one per system. Each array has shape (n_frames,).

    Examples
    --------
    >>> universes = [u_wt, u_mutant]
    >>> distances = system_dist(
    ...     universes, res1=100, res2=150, atom1='CA', atom2='CA',
    ...     end_ns=100.0, save='distances.pkl'
    ... )
    """

    repeat_dists = []

    # Analyze each system independently
    for u in systems:
        # Convert time in nanoseconds to frame index
        if end_ns is not None:
            end_frame = int(end_ns * 1000 / u.trajectory.dt)
        else:
            end_frame = None

        # Compute distances for this system
        d = simple_dist_any(u, res1, res2, atom1=atom1, atom2=atom2, end_frame=end_frame)
        repeat_dists.append(d)

    # Optionally save results to file
    if save:
        with open(save, "wb") as f:
            pickle.dump(repeat_dists, f)

    return repeat_dists