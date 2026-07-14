"""
Principal Component Analysis (PCA) for dimensionality reduction of trajectories.
"""

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import pca, align


def simple_pca(
    universe: mda.Universe,
    selection: str = "backbone",
    n_components: int = 2,
    align_first: bool = True,
) -> dict:
    """
    Run Principal Component Analysis on a single MDAnalysis trajectory.

    Performs dimensionality reduction to identify the major modes of
    structural motion in the trajectory.

    Parameters
    ----------
    universe : mda.Universe
        The molecular dynamics universe/trajectory to analyze.
    selection : str, optional
        MDAnalysis atom selection string (default: "backbone").
        Common examples:
            - "backbone": all backbone atoms (N, CA, C, O)
            - "name CA": only C-alpha atoms
            - "backbone and resid 1-100": subset of residues
            - "protein and not name H*": all protein atoms except hydrogens
    n_components : int, optional
        Number of principal components to compute (default: 2).
        Higher values retain more information but reduce interpretability.
    align_first : bool, optional
        Whether to align the trajectory to its first frame before PCA
        (default: True). Removes global rotation/translation artifacts.
        Highly recommended for meaningful results.

    Returns
    -------
    dict
        Dictionary containing:
        - "pca" : MDAnalysis.analysis.pca.PCA object with full analysis results
        - "transformed" : np.ndarray, shape (n_frames, n_components)
          Trajectory coordinates projected onto principal components

    Examples
    --------
    >>> result = simple_pca(universe, selection="name CA", n_components=3)
    >>> coords = result["transformed"]  # Plot or analyze PC space
    >>> pca_obj = result["pca"]         # Access variance info

    Notes
    -----
    - Alignment uses C-alpha atoms if not otherwise specified for selection.
    - PCA is performed on Cartesian coordinates of selected atoms.
    - Results should be interpreted in context of protein structure.
    """

    # Optionally align trajectory to remove rotational/translational motion
    if align_first:
        # Reload fresh reference universe from first frame
        ref = mda.Universe(universe.filename)
        align.AlignTraj(universe, ref, select=selection, in_memory=True).run()

    # Perform PCA on the selected atoms
    pc_analysis = pca.PCA(universe, select=selection, n_components=n_components)
    pc_analysis.run()

    # Project trajectory onto principal components
    ag = universe.select_atoms(selection)
    transformed = pc_analysis.transform(ag, n_components=n_components)

    return {
        "pca": pc_analysis,
        "transformed": transformed
    }