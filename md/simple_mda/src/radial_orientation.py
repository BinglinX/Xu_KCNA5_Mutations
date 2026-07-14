"""
Analysis of sidechain orientation relative to pore axis (ion channels).
"""

import numpy as np


def radial_orientation(u, resid,
                       axis_selection,
                       pore_selection,
                       return_angle=True):
    """
    Compute radial inward/outward orientation of a residue sidechain.

    Analyzes how a residue's sidechain is oriented relative to the pore
    axis of an ion channel. The pore axis is defined by the principal
    component of selected atoms (typically S6 helices), and sidechain
    orientation is measured relative to the radial direction at each frame.

    Parameters
    ----------
    u : MDAnalysis.Universe
        The molecular dynamics universe/trajectory.
    resid : int
        Residue ID of the residue to analyze.
    axis_selection : str
        MDAnalysis selection string for atoms defining the pore axis
        (e.g., "segid PROA and resnum 318 to 325 and name CA").
    pore_selection : str
        MDAnalysis selection string for atoms marking the pore center
        at each axial level (e.g., "segid PROA and name CA").
    return_angle : bool, optional
        If True (default), return angle in degrees (0-180).
        If False, return inwardness metric (-cosθ, range -1 to 1).

    Returns
    -------
    np.ndarray
        Array of orientation values per frame. Shape (n_frames,).
        If return_angle=True: angles in degrees (0-180).
        If return_angle=False: inwardness values (-1 to 1, where
        -1 = pointing outward, 0 = perpendicular, 1 = pointing inward).

    Notes
    -----
    - The pore axis is computed from the first frame only (static).
    - Frames where sidechain centroid coincides with CA are skipped.
    - Sidechain = all atoms except backbone (N, C, O, CA, H).
    """

    # Compute pore axis from first frame (treated as static)
    u.trajectory[0]

    # Select atoms defining the pore axis and compute principal component
    axis_atoms = u.select_atoms(axis_selection)
    coords = axis_atoms.positions

    # Center coordinates at origin
    center = coords.mean(axis=0)
    coords_centered = coords - center

    # Compute covariance matrix and eigendecomposition
    cov = np.cov(coords_centered.T)
    eigvals, eigvecs = np.linalg.eigh(cov)

    # Principal component (highest variance) defines the pore axis
    axis = eigvecs[:, np.argmax(eigvals)]
    axis /= np.linalg.norm(axis)

    # Select target residue, its CA, and sidechain atoms
    res = u.select_atoms(f"resid {resid}").residues[0]
    ca = res.atoms.select_atoms("name CA")
    sidechain = res.atoms.select_atoms("not name N C O CA H*")

    # Select atoms defining pore center at each frame
    pore_atoms = u.select_atoms(pore_selection)

    results = []

    # Analyze each frame in the trajectory
    for ts in u.trajectory:
        # Pore center at this axial level (mean of pore atom positions)
        pore_center = pore_atoms.positions.mean(axis=0)

        # Radial vector from pore center to CA (remove axial component)
        r = ca.positions[0] - pore_center
        r_parallel = np.dot(r, axis) * axis
        r_radial = r - r_parallel
        r_norm = np.linalg.norm(r_radial)
        
        # Skip frame if CA is on axis (no well-defined radial direction)
        if r_norm == 0:
            continue
        r_radial /= r_norm

        # Sidechain vector (from CA to sidechain centroid)
        sc_centroid = sidechain.center_of_mass()
        sc_vec = sc_centroid - ca.positions[0]

        # Remove axial component (keep only radial)
        sc_parallel = np.dot(sc_vec, axis) * axis
        sc_radial = sc_vec - sc_parallel
        sc_norm = np.linalg.norm(sc_radial)
        
        # Skip frame if sidechain centroid is on axis
        if sc_norm == 0:
            continue
        sc_radial /= sc_norm

        # Compute angle between radial direction and sidechain vector
        cos_theta = np.dot(sc_radial, r_radial)
        cos_theta = np.clip(cos_theta, -1.0, 1.0)

        # Convert to requested output format
        if return_angle:
            # Angle in degrees (0 = aligned, 90 = perpendicular, 180 = opposite)
            value = np.degrees(np.arccos(cos_theta))
        else:
            # Inwardness: -cosθ (1 = inward, 0 = perpendicular, -1 = outward)
            value = -cos_theta

        results.append(value)

    return np.array(results)