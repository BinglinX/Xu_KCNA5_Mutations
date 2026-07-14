"""
Analysis of helix geometry and properties in trajectories.
"""

from MDAnalysis.analysis import helix_analysis as hel


def simple_helix(u, filter, axis=[0, 0, 1]):
    """
    Compute helix geometry parameters (rise, twist, radius) for selected atoms.

    Analyzes the geometry of alpha-helical regions using MDAnalysis' HELANAL
    method. Computes helix parameters like rise per residue, twist per residue,
    and radius across the trajectory.

    Parameters
    ----------
    u : MDAnalysis.Universe
        The molecular dynamics universe/trajectory.
    filter : str
        MDAnalysis selection string for the helix region.
        Example: "segid PROA and resnum 318 to 325"
    axis : array-like, optional
        Reference axis for helix orientation (default: [0, 0, 1]).
        Typically the Z-axis for membrane proteins.

    Returns
    -------
    dict
        Results dictionary from HELANAL containing:
        - 'radius': helix radius (Angstrom)
        - 'rise': rise per residue (Angstrom)
        - 'twist': twist per residue (degrees)
        - 'n_residues': number of residues
        Plus additional geometric parameters per frame.

    Notes
    -----
    - Automatically selects C-alpha atoms from the provided filter.
    - Helix geometry is analyzed independently for each frame.
    """

    # Run HELANAL on selected atoms (automatically uses CA)
    h = hel.HELANAL(u, select=f"{filter} and name CA", ref_axis=axis).run()

    return h.results