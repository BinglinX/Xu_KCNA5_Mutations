"""
Utilities for loading MD data from NPZ files.
"""

import numpy as np


def load_md_npz(filename):
    """
    Load RMSF data saved in flat-key npz format.

    Parses NPZ files containing per-chain RMSF values organized by system
    and chain identifiers. File keys should follow the naming convention
    'SYSTEM_CHAIN' (e.g., 'WT_A', 'G128R_B') with optional 'resids' key
    for residue indices.

    Parameters
    ----------
    filename : str or Path
        Path to the NPZ file to load.

    Returns
    -------
    md_dict : dict
        Nested dictionary with structure {system: {chain: rmsf_array}}.
        Example: {'WT': {'A': array([...]), 'B': array([...])}, ...}
    resids : np.ndarray or None
        Residue indices if 'resids' key is present in the file, else None.

    Examples
    --------
    >>> md_dict, resids = load_md_npz('data.npz')
    >>> wt_chain_a_rmsf = md_dict['WT']['A']
    """

    data = np.load(filename)

    md_dict = {}
    resids = None

    # Parse each key in the NPZ file
    for key in data.files:
        if key == "resids":
            # Extract residue indices if present
            resids = data[key]
            continue

        # Split key into system and chain components
        system, chain = key.split("_")

        # Initialize system entry if not present
        if system not in md_dict:
            md_dict[system] = {}

        # Store the RMSF array for this system/chain
        md_dict[system][chain] = data[key]

    return md_dict, resids