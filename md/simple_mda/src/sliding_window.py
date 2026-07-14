"""
Smoothing utilities for analysis data using sliding window averaging.
"""

import numpy as np


def sliding_window_average(rmsf, window=5):
    """
    Compute centered sliding-window average of per-residue values.

    Smooths noisy per-residue data (e.g., RMSF, B-factors) to highlight
    regional trends while reducing high-frequency noise.

    Parameters
    ----------
    rmsf : array-like
        Values per residue (length N). Typically RMSF, B-factor, or similar
        per-residue metric. Will be converted to numpy array.
    window : int, optional
        Window size in residues (default: 5).
        Must be odd to maintain symmetry (recommended: 3, 5, 7, ...).
        Larger values produce more smoothing but lose local resolution.

    Returns
    -------
    residues : np.ndarray
        Residue indices corresponding to averaged values (1-based indexing).
        Shape: (N - window + 1,) where N = len(rmsf).
    rmsf_avg : np.ndarray
        Sliding-window averaged values. Shape: same as residues.

    Examples
    --------
    >>> rmsf_values = np.array([1.2, 1.5, 1.4, 1.8, 1.6, 1.5, 1.4])
    >>> res, avg = sliding_window_average(rmsf_values, window=3)
    >>> res  # [2, 3, 4, 5, 6]
    >>> avg  # [1.37, 1.57, 1.60, 1.63, 1.50]

    Notes
    -----
    - Uses uniform kernel (simple moving average).
    - Output is shorter by (window - 1) points due to 'valid' convolution.
    - Residue numbering is 1-based (starts at 1, not 0).
    - For window=5: first output residue is index 1+2=3, last is N-2.
    """

    # Ensure input is numpy array
    rmsf = np.asarray(rmsf)

    # Create uniform kernel for averaging
    kernel = np.ones(window) / window
    
    # Apply convolution with 'valid' mode (no zero-padding)
    rmsf_avg = np.convolve(rmsf, kernel, mode="valid")

    # Compute residue indices (1-based, centered on window)
    half = window // 2
    residues = np.arange(1 + half, len(rmsf) - half + 1)

    return residues, rmsf_avg