"""Generic gap-tolerant sliding-window primitives, shared across the
per-residue-score-based epitope predictors in this plugin (ScanNet,
DiscoTope-3.0): each tool's own '<tool>_epitope_mapping.py' wraps these
with its own residue-grouping shape (ScanNet: per-chain; DiscoTope-3.0:
single-chain per run) rather than duplicating the sliding-window logic
itself.
"""

from typing import List, Tuple


def find_valid_windows(
    scores: List[float], threshold: float, window_size: int, max_gap_residues: int
) -> List[Tuple[int, int]]:
    """Slide a ``window_size`` window (step=1) and return the valid ranges.

    A window ``[i, i + window_size - 1]`` (0-indexed, inclusive) is valid if,
    at once: (a) at most ``max_gap_residues`` of its residues have an
    individual score below ``threshold``, and (b) the window's mean score is
    ``>= threshold``.
    """
    n = len(scores)
    valid_windows = []
    for i in range(0, n - window_size + 1):
        window = scores[i: i + window_size]
        below_count = sum(1 for score in window if score < threshold)
        if below_count <= max_gap_residues and (sum(window) / window_size) >= threshold:
            valid_windows.append((i, i + window_size - 1))
    return valid_windows


def merge_overlapping_windows(windows: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Merge overlapping/adjacent valid windows into contiguous regions.

    Assumes ``windows`` is ordered by start position (guaranteed by
    :func:`find_valid_windows`'s sequential sweep).
    """
    if not windows:
        return []

    merged = [windows[0]]
    for start, end in windows[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end + 1:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged
