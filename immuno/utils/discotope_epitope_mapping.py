"""Generic gap-tolerant sliding-window epitope region mapping.

Generic gap-tolerant sliding-window mapping logic, shared in principle
across BepiPred-3.0/EpiDope/DiscoTope-3.0/ScanNet-style per-residue score
outputs, applied here on DiscoTope-3.0's per-residue 'calibrated_score'
for a single-chain PDB per run.
"""

from typing import List, Tuple

import pandas as pd


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


def extract_epitope_regions(
    scores: List[float], residues: List[str], threshold: float, minLength: int,
    windowSize: int, maxGapResidues: int,
) -> pd.DataFrame:
    """Map epitope regions with a gap-tolerant sliding window over a single chain.

    Residue position is derived from list order (1-indexed), not from a
    position column.

    Returns:
        DataFrame with columns ``start``, ``end``, ``length``,
        ``mean_score``, ``max_score``, ``sequence``.
    """
    validWindows = find_valid_windows(scores, threshold, windowSize, maxGapResidues)
    mergedRegions = merge_overlapping_windows(validWindows)

    records = []
    for start, end in mergedRegions:
        length = end - start + 1
        if length < minLength:
            continue

        blockScores = scores[start: end + 1]
        blockResidues = residues[start: end + 1]
        records.append({
            'start': start + 1,
            'end': end + 1,
            'length': length,
            'mean_score': sum(blockScores) / len(blockScores),
            'max_score': max(blockScores),
            'sequence': ''.join(blockResidues),
        })

    return pd.DataFrame.from_records(
        records, columns=['start', 'end', 'length', 'mean_score', 'max_score', 'sequence'],
    )
