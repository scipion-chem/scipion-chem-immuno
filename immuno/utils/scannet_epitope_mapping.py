"""Generic gap-tolerant sliding-window epitope region mapping.

Generic gap-tolerant sliding-window mapping logic, shared in principle
across BepiPred-3.0/EpiDope/DiscoTope-3.0/ScanNet-style per-residue score
outputs, applied here on ScanNet's per-residue scores grouped by CHAIN (a
single PDB can contain several chains, so this cannot assume a single
sequence per run the way a purely linear predictor would).
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
    scoresDf: pd.DataFrame,
    groupCol: str,
    scoreCol: str,
    residueCol: str,
    threshold: float,
    minLength: int,
    windowSize: int,
    maxGapResidues: int,
) -> pd.DataFrame:
    """Map epitope regions with a gap-tolerant sliding window, per ``groupCol`` value.

    Residue position is derived from row order within each group
    (1-indexed), not from a position column.

    Returns:
        DataFrame with columns ``group``, ``start``, ``end``, ``length``,
        ``mean_score``, ``max_score``, ``sequence``.
    """
    records = []
    for groupValue, group in scoresDf.groupby(groupCol, sort=False):
        group = group.reset_index(drop=True)
        scores = group[scoreCol].tolist()

        validWindows = find_valid_windows(scores, threshold, windowSize, maxGapResidues)
        mergedRegions = merge_overlapping_windows(validWindows)

        for start, end in mergedRegions:
            length = end - start + 1
            if length < minLength:
                continue

            block = group.iloc[start: end + 1]
            records.append({
                'group': groupValue,
                'start': start + 1,
                'end': end + 1,
                'length': length,
                'mean_score': float(block[scoreCol].mean()),
                'max_score': float(block[scoreCol].max()),
                'sequence': ''.join(block[residueCol].astype(str)),
            })

    return pd.DataFrame.from_records(
        records, columns=['group', 'start', 'end', 'length', 'mean_score', 'max_score', 'sequence'],
    )
