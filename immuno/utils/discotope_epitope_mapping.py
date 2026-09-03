"""Generic gap-tolerant sliding-window epitope region mapping.

Generic gap-tolerant sliding-window mapping logic, shared in principle
across BepiPred-3.0/EpiDope/DiscoTope-3.0/ScanNet-style per-residue score
outputs, applied here on DiscoTope-3.0's per-residue 'calibrated_score'
for a single-chain PDB per run.
"""

from typing import List

import pandas as pd

from .sliding_window_mapping import find_valid_windows, merge_overlapping_windows


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
