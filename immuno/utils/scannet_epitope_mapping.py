"""Generic gap-tolerant sliding-window epitope region mapping.

Generic gap-tolerant sliding-window mapping logic, shared in principle
across BepiPred-3.0/EpiDope/DiscoTope-3.0/ScanNet-style per-residue score
outputs, applied here on ScanNet's per-residue scores grouped by CHAIN (a
single PDB can contain several chains, so this cannot assume a single
sequence per run the way a purely linear predictor would).
"""

import pandas as pd

from .sliding_window_mapping import find_valid_windows, merge_overlapping_windows


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
