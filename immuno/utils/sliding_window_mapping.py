"""Generic gap-tolerant sliding-window epitope region mapping, shared across
the per-residue-score-based epitope predictors in this plugin (EpiDope,
DiscoTope-3.0, ScanNet).

Two shapes of consumer:

* a single sequence/chain per run (EpiDope on a protein sequence,
  DiscoTope-3.0 on a single-chain PDB) uses
  :func:`extractLinearEpitopeRegions` directly -- the mapping is identical
  for both, so neither gets a per-tool wrapper that would only duplicate it;
* several chains in one run (ScanNet, where a PDB can hold more than one
  chain) keeps its own 'scannet_epitope_mapping.py' wrapper, because the
  per-chain grouping is real tool-specific shape and not a copy of the
  sliding-window logic.
"""

from typing import List, Tuple

import pandas as pd

REGION_COLUMNS = ['start', 'end', 'length', 'mean_score', 'max_score', 'sequence']


def findValidWindows(
    scores: List[float], threshold: float, windowSize: int, maxGapResidues: int
) -> List[Tuple[int, int]]:
    """Slide a ``windowSize`` window (step=1) and return the valid ranges.

    A window ``[i, i + windowSize - 1]`` (0-indexed, inclusive) is valid if,
    at once: (a) at most ``maxGapResidues`` of its residues have an
    individual score below ``threshold``, and (b) the window's mean score is
    ``>= threshold``.
    """
    n = len(scores)
    validWindows = []
    for i in range(0, n - windowSize + 1):
        window = scores[i: i + windowSize]
        belowCount = sum(1 for score in window if score < threshold)
        if belowCount <= maxGapResidues and (sum(window) / windowSize) >= threshold:
            validWindows.append((i, i + windowSize - 1))
    return validWindows


def mergeOverlappingWindows(windows: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Merge overlapping/adjacent valid windows into contiguous regions.

    Assumes ``windows`` is ordered by start position (guaranteed by
    :func:`findValidWindows`'s sequential sweep).
    """
    if not windows:
        return []

    merged = [windows[0]]
    for start, end in windows[1:]:
        lastStart, lastEnd = merged[-1]
        if start <= lastEnd + 1:
            merged[-1] = (lastStart, max(lastEnd, end))
        else:
            merged.append((start, end))
    return merged


def extractLinearEpitopeRegions(
    scores: List[float], residues: List[str], threshold: float, minLength: int,
    windowSize: int, maxGapResidues: int,
) -> pd.DataFrame:
    """Maps epitope regions with a gap-tolerant sliding window over ONE sequence.

    Used as-is by EpiDope (per-residue score of a protein sequence) and by
    DiscoTope-3.0 (per-residue 'calibrated_score' of a single-chain PDB).
    Residue position is derived from list order (1-indexed), never from a
    tool's own position column.

    Returns:
        DataFrame with columns ``start``, ``end``, ``length``,
        ``mean_score``, ``max_score``, ``sequence``.
    """
    validWindows = findValidWindows(scores, threshold, windowSize, maxGapResidues)
    mergedRegions = mergeOverlappingWindows(validWindows)

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

    return pd.DataFrame.from_records(records, columns=REGION_COLUMNS)
