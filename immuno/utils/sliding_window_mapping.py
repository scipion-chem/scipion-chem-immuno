"""Generic gap-tolerant sliding-window primitives, shared across the
per-residue-score-based epitope predictors in this plugin (ScanNet,
DiscoTope-3.0): each tool's own '<tool>_epitope_mapping.py' wraps these
with its own residue-grouping shape (ScanNet: per-chain; DiscoTope-3.0:
single-chain per run) rather than duplicating the sliding-window logic
itself.
"""

from typing import List, Tuple


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
