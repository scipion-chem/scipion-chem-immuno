"""N-glycosylation sequon scanning, dataset.txt construction and prediction parsing.

StackGlyEmbed's embeddings (ProteinBERT global + ESM-2 windowed + ProtT5
residue-point) are computed on the FULL parent protein sequence, not on an
isolated peptide fragment: the windowing (+-15 aa around each site) needs
real flanking sequence context, and a short candidate peptide would give a
different (wrong) embedding than the same site evaluated in its native
protein context. This is why this plugin groups sites by parent sequence
and evaluates each unique parent sequence's full-length embedding only
once, even when several input ROIs share the same parent.
"""

import re
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd

from ..constants import STACKGLYEMBED_SEQUON_PATTERN as SEQUON_PATTERN

_SEQUON_RE = re.compile(SEQUON_PATTERN)


def scanSequons(seq: str) -> List[int]:
    """Find every N-X-[S/T] sequon in ``seq``.

    Lookahead-based (not a plain ``finditer``) so overlapping sequons are
    all found (e.g. 'NPNSTPNST' reports N at position 1 AND 7).

    Returns:
        1-indexed positions of the Asn ('N') residue, relative to ``seq``.
    """
    return [m.start() + 1 for m in _SEQUON_RE.finditer(seq)]


def buildDataset(parentSites: Dict[str, List[int]], datasetPath: Path) -> List[Tuple[str, int]]:
    """Write ``dataset.txt`` (StackGlyEmbed's original format).

    Args:
        parentSites: ``{parent_sequence: [absolute_site_position, ...]}``,
            one entry per UNIQUE parent sequence (already deduplicated by
            the caller), each list already sorted/deduplicated.
        datasetPath: Output path for ``dataset.txt``.

    Returns:
        Flat, ordered list of ``(parent_sequence, site_position)`` matching
        the exact row order ``predict_local.py`` will produce in
        ``predicted_values.csv`` -- needed to map predictions back.
    """
    order: List[Tuple[str, int]] = []
    with open(datasetPath, 'w') as fh:
        for i, (seq, sites) in enumerate(parentSites.items()):
            proteinId = f'protein_{i}'
            fh.write(f'{proteinId},' + ','.join(str(s) for s in sites) + '\n')
            fh.write(f'{seq}\n')
            order.extend((seq, s) for s in sites)
    return order


class StackGlyEmbedParseError(Exception):
    """The StackGlyEmbed prediction output does not match the expected format."""


def parsePredictions(predictedPath: Path, order: List[Tuple[str, int]]) -> pd.DataFrame:
    """Parse ``predicted_values.csv`` and attach it back to ``order``.

    Returns:
        DataFrame with columns ``sequence`` (parent), ``site_position``
        (absolute, 1-indexed Asn position), ``stackglyembed_verdict``
        (``'Glycosylated'``/``'Not glycosylated'``), ``stackglyembed_score``.
    """
    try:
        raw = pd.read_csv(predictedPath)
    except Exception as exc:
        raise StackGlyEmbedParseError(f"Could not parse StackGlyEmbed output at '{predictedPath}': {exc}") from exc

    if {'prediction', 'probability'} - set(raw.columns):
        raise StackGlyEmbedParseError(
            f"StackGlyEmbed output format does not match what was expected: columns found {list(raw.columns)}."
        )
    if len(raw) != len(order):
        raise StackGlyEmbedParseError(
            f"StackGlyEmbed returned {len(raw)} row(s) but {len(order)} site(s) were submitted."
        )

    return pd.DataFrame({
        'sequence': [seq for seq, _ in order],
        'site_position': [pos for _, pos in order],
        'stackglyembed_verdict': raw['prediction'].apply(lambda p: 'Glycosylated' if int(p) == 1 else 'Not glycosylated'),
        'stackglyembed_score': raw['probability'].to_numpy(),
    })
