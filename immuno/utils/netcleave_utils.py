"""Parsing of NetCleave output and C-terminal cleavage cross-reference.

NetCleave is a complementary signal, not a filter: a peptide can bind
MHC-I/II strongly but never actually be generated via real antigen
processing if the proteasome does not cut exactly where the accepted
binding core ends.
"""

from pathlib import Path
from typing import Optional

import pandas as pd

_RAW_REQUIRED_COLUMNS = {
    "Cleavage site", "Cleavage site after position",
    "Cleavage site after residue", "Cleavage site prediction score",
}

_OUTPUT_COLUMNS = ["sequence_window", "cleavage_position", "cleavage_residue", "cleavage_score"]


class NetCleaveParseError(Exception):
    """The NetCleave output .xlsx does not match the expected format."""


def parseOutput(xlsxPath: Path) -> pd.DataFrame:
    """Parse a single NetCleave ``--score_fasta`` .xlsx output file.

    Returns:
        DataFrame with columns ``sequence_window`` (7-residue window around
        the cut, format ``XXXX|XXX``), ``cleavage_position`` (1-indexed,
        residue IMMEDIATELY AFTER the cut), ``cleavage_residue``,
        ``cleavage_score`` (0-1, raw NetCleave score, higher = more likely cut).
    """
    try:
        raw = pd.read_excel(xlsxPath)
    except Exception as exc:
        raise NetCleaveParseError(f"Could not parse NetCleave output at '{xlsxPath}': {exc}") from exc

    if not _RAW_REQUIRED_COLUMNS.issubset(raw.columns):
        raise NetCleaveParseError(
            f"NetCleave output .xlsx format does not match what was expected: missing columns "
            f"{_RAW_REQUIRED_COLUMNS - set(raw.columns)}. Columns found: {list(raw.columns)}."
        )

    return pd.DataFrame({
        'sequence_window': raw['Cleavage site'],
        'cleavage_position': raw['Cleavage site after position'],
        'cleavage_residue': raw['Cleavage site after residue'],
        'cleavage_score': raw['Cleavage site prediction score'],
    })[_OUTPUT_COLUMNS]


def findCTermMatch(cleavageDf: pd.DataFrame, coreSeq: str, windowSeq: str) -> Optional[float]:
    """Check whether ``cleavageDf`` has a cut exactly after ``coreSeq`` inside ``windowSeq``.

    Args:
        cleavageDf: Output of :func:`parse_output`, evaluated on ``windowSeq``.
        coreSeq: The accepted MHC-binding core (may equal ``windowSeq`` itself
            if the caller has no narrower core, e.g. no upstream
            ``_core9aa`` attribute available).
        windowSeq: The full peptide NetCleave was run on (needs real
            flanking context around the cut, not just the isolated core).

    Returns:
        The best (highest) ``cleavage_score`` among exact C-terminal
        matches, or ``None`` if no cut falls exactly one residue after
        ``coreSeq``'s last residue within ``windowSeq``.
    """
    if cleavageDf.empty:
        return None

    offset = windowSeq.find(coreSeq)
    if offset == -1:
        return None
    cTermPosition = offset + len(coreSeq) + 1

    matches = cleavageDf[cleavageDf['cleavage_position'] == cTermPosition]
    if matches.empty:
        return None
    return float(matches['cleavage_score'].max())
