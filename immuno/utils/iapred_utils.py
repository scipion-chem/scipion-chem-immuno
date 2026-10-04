"""Output parsing for IApred, including the real <20 aa 'too short' case."""

from pathlib import Path

import pandas as pd

from ..constants import IAPRED_TOO_SHORT_CATEGORY

_RAW_REQUIRED_COLUMNS = {'Header', 'Sequence_Length', 'Intrinsic_Antigenicity_Score',
                         'Antigenicity_Category'}


class IApredParseError(Exception):
    """The IApred output CSV does not match the expected format."""


def parseOutput(csvPath: Path, nExpected: int) -> pd.DataFrame:
    """Parse IApred's raw output CSV.

    For sequences shorter than IApred's own 20 aa floor it writes
    the literal text 'Sequence too short' in the score column instead of a
    number, so the score is coerced to NaN and the category is replaced by
    an explicit string rather than left as a silent unexplained NaN.

    Returns:
        DataFrame with columns ``iapred_score`` (NaN if too short) and
        ``iapred_category``.
    """
    try:
        raw = pd.read_csv(csvPath)
    except Exception as exc:
        raise IApredParseError(f"Could not parse IApred output at '{csvPath}': {exc}") from exc

    missing = _RAW_REQUIRED_COLUMNS - set(raw.columns)
    if missing:
        raise IApredParseError(
            f"IApred output CSV does not contain the expected columns {sorted(missing)}. "
            f"Columns found: {list(raw.columns)}."
        )
    if len(raw) != nExpected:
        raise IApredParseError(f"IApred returned {len(raw)} row(s), {nExpected} were expected.")

    scores = pd.to_numeric(raw['Intrinsic_Antigenicity_Score'], errors='coerce')
    categories = [
        IAPRED_TOO_SHORT_CATEGORY if isShort and pd.isna(category) else category
        for isShort, category in zip(scores.isna(), raw['Antigenicity_Category'].tolist())
    ]

    return pd.DataFrame({'iapred_score': scores.tolist(), 'iapred_category': categories})
