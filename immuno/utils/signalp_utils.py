"""Parsing of SignalP-6.0 output, including a real comment-line-count bug fix."""

from pathlib import Path

import pandas as pd

_COLUMNS = ['ID', 'Prediction', 'OTHER', 'SP', 'LIPO', 'TAT', 'TATLIPO', 'PILIN', 'CS_Position']


class SignalPParseError(Exception):
    """The SignalP-6.0 output does not match the expected format."""


def parse_output(results_path: Path, n_expected: int) -> pd.DataFrame:
    """Parse SignalP-6.0's ``prediction_results.txt``.

    ``prediction_results.txt`` carries a VARIABLE number of ``#``-prefixed
    comment lines before the data rows (a version header + a column-header
    line, 2 in practice, but a fixed ``skiprows=1`` misread the second
    comment line as a data row -- a real bug found this way). Using
    ``comment='#'`` skips ALL comment lines regardless of count, instead of
    assuming a fixed number.

    Returns:
        DataFrame with columns ``signalp_prediction`` (``'OTHER'`` if no
        signal peptide, ``'SP(Sec/SPI)'``/``'LIPO(Sec/SPII)'``/etc.
        otherwise), ``signalp_prob_other``, ``signalp_prob_sp``,
        ``signalp_cs_position`` (cleavage site position, empty if
        prediction is ``'OTHER'``).
    """
    try:
        raw = pd.read_csv(results_path, sep='\t', comment='#', header=None)
    except Exception as exc:
        raise SignalPParseError(f"Could not parse SignalP-6.0 output at '{results_path}': {exc}") from exc

    if raw.shape[1] != len(_COLUMNS):
        raise SignalPParseError(
            f"SignalP-6.0 output format does not match what was expected: found {raw.shape[1]} column(s), "
            f"expected {len(_COLUMNS)}."
        )
    raw.columns = _COLUMNS

    if len(raw) != n_expected:
        raise SignalPParseError(f"SignalP-6.0 returned {len(raw)} prediction(s), {n_expected} were expected.")

    return pd.DataFrame({
        'signalp_prediction': raw['Prediction'].tolist(),
        'signalp_prob_other': raw['OTHER'].tolist(),
        'signalp_prob_sp': raw['SP'].tolist(),
        'signalp_cs_position': raw['CS_Position'].tolist(),
    })
