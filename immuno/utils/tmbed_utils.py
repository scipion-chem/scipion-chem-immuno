"""Parsing of TMbed's 3-line prediction output (see ``ProtTMbedPredict`` for
how this is used).

TMbed's ``--out-format 1`` emits, per protein, three lines: FASTA header,
amino acid sequence, and a same-length string of per-residue class letters:
'B'/'H' for transmembrane beta-strand/alpha-helix residues, 'S' for signal
peptide, and 'i'/'o' for non-membrane residues on the inside/outside of the
membrane. Only B/H/S are collapsed into masking regions here: 'i'/'o'
residues are already presumed solvent-accessible and are not meant to be
excluded from downstream epitope candidates.
"""

from typing import Dict, List, Tuple

from .tmbed_exceptions import TMbedParseError

# Only these classes are reported as masking regions; 'i'/'o' (non-membrane
# loop, inside/outside the membrane) are left untouched.
MASKED_CLASS_TYPES = {
    'B': 'TM_beta_strand',
    'H': 'TM_alpha_helix',
    'S': 'signal_peptide',
}


def parsePredictions(predPath: str) -> Dict[str, Tuple[str, str]]:
    """Parse a TMbed 3-line prediction file into {header: (sequence, classes)}.

    Raises TMbedParseError if the file is empty, truncated, or a
    sequence/prediction length mismatch is found (should not happen unless
    TMbed's own output contract changes).
    """
    with open(predPath) as fh:
        lines = [line.rstrip('\n') for line in fh if line.strip()]

    if not lines or len(lines) % 3 != 0:
        raise TMbedParseError(
            f"TMbed output at '{predPath}' does not have the expected 3-line-per-protein format "
            f"(header/sequence/prediction): found {len(lines)} non-empty line(s)."
        )

    records = {}
    for i in range(0, len(lines), 3):
        header, sequence, classes = lines[i], lines[i + 1], lines[i + 2]
        if not header.startswith('>'):
            raise TMbedParseError(f"Expected a FASTA header at line {i + 1} of '{predPath}', got: '{header}'.")
        if len(sequence) != len(classes):
            raise TMbedParseError(
                f"Sequence/prediction length mismatch for '{header}' in '{predPath}': "
                f"{len(sequence)} residue(s) vs {len(classes)} class letter(s)."
            )
        records[header[1:].strip()] = (sequence, classes)

    return records


def extractMaskingRegions(classes: str, minLength: int = 1) -> List[dict]:
    """Collapse a per-residue class string into contiguous masking regions.

    Only residues classified as 'B' (TM beta-strand), 'H' (TM alpha-helix)
    or 'S' (signal peptide) are turned into regions; 'i'/'o' (non-membrane
    loop) residues are skipped. Returns a list of dicts with 1-indexed
    'start', 'end' and 'type' (see MASKED_CLASS_TYPES), in sequence order,
    dropping any region shorter than minLength.
    """
    regions = []
    curType, start = None, None
    n = len(classes)
    for i in range(n):
        roiType = MASKED_CLASS_TYPES.get(classes[i])
        pos = i + 1
        if roiType != curType:
            if curType is not None and (pos - start) >= minLength:
                regions.append({'start': start, 'end': pos - 1, 'type': curType})
            curType = roiType
            start = pos if roiType is not None else None
    if curType is not None and (n + 1 - start) >= minLength:
        regions.append({'start': start, 'end': n, 'type': curType})
    return regions
