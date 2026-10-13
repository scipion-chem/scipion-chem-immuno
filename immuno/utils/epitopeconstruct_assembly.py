"""Multi-epitope construct selection and assembly logic.

Pure selection/assembly logic, operating on Scipion ``SequenceROI``
objects (already annotated by upstream protocols). See this plugin's
protocol docstring for the full design rationale (linker choices, block
order, why no cross-class fusion, etc.).

PROJECT CONVENTION: any
Scipion-Chem protocol that predicts B-cell epitope regions and outputs a
``SetOfSequenceROIs`` MUST expose ``_meanScore`` (``Float``, required) and
SHOULD expose ``_maxScore`` (``Float``, optional) on every output ROI --
this is how this protocol ranks B-cell candidates, using one shared
attribute name instead of a per-tool score column. ScanNet and
DiscoTope-3.0 already comply.
Any future B-cell protocol (e.g. BepiPred/EpiDope ported to Scipion) MUST
also expose ``_meanScore`` to be ranked correctly here -- if it does not,
this protocol falls back to ranking by ROI length instead of silently
misbehaving, but that fallback is a degraded mode, not the intended path.
"""

import re
from typing import Any, List, Optional, Tuple

from ..constants import (
    CONSTRUCT_LINKER_ADJUVANT, CONSTRUCT_LINKER_BCELL, CONSTRUCT_LINKER_CTL,
    CONSTRUCT_LINKER_HTL, CONSTRUCT_LINKER_INTERBLOCK,
)

_GLYCO_ENTRY = re.compile(r'(\d+):Glicosilado')


def getAttr(roi, name: str, default: Any = None) -> Any:
    """Read a Scipion Object-wrapped dynamic attribute (Float/Integer/Boolean/String) as a plain Python value."""
    attr = getattr(roi, name, None)
    return attr.get() if attr is not None else default


def scoreNote(roi, fields: List[str]) -> str:
    parts = []
    for field in fields:
        value = getAttr(roi, field, None)
        if value is not None:
            parts.append(f'{field.lstrip("_")}={value}')
    return ', '.join(parts)


def scoreNoteBcell(roi) -> str:
    return scoreNote(roi, ['_meanScore', '_maxScore', '_algpredVerdict'])


def scoreNoteHtl(roi) -> str:
    return scoreNote(roi, ['_nPromiscuousAlleles', '_nAllelesEvaluated', '_minRankEl'])


def scoreNoteCtl(roi) -> str:
    return scoreNote(roi, ['_nPromiscuousAlleles', '_nAllelesEvaluated', '_minRankEl',
                             '_netcleaveCTermMatch', '_netcleaveCTermScore'])


def selectBcellCandidates(rois: List, topN: int) -> List:
    """Filter by no glycosylated sequon (only if that attribute is present), rank, top-N.

    Does NOT exclude ROIs with ``_algpredVerdict == 'Allergen'``: that
    verdict is not a reliable proxy for real clinical/population
    correlate on its own, and can flag candidates independent of real
    antigen chemistry. Every ROI enters the ranking regardless of
    allergenicity verdict -- ``_algpredVerdict`` stays available on the
    ROI (already surfaced via ``scoreNoteBcell``) for an informed
    decision downstream.

    An upstream ROI that never went through AlgPred2/StackGlyEmbed
    annotation (attribute absent) is NOT excluded on that basis -- absence
    of an attribute is treated as "not evaluated", not "failed".
    """
    candidates = [roi for roi in rois if not getAttr(roi, '_hasGlycoSequon', False)]

    def sortKey(roi):
        score = getAttr(roi, '_meanScore')
        if score is not None:
            return (1, score)
        return (0, len(roi.getROISequence()))

    candidates.sort(key=sortKey, reverse=True)
    return candidates[:topN]


def extractGlycoRegions(bcellRois: List) -> List[Tuple[str, int, int]]:
    """Translate every 'Glicosilado' sequon annotation into an ABSOLUTE (parent_sequence, start, end) region.

    ``_glycoSequonSummary`` (set by the StackGlyEmbed protocol) already
    stores ABSOLUTE positions, not local to the ROI's own window -- see
    that protocol's own docstring. Each region spans the 3-residue
    N-X-[S/T] sequon starting at the annotated Asn position.
    """
    regions = []
    for roi in bcellRois:
        summary = getAttr(roi, '_glycoSequonSummary', '') or ''
        if not summary:
            continue
        parentSeq = roi._sequence.getSequence()
        for match in _GLYCO_ENTRY.finditer(summary):
            pos = int(match.group(1))
            regions.append((parentSeq, pos, pos + 2))
    return regions


def overlapsGlycoRegion(roi, glycoRegions: List[Tuple[str, int, int]]) -> bool:
    """Whether ``roi``'s [start, end] range overlaps any glycosylated region."""
    if not glycoRegions:
        return False
    parentSeq = roi._sequence.getSequence()
    start, end = roi.getROIIdx(), roi.getROIIdx2()
    for seq, gStart, gEnd in glycoRegions:
        if seq == parentSeq and gStart <= end and gEnd >= start:
            return True
    return False


def dedupeByCore(rois: List, sortKeyFn) -> List:
    """Collapse ROIs sharing the same '_core9aa' (same MHC-binding core evaluated in
    neighbouring windows is the same prediction, not distinct epitopes), keeping the
    best-ranked one per ``sortKeyFn`` (ascending)."""
    ordered = sorted(rois, key=sortKeyFn)
    seen = set()
    result = []
    for roi in ordered:
        core = getAttr(roi, '_core9aa')
        if core in seen:
            continue
        seen.add(core)
        result.append(roi)
    return result


def selectHtlCandidates(rois: List, glycoRegions: List[Tuple[str, int, int]], topN: int) -> List:
    """Exclude glycosylated windows, dedupe by core, rank by promiscuity then %Rank, top-N."""
    filtered = [r for r in rois if not overlapsGlycoRegion(r, glycoRegions)]
    if not filtered:
        return filtered

    def keyFn(roi):
        nProm = getAttr(roi, '_nPromiscuousAlleles', 0)
        minRank = getAttr(roi, '_minRankEl', float('inf'))
        return (-nProm, minRank)

    deduped = dedupeByCore(filtered, keyFn)
    deduped.sort(key=keyFn)
    return deduped[:topN]


def selectCtlCandidates(rois: List, glycoRegions: List[Tuple[str, int, int]], topN: int) -> List:
    """Same as HTL, but prioritizes a confirmed NetCleave C-terminal cleavage match first."""
    filtered = [r for r in rois if not overlapsGlycoRegion(r, glycoRegions)]
    if not filtered:
        return filtered

    def keyFn(roi):
        netcleaveMatch = getAttr(roi, '_netcleaveCTermMatch', False)
        nProm = getAttr(roi, '_nPromiscuousAlleles', 0)
        minRank = getAttr(roi, '_minRankEl', float('inf'))
        return (0 if netcleaveMatch else 1, -nProm, minRank)

    deduped = dedupeByCore(filtered, keyFn)
    deduped.sort(key=keyFn)
    return deduped[:topN]


def collectBlocks(bcellSelected: List, htlSelected: List, ctlSelected: List) -> List[Tuple]:
    """Build the ordered (label, rois, intraLinker, seqGetter, noteFn) block list,
    skipping any class with no selected candidates."""
    blockSpecs = [
        ('B-cell', bcellSelected, CONSTRUCT_LINKER_BCELL, lambda r: r.getROISequence(), scoreNoteBcell),
        ('HTL', htlSelected, CONSTRUCT_LINKER_HTL, lambda r: getAttr(r, '_core9aa'), scoreNoteHtl),
        ('CTL', ctlSelected, CONSTRUCT_LINKER_CTL, lambda r: getAttr(r, '_core9aa'), scoreNoteCtl),
    ]
    return [spec for spec in blockSpecs if spec[1]]


def roiSourceFields(sourceRoi) -> dict:
    """Segment metadata sourced from the originating ROI, or all-None for a linker/adjuvant segment."""
    if sourceRoi is None:
        return {'source_parent_id': None, 'source_start': None, 'source_end': None}
    return {
        'source_parent_id': sourceRoi._sequence.getId(),
        'source_start': sourceRoi.getROIIdx(),
        'source_end': sourceRoi.getROIIdx2(),
    }


def appendSegment(segments: List[dict], cursor: int, blockLabel: str, sequence: str,
                    sourceRoi=None, scoreNote: str = '') -> int:
    """Append one segment (epitope or linker) to ``segments`` and return the next cursor position."""
    end = cursor + len(sequence) - 1
    segments.append({
        'block': blockLabel,
        'sequence': sequence,
        'start': cursor,
        'end': end,
        'source_score_note': scoreNote,
        **roiSourceFields(sourceRoi),
    })
    return end + 1


def appendBlockSegments(segments: List[dict], cursor: int, block: Tuple) -> int:
    """Append every ROI in one class block, with its intra-class linker between consecutive ROIs."""
    label, rois, intraLinker, seqGetter, noteFn = block
    lastIdx = len(rois) - 1
    for i, roi in enumerate(rois):
        cursor = appendSegment(segments, cursor, label, seqGetter(roi), roi, noteFn(roi))
        if i < lastIdx:
            cursor = appendSegment(segments, cursor, f'Linker (intra-{label})', intraLinker)
    return cursor


def assembleConstruct(
    bcellSelected: List, htlSelected: List, ctlSelected: List,
    adjuvantSequence: Optional[str] = None,
) -> Tuple[str, List[dict]]:
    """Concatenate the selected candidates with the standard linkers into the final construct.

    Block order: B-cell -> HTL -> CTL (fixed, see module/protocol docstring
    for rationale). No adjuvant by default; if given, it is prepended at
    the N-terminal with its own rigid EAAAK linker.

    Returns:
        ``(construct_sequence, segments)``: the assembled sequence, and a
        list of one dict per SEGMENT (epitope or linker, in order) with
        keys ``block``, ``sequence``, ``start``, ``end`` (1-indexed
        position in the final construct), ``source_parent_id``,
        ``source_start``, ``source_end`` (``None`` for linker/adjuvant
        segments), ``source_score_note``. Concatenating every segment's
        ``sequence`` in order reconstructs ``construct_sequence`` exactly.
        ``("", [])`` if all 3 classes are empty and no adjuvant is given.
    """
    blocks = collectBlocks(bcellSelected, htlSelected, ctlSelected)
    if not blocks and not adjuvantSequence:
        return '', []

    segments: List[dict] = []
    cursor = 1

    if adjuvantSequence:
        cursor = appendSegment(segments, cursor, 'Adjuvant', adjuvantSequence)
        cursor = appendSegment(segments, cursor, 'Linker', CONSTRUCT_LINKER_ADJUVANT)

    lastBlockIdx = len(blocks) - 1
    for blockIdx, block in enumerate(blocks):
        cursor = appendBlockSegments(segments, cursor, block)
        if blockIdx < lastBlockIdx:
            cursor = appendSegment(segments, cursor, 'Linker (inter-block)', CONSTRUCT_LINKER_INTERBLOCK)

    constructSequence = ''.join(s['sequence'] for s in segments)
    return constructSequence, segments
