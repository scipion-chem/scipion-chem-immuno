# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Enzo Sierra (enzogael57@gmail.com)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

"""
This protocol is used to automatically assemble a multi-epitope vaccine
construct from B-cell/HTL/CTL candidates already selected by upstream
protocols.
"""

import csv
import os

from pwchem.objects import Sequence, SequenceROI, SetOfSequenceROIs
from pwem.protocols import EMProtocol
from pyworkflow.protocol import params

from ..constants import CONSTRUCT_DEFAULT_TOP_N_PER_CLASS
from ..utils.epitopeconstruct_assembly import (
    assembleConstruct, extractGlycoRegions, selectBcellCandidates,
    selectCtlCandidates, selectHtlCandidates,
)
from ..utils.epitopeconstruct_exceptions import EpitopeConstructError


class ProtEpitopeConstructAssembly(EMProtocol):
    """
    AI Generated:

    Automatically assembles a multi-epitope vaccine construct from
    already-selected B-cell/HTL/CTL candidates, concatenated with the
    field-standard linkers for multi-epitope vaccine design (multiple
    independent sources agree, not a fixed biological rule).

    Wraps NO external tool: pure selection/assembly logic. Kept as a
    separate protocol from ProtOptimizeMultiEpitope (pwchem core, uses
    genetic algorithms) rather than merged into one: the assembly method
    itself (chimeric, non-GA) genuinely differs, not just its packaging.

    Selection per class (top-N each, default 3 -- a real run easily
    produces 15-20+ valid HTL/CTL candidates alone, too many for a
    manageable construct):

    - **B-cell**: excludes any candidate with >=1 'Glicosilado' sequon
      (StackGlyEmbed) -- ONLY if that upstream annotation is actually
      present on the input ROIs (its absence is treated as "not
      evaluated", never as "failed"). Does NOT exclude 'Allergen'
      (AlgPred2) candidates: that verdict often has no real clinical
      correlate on its own; every ROI is ranked regardless of
      allergenicity, with ``_algpredVerdict`` kept visible on the ROI for
      an informed decision. Ranked by ``_meanScore``, a PROJECT-WIDE
      CONVENTION:
      every B-cell prediction protocol in this project MUST expose
      ``_meanScore`` (``Float``) on its output ROIs to be ranked correctly
      here -- ScanNet and DiscoTope-3.0 already comply. ROI length is a
      fallback for protocols that don't (a degraded mode, not the
      intended path).
    - **HTL/CTL**: collapsed by ``_core9aa`` (the same MHC-binding core
      evaluated in neighbouring windows is the same prediction, not
      distinct epitopes), excluding any window overlapping a glycosylated
      sequon from the B-cell input (glycosylation can physically block
      MHC groove binding regardless of presentation pathway, unlike
      allergenicity -- see below). Ranked by promiscuous-allele count then
      %Rank_EL; CTL additionally prioritizes a confirmed NetCleave
      C-terminal cleavage match first.
    - Allergenicity (AlgPred2) is deliberately NOT applied to HTL/CTL:
      IgE/mast-cell recognition needs the epitope to circulate intact and
      exposed -- never true for an 8-11 aa core buried in the MHC groove.

    Design decisions kept fixed (no strong literature consensus otherwise,
    or explicit project decision): block order B-cell -> HTL -> CTL;
    ``core_9aa`` (not the full evaluated window) is what gets inserted for
    HTL/CTL; no fusion of overlapping candidates from DIFFERENT classes
    (would break linker semantics -- intra-class overlap is already
    resolved upstream, before classes separate); no adjuvant by default
    (``adjuvantSequence`` param available to add one with its own rigid
    EAAAK linker, Arai et al. 2001, without redesigning anything).

    Output
    ------
    outputROIs: SetOfSequenceROIs with EXACTLY ONE SequenceROI -- the full
    assembled construct (its own parent sequence, spanning the whole
    length) -- so this same output can be fed directly as ``inputROIs``
    into a construct-level check protocol (AlgPred2, ToxinPred2, IApred,
    SignalP-6.0), which all expect a SetOfSequenceROIs. The full
    per-segment breakdown (which peptide/linker contributed each stretch,
    in what order, from where) is persisted to
    ``extra/construct_metadata.csv``.
    """

    _label = 'epitope construct assembly'
    _CONSTRUCT_FASTA_FILENAME = 'construct.fasta'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('bcellROIs', params.PointerParam, pointerClass='SetOfSequenceROIs',
                       label='B-cell candidates: ', allowsNull=True,
                       help='Output of an upstream B-cell prediction protocol (BepiPred/EpiDope/'
                            'ScanNet/DiscoTope-3.0), ideally already annotated by AlgPred2 and '
                            'StackGlyEmbed (both annotations are optional, but recommended).')
        form.addParam('htlROIs', params.PointerParam, pointerClass='SetOfSequenceROIs',
                       label='HTL (MHC-II) candidates: ', allowsNull=True,
                       help='Output of scipion-chem-netmhciipan (must carry _core9aa).')
        form.addParam('ctlROIs', params.PointerParam, pointerClass='SetOfSequenceROIs',
                       label='CTL (MHC-I) candidates: ', allowsNull=True,
                       help='Output of scipion-chem-netmhcpan, ideally further annotated by '
                            'scipion-chem-netcleave (must carry _core9aa).')
        form.addParam('topNPerClass', params.IntParam, default=CONSTRUCT_DEFAULT_TOP_N_PER_CLASS,
                       label='Max. epitopes per class: ')
        form.addParam('adjuvantSequence', params.StringParam, default='',
                       label='Adjuvant sequence (optional): ',
                       help='If given, prepended at the N-terminal with its own rigid EAAAK linker '
                            '(Arai et al. 2001). Empty by default -- no adjuvant is chosen '
                            'automatically, that choice needs pathogen/host-specific biological '
                            'judgement outside this protocol\'s scope.')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.assembleStep)
        self._insertFunctionStep(self.createOutputStep)

    # ---------------------------------- Steps -----------------------------------

    def _getRois(self, pointerParam):
        pointer = getattr(self, pointerParam)
        if pointer.get() is None:
            return []
        # Iterating a Scipion SetOfXXX reuses the same Python object per
        # row (the underlying sqlite cursor): each item must be cloned when
        # materialized into a list, or all N references end up pointing to
        # the cursor's last state.
        return [roi.clone() for roi in pointer.get()]

    def assembleStep(self):
        bcellRois = self._getRois('bcellROIs')
        htlRois = self._getRois('htlROIs')
        ctlRois = self._getRois('ctlROIs')

        topN = self.topNPerClass.get()
        glycoRegions = extractGlycoRegions(bcellRois)

        bcellSelected = selectBcellCandidates(bcellRois, topN)
        htlSelected = selectHtlCandidates(htlRois, glycoRegions, topN)
        ctlSelected = selectCtlCandidates(ctlRois, glycoRegions, topN)

        adjuvant = self.adjuvantSequence.get().strip() or None
        constructSequence, segments = assembleConstruct(bcellSelected, htlSelected, ctlSelected, adjuvant)

        if not constructSequence:
            raise EpitopeConstructError(
                'No candidate survived selection in any class (B-cell/HTL/CTL) and no adjuvant '
                'sequence was given: there is nothing to assemble a construct from.'
            )

        with open(self._getExtraPath('construct_metadata.csv'), 'w', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=[
                'block', 'sequence', 'start', 'end', 'source_parent_id', 'source_start',
                'source_end', 'source_score_note',
            ])
            writer.writeheader()
            writer.writerows(segments)

        with open(self._getExtraPath(self._CONSTRUCT_FASTA_FILENAME), 'w') as fh:
            fh.write(f'>construct\n{constructSequence}\n')

    def createOutputStep(self):
        fastaPath = self._getExtraPath(self._CONSTRUCT_FASTA_FILENAME)
        if not os.path.isfile(fastaPath):
            return

        with open(fastaPath) as fh:
            fh.readline()
            constructSequence = fh.readline().strip()
        if not constructSequence:
            return

        # Single-ROI output spanning the whole construct: this is exactly
        # the shape a construct-level check protocol (AlgPred2, ToxinPred2,
        # IApred, SignalP-6.0) expects as 'inputROIs'.
        parentSeq = Sequence(sequence=constructSequence, name='construct', id='construct',
                              description='Assembled multi-epitope construct')
        roiSeq = Sequence(sequence=constructSequence, name='construct_roi', id='construct_roi',
                           description='Assembled multi-epitope construct')
        constructRoi = SequenceROI(sequence=parentSeq, seqROI=roiSeq, roiIdx=1, roiIdx2=len(constructSequence))

        outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
        outROIs.append(constructRoi)
        self._defineOutputs(outputROIs=outROIs)

        for pointerParam in ('bcellROIs', 'htlROIs', 'ctlROIs'):
            pointer = getattr(self, pointerParam)
            if pointer.get() is not None:
                self._defineSourceRelation(pointer, outROIs)

    # ---------------------------------- Validation -------------------------------

    def _validate(self):
        errors = []
        if self.bcellROIs.get() is None and self.htlROIs.get() is None and self.ctlROIs.get() is None:
            errors.append('At least one of B-cell/HTL/CTL candidate sets must be provided.')
        return errors

    def _summary(self):
        summary = []
        fastaPath = self._getExtraPath(self._CONSTRUCT_FASTA_FILENAME)
        if self.isFinished() and os.path.isfile(fastaPath):
            with open(fastaPath) as fh:
                fh.readline()
                seq = fh.readline().strip()
            summary.append(f'Construct length: {len(seq)} aa.')
        return summary
