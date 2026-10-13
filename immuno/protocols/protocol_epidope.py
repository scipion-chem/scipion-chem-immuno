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
This protocol is used to predict linear (sequence-based) B-cell epitopes
from a protein sequence with a local EpiDope installation.
"""

from pathlib import Path

import pandas as pd
from pwchem.objects import Sequence, SequenceROI, SetOfSequenceROIs
from pwem.protocols import EMProtocol
from pyworkflow.object import Float
from pyworkflow.protocol import params
from pyworkflow.utils import Message

from immuno import Plugin as immunoPlugin
from ..constants import (
    EPIDOPE_DEFAULT_THRESHOLD as DEFAULT_THRESHOLD,
    EPIDOPE_RESIDUE_COLUMN as RESIDUE_COLUMN,
    EPIDOPE_SCORE_COLUMN as SCORE_COLUMN,
)
from ..utils.epidope_utils import loadRawScores
from ..utils.sliding_window_mapping import extractLinearEpitopeRegions


class ProtEpiDopePrediction(EMProtocol):
    """
    AI Generated:

    Predicts linear (sequence-based) B-cell epitopes along a protein
    sequence using a local EpiDope installation (legacy conda-env runtime,
    see plugin ``README.rst``), and maps the resulting per-residue scores
    into contiguous epitope regions via a gap-tolerant sliding window (the
    same algorithm this plugin applies to DiscoTope-3.0 and ScanNet).

    EpiDope scores every residue independently, so a single weak residue
    inside an otherwise strong stretch would break the region in two: the
    window therefore tolerates up to 'maxGapResidues' below-threshold
    residues before it stops extending the current candidate.

    Output
    ------
    outputROIs: SetOfSequenceROIs, one ROI per epitope region found, each
    annotated with its mean EpiDope score (_meanScore).
    """

    _label = 'epidope epitope prediction'

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSequence', params.PointerParam, pointerClass='Sequence',
                      label='Input protein sequence: ',
                      help='Protein sequence to scan for linear B-cell epitopes.')

        gGroup = form.addGroup('Epitope mapping')
        gGroup.addParam('threshold', params.FloatParam, default=DEFAULT_THRESHOLD,
                        label='Score threshold: ',
                        help='EpiDope score above which a residue is considered part of a '
                             'candidate epitope.')
        gGroup.addParam('windowSize', params.IntParam, default=9, label='Window size (aa): ',
                        help='Length (in residues) of the sliding window used to collapse '
                             'above-threshold residues into contiguous candidate epitope regions. '
                             'The default is the minimum B-cell recognition footprint.')
        gGroup.addParam('maxGapResidues', params.IntParam, default=2,
                        label='Max. below-threshold residues per window: ',
                        help='How many below-threshold residues are tolerated inside a sliding '
                             'window before it stops extending the current candidate region.')
        gGroup.addParam('minLength', params.IntParam, default=9, label='Min. epitope length (aa): ',
                        help='Minimum length (in residues) a mapped region must reach to be kept '
                             'as a candidate epitope; shorter regions are discarded.')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.epidopeStep)
        self._insertFunctionStep(self.createOutputStep)

    # ---------------------------------- Steps -----------------------------------

    def epidopeStep(self):
        fastaPath = Path(self._getExtraPath('inputSequence.fa'))
        self.inputSequence.get().exportToFile(str(fastaPath))

        resultDir = Path(self._getExtraPath('epidope_out'))
        resultDir.mkdir(parents=True, exist_ok=True)

        args = (
            f'-i {fastaPath.resolve()} '
            f'-o {resultDir.resolve()} '
            f'-t {self.threshold.get()}'
        )
        immunoPlugin.runEpiDope(self, args)

        loadRawScores(resultDir).to_csv(self._getRawScoresPath(), index=False)

    def createOutputStep(self):
        rawDf = pd.read_csv(self._getRawScoresPath())
        if rawDf.empty:
            return

        residues = rawDf[RESIDUE_COLUMN].astype(str).tolist()
        regionsDf = extractLinearEpitopeRegions(
            scores=rawDf[SCORE_COLUMN].tolist(), residues=residues,
            threshold=self.threshold.get(), minLength=self.minLength.get(),
            windowSize=self.windowSize.get(), maxGapResidues=self.maxGapResidues.get(),
        )
        if regionsDf.empty:
            return

        inputSeq = self.inputSequence.get()
        outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
        for row in regionsDf.itertuples(index=False):
            roiId = f'ROI_{row.start}-{row.end}'
            roiSeq = Sequence(sequence=row.sequence, name=roiId, id=roiId,
                              description='EpiDope linear epitope region')
            newRoi = SequenceROI(sequence=inputSeq, seqROI=roiSeq, roiIdx=row.start, roiIdx2=row.end)
            newRoi._meanScore = Float(row.mean_score)
            outROIs.append(newRoi)

        self._defineOutputs(outputROIs=outROIs)
        self._defineSourceRelation(self.inputSequence, outROIs)

    # ---------------------------------- Utils -----------------------------------

    def _getRawScoresPath(self):
        return self._getExtraPath('epidope_scores.csv')

    # ---------------------------------- Info ------------------------------------

    def _summary(self):
        summary = []
        if self.isFinished():
            outROIs = getattr(self, 'outputROIs', None)
            nROIs = len(outROIs) if outROIs is not None else 0
            summary.append(f'{nROIs} epitope region(s) found above threshold {self.threshold.get()}.')
        return summary

    def _validate(self):
        return immunoPlugin.validateEpiDopeInstallation()
