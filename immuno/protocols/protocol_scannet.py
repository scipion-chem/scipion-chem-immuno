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
This protocol is used to predict structure-based (conformational) B-cell
epitopes from a PDB structure with a local ScanNet installation.
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd
from pwchem.objects import Sequence, SequenceROI, SetOfSequenceROIs
from pwem.protocols import EMProtocol
from pyworkflow.object import Float
from pyworkflow.protocol import params

from immuno import Plugin as immunoPlugin
from ..constants import (
    SCANNET_RAW_CHAIN_COLUMN as RAW_CHAIN_COLUMN,
    SCANNET_RAW_RESIDUE_COLUMN as RAW_RESIDUE_COLUMN,
    SCANNET_RAW_SCORE_COLUMN as RAW_SCORE_COLUMN,
)
from ..utils.scannet_epitope_mapping import extractEpitopeRegions
from ..utils.scannet_exceptions import ScanNetExecutionError
from ..utils.scannet_utils import buildCommand, loadRawScores


class ProtScanNetPrediction(EMProtocol):
    """
    AI Generated:

    Predicts structure-based (conformational) B-cell epitopes from a PDB
    structure using a local ScanNet installation (legacy conda-env
    runtime, see plugin ``README.rst``), and maps the resulting per-residue
    scores into contiguous linear epitope regions via a gap-tolerant
    sliding window (same algorithm this project already uses for BepiPred/
    EpiDope/DiscoTope-3.0).

    ScanNet has no absolute score threshold published by its authors
    (unlike DiscoTope-3.0's 'calibrated_score'): by default the threshold
    is computed ADAPTIVELY, separately for EACH chain, as a percentile of
    that chain's own scores -- the same per-antigen normalization principle
    DiscoTope-3.0 applies internally, applied here on this side since
    ScanNet does not do it on its own.

    Biological note: ScanNet scores CONFORMATIONAL epitopes (potentially
    sequence-discontinuous 3D patches); collapsing to contiguous linear
    regions via sliding window is a deliberate simplification to keep
    downstream phases operating on synthesizable linear peptides.

    Output
    ------
    outputROIs: SetOfSequenceROIs, one ROI per mapped epitope region
    (possibly several chains, one Sequence object per chain built from
    ScanNet's own per-residue output).
    """

    _label = 'scannet epitope prediction'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStructure', params.PointerParam, pointerClass='AtomStruct',
                       label='Input structure: ',
                       help='PDB structure to scan for conformational B-cell epitopes.')
        form.addParam('timeoutSeconds', params.IntParam, label='Timeout (s): ', default=1800,
                       expertLevel=params.LEVEL_ADVANCED,
                       help='Maximum time the ScanNet run is allowed to take before the step is '
                            'aborted as failed.')

        gGroup = form.addGroup('Epitope mapping')
        gGroup.addParam('thresholdPercentile', params.FloatParam, default=90.0,
                         label='Adaptive threshold percentile: ',
                         help='Per-chain percentile (0-100) of ScanNet scores used as the adaptive '
                              'threshold (default 90: the top 10%% of each chain\'s residues count as '
                              '"high").')
        gGroup.addParam('windowSize', params.IntParam, default=9, label='Window size (aa): ',
                         help='Length (in residues) of the sliding window used to collapse '
                              'above-threshold residues into contiguous candidate epitope regions.')
        gGroup.addParam('maxGapResidues', params.IntParam, default=2,
                         label='Max. below-threshold residues per window: ',
                         help='How many below-threshold residues are tolerated inside a sliding '
                              'window before it stops extending the current candidate region.')
        gGroup.addParam('minLength', params.IntParam, default=9, label='Min. epitope length (aa): ',
                         help='Minimum length (in residues) a mapped region must reach to be kept '
                              'as a candidate epitope; shorter regions are discarded.')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.scannetStep)
        self._insertFunctionStep(self.createOutputStep)

    # ---------------------------------- Steps -----------------------------------

    def _getRawScoresPath(self):
        return self._getExtraPath('raw_scores.csv')

    def scannetStep(self):
        pdbPath = Path(self.inputStructure.get().getFileName())
        resultDir = Path(self._getExtraPath('scannet_out'))
        resultDir.mkdir(parents=True, exist_ok=True)
        name = pdbPath.stem

        home = immunoPlugin.getScanNetDir()
        args = buildCommand(pdbPath, resultDir, name)
        immunoPlugin.runScanNet(self, args, cwd=home)

        rawDf = loadRawScores(resultDir)
        rawDf.to_csv(self._getRawScoresPath(), index=False)

    def createOutputStep(self):
        if not os.path.isfile(self._getRawScoresPath()):
            return

        rawDf = pd.read_csv(self._getRawScoresPath())
        if rawDf.empty:
            return

        frames = []
        for chain, group in rawDf.groupby(RAW_CHAIN_COLUMN, sort=False):
            adaptiveThreshold = float(np.percentile(group[RAW_SCORE_COLUMN], self.thresholdPercentile.get()))
            regionsDf = extractEpitopeRegions(
                group, groupCol=RAW_CHAIN_COLUMN, scoreCol=RAW_SCORE_COLUMN, residueCol=RAW_RESIDUE_COLUMN,
                threshold=adaptiveThreshold, minLength=self.minLength.get(),
                windowSize=self.windowSize.get(), maxGapResidues=self.maxGapResidues.get(),
            )
            frames.append(regionsDf)
        epitopesDf = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

        if epitopesDf.empty:
            return

        # One parent Sequence per chain, built from ScanNet's own
        # per-residue output (not re-derived from the PDB) -- guarantees
        # the ROI substring matches exactly what ScanNet actually scored.
        chainSequences = {
            chain: ''.join(group[RAW_RESIDUE_COLUMN].astype(str))
            for chain, group in rawDf.groupby(RAW_CHAIN_COLUMN, sort=False)
        }

        outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
        for row in epitopesDf.itertuples(index=False):
            chainSeq = chainSequences[row.group]
            parentSeq = Sequence(sequence=chainSeq, name=f'chain_{row.group}', id=f'chain_{row.group}',
                                  description=f'ScanNet input, chain {row.group}')
            roiId = f'ROI_{row.group}_{row.start}-{row.end}'
            roiSeq = Sequence(sequence=row.sequence, name=roiId, id=roiId,
                               description='ScanNet conformational epitope region')
            newRoi = SequenceROI(sequence=parentSeq, seqROI=roiSeq, roiIdx=row.start, roiIdx2=row.end)
            newRoi._meanScore = Float(row.mean_score)
            newRoi._maxScore = Float(row.max_score)
            outROIs.append(newRoi)

        if len(outROIs) > 0:
            self._defineOutputs(outputROIs=outROIs)
            self._defineSourceRelation(self.inputStructure, outROIs)

    # ---------------------------------- Validation -------------------------------

    def _validate(self):
        return immunoPlugin.validateScanNetInstallation()

    def _summary(self):
        summary = []
        if self.isFinished():
            outROIs = getattr(self, 'outputROIs', None)
            n = len(outROIs) if outROIs is not None else 0
            summary.append(f'{n} conformational epitope region(s) found.')
        return summary
