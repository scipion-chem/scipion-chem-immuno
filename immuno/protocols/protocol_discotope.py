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
epitopes from a SINGLE-CHAIN PDB structure with a local DiscoTope-3.0
installation.
"""

import glob
import os
from pathlib import Path

import pandas as pd
from pwchem.objects import Sequence, SequenceROI, SetOfSequenceROIs
from pwem.protocols import EMProtocol
from pyworkflow.object import Float
from pyworkflow.protocol import params

from immuno import Plugin as immunoPlugin
from ..constants import (
    DISCOTOPE_DEFAULT_THRESHOLD as DEFAULT_THRESHOLD,
    DISCOTOPE_RAW_RESIDUE_COLUMN as RAW_RESIDUE_COLUMN,
    DISCOTOPE_RAW_SCORE_COLUMN as RAW_SCORE_COLUMN,
)
from ..utils.discotope_epitope_mapping import extractEpitopeRegions
from ..utils.discotope_exceptions import DiscoTopeExecutionError


class ProtDiscoTopePrediction(EMProtocol):
    """
    AI Generated:

    Predicts structure-based (conformational) B-cell epitopes from a
    SINGLE-CHAIN PDB structure using a local DiscoTope-3.0 installation,
    and maps the resulting per-residue 'calibrated_score' into contiguous
    linear epitope regions via a gap-tolerant sliding window (same
    algorithm this project already uses for BepiPred/EpiDope/ScanNet).

    Unlike ScanNet, DiscoTope-3.0 expects (and only supports) a
    single-chain PDB per run -- a multi-chain input would produce more
    than one output CSV with no reliable way for this protocol to know
    which corresponds to which chain, so it is rejected up front.

    Uses 'calibrated_score' (not the raw 'DiscoTope-3.0_score'):
    DiscoTope-3.0's own authors publish reference thresholds with expected
    recall for this normalized column (see plugin ``constants.py``),
    unlike the raw score, which would need a threshold hand-calibrated
    with no such backing.

    Biological note: DiscoTope-3.0 scores CONFORMATIONAL epitopes
    (potentially sequence-discontinuous 3D patches); collapsing to
    contiguous linear regions via sliding window is a deliberate
    simplification to keep downstream phases operating on synthesizable
    linear peptides.

    Output
    ------
    outputROIs: SetOfSequenceROIs, one ROI per mapped epitope region.
    """

    _label = 'discotope epitope prediction'

    def _defineParams(self, form):
        form.addHidden(params.USE_GPU, params.BooleanParam, default=True,
                       label='Use GPU: ',
                       help='Whether to use GPU or not. (Unable to choose the GPU id).')
        form.addHidden(params.GPU_LIST, params.StringParam, default='0', label='Choose GPU IDs',
                       help='Add a list of GPU devices that can be used')

        form.addSection(label='Input')
        form.addParam('inputStructure', params.PointerParam, pointerClass='AtomStruct',
                       label='Input structure (single chain): ',
                       help='Single-chain PDB structure to scan for conformational B-cell epitopes.')
        form.addParam('strucType', params.EnumParam, choices=['solved', 'alphafold'], default=0,
                       label='Structure type: ',
                       help='"solved" for experimentally solved structures, "alphafold" for '
                            'AlphaFold-predicted models (affects internal pLDDT handling).')
        form.addParam('timeoutSeconds', params.IntParam, label='Timeout (s): ', default=1800,
                       expertLevel=params.LEVEL_ADVANCED,
                       help='Maximum time the DiscoTope-3.0 run is allowed to take before the step '
                            'is aborted as failed.')

        gGroup = form.addGroup('Epitope mapping')
        gGroup.addParam('threshold', params.FloatParam, default=DEFAULT_THRESHOLD,
                         label='calibrated_score threshold: ',
                         help='DiscoTope-3.0\'s own published reference: ~0.40 ("low", ~70%% recall), '
                              '~0.90 ("moderate", default), ~1.51 ("higher", lower recall/higher '
                              'precision).')
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
        self._insertFunctionStep(self.discotopeStep)
        self._insertFunctionStep(self.createOutputStep)

    # ---------------------------------- Steps -----------------------------------

    def _getRawScoresPath(self):
        return self._getExtraPath('raw_scores.csv')

    def discotopeStep(self):
        pdbPath = Path(self.inputStructure.get().getFileName())
        resultDir = Path(self._getExtraPath('discotope_out'))
        resultDir.mkdir(parents=True, exist_ok=True)

        strucType = self.getEnumText('strucType')
        args = (
            f'--pdb_or_zip_file {pdbPath.resolve()} '
            f'--out_dir {resultDir.resolve()} '
            f'--struc_type {strucType} '
            f'--models_dir {immunoPlugin.getDiscoTopeModelsDir()}'
        )
        if not self.useGpu.get():
            args += ' --cpu_only'
        immunoPlugin.runDiscoTope(self, args)

        csvFiles = sorted(Path(p) for p in glob.glob(str(resultDir / '**' / '*_discotope3.csv'), recursive=True))
        if not csvFiles:
            raise DiscoTopeExecutionError(
                f"No '*_discotope3.csv' output file found in '{resultDir}'."
            )
        if len(csvFiles) > 1:
            raise DiscoTopeExecutionError(
                f"More than one output CSV found in '{resultDir}' ({[c.name for c in csvFiles]}): the "
                "input PDB must contain a single chain."
            )

        rawDf = pd.read_csv(csvFiles[0])
        missing = {RAW_RESIDUE_COLUMN, RAW_SCORE_COLUMN} - set(rawDf.columns)
        if missing:
            raise DiscoTopeExecutionError(
                f"Output CSV '{csvFiles[0]}' does not contain the expected columns {sorted(missing)}. "
                f"Columns found: {list(rawDf.columns)}."
            )
        rawDf[[RAW_RESIDUE_COLUMN, RAW_SCORE_COLUMN]].to_csv(self._getRawScoresPath(), index=False)

    def createOutputStep(self):
        if not os.path.isfile(self._getRawScoresPath()):
            return

        rawDf = pd.read_csv(self._getRawScoresPath())
        if rawDf.empty:
            return

        chainSeq = ''.join(rawDf[RAW_RESIDUE_COLUMN].astype(str))
        regionsDf = extractEpitopeRegions(
            scores=rawDf[RAW_SCORE_COLUMN].tolist(), residues=rawDf[RAW_RESIDUE_COLUMN].astype(str).tolist(),
            threshold=self.threshold.get(), minLength=self.minLength.get(),
            windowSize=self.windowSize.get(), maxGapResidues=self.maxGapResidues.get(),
        )
        if regionsDf.empty:
            return

        stem = Path(self.inputStructure.get().getFileName()).stem
        parentSeq = Sequence(sequence=chainSeq, name=stem, id=stem,
                              description='DiscoTope-3.0 input structure')

        outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
        for row in regionsDf.itertuples(index=False):
            roiId = f'ROI_{row.start}-{row.end}'
            roiSeq = Sequence(sequence=row.sequence, name=roiId, id=roiId,
                               description='DiscoTope-3.0 conformational epitope region')
            newRoi = SequenceROI(sequence=parentSeq, seqROI=roiSeq, roiIdx=row.start, roiIdx2=row.end)
            newRoi._meanScore = Float(row.mean_score)
            newRoi._maxScore = Float(row.max_score)
            outROIs.append(newRoi)

        if len(outROIs) > 0:
            self._defineOutputs(outputROIs=outROIs)
            self._defineSourceRelation(self.inputStructure, outROIs)

    # ---------------------------------- Validation -------------------------------

    def _validate(self):
        return immunoPlugin.validateDiscoTopeInstallation()

    def _summary(self):
        summary = []
        if self.isFinished():
            outROIs = getattr(self, 'outputROIs', None)
            n = len(outROIs) if outROIs is not None else 0
            summary.append(f'{n} conformational epitope region(s) found.')
        return summary
