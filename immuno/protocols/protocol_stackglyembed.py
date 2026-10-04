# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Enzo Sierra (enzogael57@gmail.com)
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
This protocol is used to predict N-linked glycosylation sequons of a set
of peptide candidates with a local StackGlyEmbed installation.
"""

import os

from pwchem.objects import SetOfSequenceROIs
from pwem.protocols import EMProtocol
from pyworkflow.object import Boolean, String
from pyworkflow.protocol import params

from immuno import Plugin as immunoPlugin
from ..constants import STACKGLYEMBED_DIC
from ..utils.stackglyembed_utils import buildDataset, parsePredictions, scanSequons


class ProtStackGlyEmbedPrediction(EMProtocol):
    """
    AI Generated:

    Scans every input ROI's peptide for N-X-[S/T] N-glycosylation sequons
    and predicts each one with a local StackGlyEmbed installation, then
    annotates (does NOT filter) each ROI with the result.

    Overview
    --------
    1. For each input ROI, scan its own peptide substring for sequons
       (N-X-[S/T], Xaa != Pro), and map each hit to its ABSOLUTE position on
       the ROI's parent sequence.
    2. Group sites by unique parent sequence: StackGlyEmbed's embeddings
       (ProteinBERT + ESM-2 + ProtT5) are computed on the FULL protein, not
       an isolated fragment, so several ROIs sharing the same parent only
       trigger one embedding pass for that parent, not one per ROI.
    3. Run StackGlyEmbed once on the deduplicated parent/site list
       (stackglyembedStep), persisting the raw predictions.
    4. createOutputStep re-attaches results per ROI: ``_hasGlycoSequon``
       (True if >=1 of its own sequons predicts 'Glycosylated') and
       ``_glycoSequonSummary`` (human-readable detail, one entry per sequon
       found in that ROI, empty string if the ROI has no sequon at all).

    Output
    ------
    outputROIs: the same SetOfSequenceROIs as the input, annotated. Does not
    filter -- e.g. the construct assembly protocol decides to exclude
    B-cell candidates with >=1 'Glycosylated' sequon, but that decision does
    not belong here.
    """

    _label = 'stackglyembed n-glyco prediction'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputROIs', params.PointerParam, pointerClass='SetOfSequenceROIs',
                       label='Sequence ROIs: ',
                       help='Peptide candidates to scan for N-glycosylation sequons.')
        form.addParam('timeoutSeconds', params.IntParam, label='Timeout (s): ', default=900,
                       expertLevel=params.LEVEL_ADVANCED,
                       help='StackGlyEmbed loads 3 protein language models (ProteinBERT, ESM-2 650M, '
                            'ProtT5) on a cold start, so a generous timeout is needed.')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.stackglyembedStep)
        self._insertFunctionStep(self.createOutputStep)

    # ---------------------------------- Steps -----------------------------------

    def _getRois(self):
        # Iterating a Scipion SetOfXXX reuses the same Python object per row
        # (the underlying sqlite cursor): each item must be cloned when
        # materialized into a list, or all N references end up pointing to
        # the cursor's last state.
        return [roi.clone() for roi in self.inputROIs.get()]

    def _getRoiSequons(self, roi):
        """Returns [(absolute_site_position, ...)] for the sequons found inside this ROI's own peptide."""
        parentSeq = roi._sequence.getSequence()
        roiStart = roi.getROIIdx()
        localHits = scanSequons(roi.getROISequence())
        return parentSeq, [roiStart + local - 1 for local in localHits]

    def _getPredictedPath(self):
        return self._getExtraPath('predicted_values.csv')

    def stackglyembedStep(self):
        rois = self._getRois()

        parentSites = {}
        for roi in rois:
            parentSeq, sites = self._getRoiSequons(roi)
            if not sites:
                continue
            parentSites.setdefault(parentSeq, set()).update(sites)

        if not parentSites:
            return

        parentSites = {seq: sorted(sites) for seq, sites in parentSites.items()}

        datasetPath = self._getExtraPath('dataset.txt')
        order = buildDataset(parentSites, datasetPath)
        # Persist evaluation order alongside the raw predictions: createOutputStep
        # runs in a SEPARATE subprocess, so it must know the exact row order
        # predict_local.py produced to zip results back to the right (parent,
        # site) pair. The full sequence string is written verbatim (not a
        # Python hash()): CPython's built-in hash() is randomized per-process
        # by default (PYTHONHASHSEED), so a hash computed here would not
        # match one recomputed in the next step's own subprocess.
        with open(self._getExtraPath('order.txt'), 'w') as fh:
            for seq, pos in order:
                fh.write(f'{pos}\t{seq}\n')

        args = (
            f'--dataset {datasetPath} --output-dir {self._getExtraPath()} '
            f'--models-dir {immunoPlugin.getStackGlyEmbedModelsDir()} '
            f'--t5-model-path {immunoPlugin.getStackGlyEmbedT5ModelName()} '
            f'--esm-model-name {immunoPlugin.getVar(STACKGLYEMBED_DIC["esm_model_name"])} '
            f'--proteinbert-dir {immunoPlugin.getStackGlyEmbedProteinBertDir()}'
        )
        scriptPath = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'scripts', 'stackglyembed_predict_local.py')
        immunoPlugin.runStackGlyEmbed(self, scriptPath, args)

    def createOutputStep(self):
        rois = self._getRois()

        predictedPath = self._getPredictedPath()
        if not os.path.isfile(predictedPath):
            # No ROI had any sequon: stackglyembedStep never ran the program.
            resultsBySite = {}
        else:
            order = []
            with open(self._getExtraPath('order.txt')) as fh:
                for line in fh:
                    pos, seq = line.rstrip('\n').split('\t', 1)
                    order.append((seq, int(pos)))
            resultDf = parsePredictions(predictedPath, order)
            resultsBySite = {
                (row.sequence, row.site_position): (row.stackglyembed_verdict, row.stackglyembed_score)
                for row in resultDf.itertuples(index=False)
            }

        outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
        for roi in rois:
            parentSeq, sites = self._getRoiSequons(roi)
            details = []
            hasGlyco = False
            for pos in sites:
                verdict, score = resultsBySite.get((parentSeq, pos), ('Not evaluated', float('nan')))
                if verdict == 'Glycosylated':
                    hasGlyco = True
                details.append(f'{pos}:{verdict}({score:.3f})' if score == score else f'{pos}:{verdict}')

            roi._hasGlycoSequon = Boolean(hasGlyco)
            roi._glycoSequonSummary = String(';'.join(details))
            outROIs.append(roi)

        if len(outROIs) > 0:
            self._defineOutputs(outputROIs=outROIs)
            self._defineSourceRelation(self.inputROIs, outROIs)

    # ---------------------------------- Validation -------------------------------

    def _validate(self):
        return immunoPlugin.validateStackGlyEmbedInstallation()

    def _summary(self):
        summary = []
        if self.isFinished():
            outROIs = getattr(self, 'outputROIs', None)
            if outROIs is not None:
                nGlyco = sum(1 for roi in outROIs if roi._hasGlycoSequon.get())
                summary.append(f'{nGlyco}/{len(outROIs)} candidate(s) with >= 1 glycosylated sequon.')
        return summary
