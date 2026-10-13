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
This protocol runs TMbed on a protein sequence to predict transmembrane
helices/strands and the signal peptide, for masking those regions out of
downstream B-cell epitope candidates: transmembrane residues are not
solvent-accessible and signal peptides are cleaved off the mature protein,
so antigenicity predictions falling there are not biologically meaningful.
"""

import os

from pwchem.objects import Sequence, SequenceROI, SetOfSequenceROIs
from pwem.protocols import EMProtocol
from pyworkflow.protocol import params

from immuno import Plugin as immunoPlugin
from ..constants import TMBED_OUT_FORMAT
from ..utils.tmbed_exceptions import TMbedExecutionError
from ..utils.tmbed_utils import extractMaskingRegions, parsePredictions


class ProtTMbedPredict(EMProtocol):
    """
    AI Generated:

    Predicts transmembrane beta-strands/alpha-helices and the signal peptide
    of a protein sequence with a local TMbed installation, and reports each
    contiguous region as a masking ROI (type 'TM_beta_strand',
    'TM_alpha_helix' or 'signal_peptide'). Non-membrane residues ('i'/'o' in
    TMbed's own class scheme) are not reported: they are presumed already
    solvent-accessible.

    Downstream protocols can subtract these ROIs from any B-cell epitope
    candidate set (e.g. BepiPred/EpiDope/DiscoTope/ScanNet output) to
    discard regions that are not actually antibody-accessible on the
    mature, membrane-bound protein.
    """

    _label = 'tmbed predict'

    def _defineParams(self, form):
        form.addHidden(params.USE_GPU, params.BooleanParam, default=True,
                       label='Use GPU: ',
                       help='Whether to use GPU or not. (Unable to choose the GPU id).')
        form.addHidden(params.GPU_LIST, params.StringParam, default='0', label='Choose GPU IDs',
                       help='Add a list of GPU devices that can be used')
        form.addHidden('numThreads', params.IntParam, default=4,
                       help='Number of PyTorch CPU threads to use. Ignored when running on GPU.')

        form.addSection(label='Input')
        form.addParam('inputSequence', params.PointerParam, pointerClass='Sequence',
                       label='Input protein sequence: ',
                       help='Protein sequence to scan for transmembrane/signal regions.')

        gGroup = form.addGroup('Prediction')
        gGroup.addParam('minRegionLength', params.IntParam, label='Min. region length (aa): ', default=1,
                         help='Discard masking regions shorter than this after collapsing consecutive '
                              'same-class residues.')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.tmbedStep)
        self._insertFunctionStep(self.createOutputStep)

    # ---------------------------------- Steps -----------------------------------

    def _writeInputFasta(self):
        faFile = self._getExtraPath('inputSequence.fa')
        self.inputSequence.get().exportToFile(faFile)
        return os.path.abspath(faFile)

    def _getPredictionsPath(self):
        return self._getExtraPath('predictions.pred')

    def tmbedStep(self):
        faFile = self._writeInputFasta()
        predPath = self._getPredictionsPath()
        modelDir = immunoPlugin.getTMbedModelDir()

        args = (
            f'predict --fasta {faFile} --predictions {predPath} '
            f'--out-format {TMBED_OUT_FORMAT} --model-dir {modelDir} --threads {self.numThreads.get()} '
        )
        args += '--use-gpu' if self.useGpu.get() else '--no-use-gpu'
        immunoPlugin.runTMbed(self, args)

        if not os.path.isfile(predPath):
            raise TMbedExecutionError(
                f"TMbed finished without generating the expected output file at '{predPath}'."
            )

    def createOutputStep(self):
        predPath = self._getPredictionsPath()
        if not os.path.isfile(predPath):
            return

        records = parsePredictions(predPath)
        if not records:
            return

        # Single-sequence input: exactly one record is expected.
        _, classes = next(iter(records.values()))
        regions = extractMaskingRegions(classes, minLength=self.minRegionLength.get())

        inpSeq = self.inputSequence.get()
        fullSequence = inpSeq.getSequence()
        outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
        for region in regions:
            idxs = [region['start'], region['end']]
            roiSeq = Sequence(sequence=fullSequence[idxs[0] - 1:idxs[1]],
                               name='ROI_{}-{}'.format(*idxs), id='ROI_{}-{}'.format(*idxs),
                               description=f"TMbed {region['type']}")
            seqROI = SequenceROI(sequence=inpSeq, seqROI=roiSeq, roiIdx=idxs[0], roiIdx2=idxs[1])
            seqROI.setType(region['type'])
            outROIs.append(seqROI)

        if len(outROIs) > 0:
            self._defineOutputs(outputROIs=outROIs)
            self._defineSourceRelation(self.inputSequence, outROIs)

    # ---------------------------------- Validation -------------------------------

    def _validate(self):
        return immunoPlugin.validateTMbedInstallation()

    def _summary(self):
        summary = []
        if self.isFinished():
            outROIs = getattr(self, 'outputROIs', None)
            n = len(outROIs) if outROIs is not None else 0
            summary.append(f'{n} transmembrane/signal-peptide masking region(s) found.')
        return summary
