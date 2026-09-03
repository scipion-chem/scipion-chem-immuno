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
This protocol is used to predict signal peptides of a full-length
protein/construct with a local SignalP-6.0 installation.
"""

from pathlib import Path

import pandas as pd
from pwchem.objects import SetOfSequenceROIs
from pwem.protocols import EMProtocol
from pyworkflow.object import Float, String
from pyworkflow.protocol import params

from immuno import Plugin as immunoPlugin
from ..constants import SIGNALP_DEFAULT_ORGANISM as DEFAULT_ORGANISM, SIGNALP_DIC
from ..utils.signalp_exceptions import SignalPExecutionError
from ..utils.signalp_utils import parse_output


class ProtSignalPPrediction(EMProtocol):
    """
    AI Generated:

    Predicts N-terminal signal peptides of a set of sequences using a
    local SignalP-6.0 installation, and annotates every input ROI with the
    result (does NOT filter).

    Purpose in the construct-level check (Fase 8): confirm that the final
    assembled multi-epitope construct does NOT have a predicted signal
    peptide at its N-terminal. A synthetic fusion construct meant for
    standard recombinant expression should not have one -- if SignalP
    predicts one, it signals that fragment junctions (or the first B-cell
    epitope itself) accidentally created a motif with that shape, worth
    reviewing before considering the construct valid. Purely informative
    here, not an automatic filter.

    Runs in 'slow-sequential' mode (~9.2GB weights): the full model run
    sequentially rather than in parallel, same RAM footprint as 'fast' but
    ~6x slower -- meant for CPU-only machines with limited RAM (the
    parallel 'slow' mode needs >14GB free).

    Output
    ------
    outputROIs: the same SetOfSequenceROIs as the input, with each ROI
    annotated with ``_signalpPrediction`` (``'OTHER'``/``'SP(Sec/SPI)'``/
    etc.), ``_signalpProbOther``, ``_signalpProbSp``,
    ``_signalpCsPosition`` (cleavage site position, empty string if
    prediction is ``'OTHER'``).
    """

    _label = 'signalp signal peptide'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputROIs', params.PointerParam, pointerClass='SetOfSequenceROIs',
                       label='Sequence ROIs: ',
                       help='Sequences to evaluate for an N-terminal signal peptide (typically the '
                            'single assembled multi-epitope construct).')
        form.addParam('organism', params.StringParam, default=DEFAULT_ORGANISM,
                       label='Organism group: ',
                       help="SignalP-6.0's '--organism' flag (e.g. 'other', 'eukarya').")
        form.addParam('timeoutSeconds', params.IntParam, label='Timeout (s): ', default=300,
                       expertLevel=params.LEVEL_ADVANCED)

    def _insertAllSteps(self):
        self._insertFunctionStep(self.signalpStep)
        self._insertFunctionStep(self.createOutputStep)

    # ---------------------------------- Steps -----------------------------------

    def _getRawResultsPath(self):
        return self._getExtraPath('signalp_out', 'prediction_results.txt')

    def _getRois(self):
        # Iterating a Scipion SetOfXXX reuses the same Python object per row
        # (the underlying sqlite cursor): each item must be cloned when
        # materialized into a list, or all N references end up pointing to
        # the cursor's last state.
        return [roi.clone() for roi in self.inputROIs.get()]

    def signalpStep(self):
        rois = self._getRois()
        sequences = [roi.getROISequence() for roi in rois]
        if not sequences:
            return

        # Resolved to absolute: the subprocess must receive absolute paths
        # regardless of its own cwd.
        fastaPath = Path(self._getExtraPath('candidates.fasta')).resolve()
        with open(fastaPath, 'w') as fh:
            for i, seq in enumerate(sequences):
                fh.write(f'>candidate_{i}\n{seq}\n')
        outputDir = Path(self._getExtraPath('signalp_out')).resolve()

        binary = immunoPlugin.getSignalPBinaryPath()
        args = (
            f'--fastafile {fastaPath} --output_dir {outputDir} --format none '
            f'--mode slow-sequential --organism {self.organism.get()} '
            f'--model_dir {immunoPlugin.getVar(SIGNALP_DIC["model_dir"])}'
        )
        self.runJob(binary, args)

        resultsPath = Path(self._getRawResultsPath())
        if not resultsPath.is_file():
            raise SignalPExecutionError(
                f"SignalP-6.0 finished without error but did not generate '{resultsPath}'."
            )

    def createOutputStep(self):
        rois = self._getRois()
        sequences = [roi.getROISequence() for roi in rois]
        if not sequences:
            return

        resultDf = parse_output(self._getRawResultsPath(), n_expected=len(sequences))

        outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
        for roi, row in zip(rois, resultDf.itertuples(index=False)):
            roi._signalpPrediction = String(row.signalp_prediction)
            roi._signalpProbOther = Float(row.signalp_prob_other)
            roi._signalpProbSp = Float(row.signalp_prob_sp)
            roi._signalpCsPosition = String(str(row.signalp_cs_position) if pd.notna(row.signalp_cs_position) else '')
            outROIs.append(roi)

        if len(outROIs) > 0:
            self._defineOutputs(outputROIs=outROIs)
            self._defineSourceRelation(self.inputROIs, outROIs)

    # ---------------------------------- Validation -------------------------------

    def _validate(self):
        return immunoPlugin.validateSignalPInstallation()

    def _summary(self):
        summary = []
        if self.isFinished():
            outROIs = getattr(self, 'outputROIs', None)
            if outROIs is not None:
                nWithSp = sum(1 for roi in outROIs if roi._signalpPrediction.get() != 'OTHER')
                summary.append(f'{nWithSp}/{len(outROIs)} candidate(s) with a predicted signal peptide.')
        return summary
