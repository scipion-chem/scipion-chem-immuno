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
This protocol is used to predict intrinsic antigenicity of a full-length
protein or construct with a local IApred installation.
"""

import math
from pathlib import Path

from pwchem.objects import SetOfSequenceROIs
from pwem.protocols import EMProtocol
from pyworkflow.object import Float, String
from pyworkflow.protocol import params

from immuno import Plugin as immunoPlugin
from ..utils.iapred_utils import parseOutput


class ProtIApredPrediction(EMProtocol):
    """
    AI Generated:

    Predicts intrinsic antigenicity of a set of sequences using a local
    IApred installation (SVM over aggregated physicochemical features,
    scoring each FULL sequence at once rather than per residue) and
    annotates every input ROI with the resulting score and verdict. It
    does NOT filter.

    It replaces VaxiJen, the historical reference for this task, which is
    not open-source and offers neither a downloadable standalone binary
    nor a documented public API.

    Output
    ------
    outputROIs: the same SetOfSequenceROIs as the input, with each ROI
    annotated with ``_iapredScore`` (NaN for sequences below IApred's own
    20 aa floor) and ``_iapredCategory``.
    """

    _label = 'iapred antigenicity'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputROIs', params.PointerParam, pointerClass='SetOfSequenceROIs',
                      label='Sequence ROIs: ',
                      help='Sequences to evaluate for intrinsic antigenicity (typically the single '
                           'assembled multi-epitope construct).')
        form.addParam('timeoutSeconds', params.IntParam, label='Timeout (s): ', default=300,
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Maximum time IApred is allowed to run before the step is aborted as '
                           'failed.')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.iapredStep)
        self._insertFunctionStep(self.createOutputStep)

    # ---------------------------------- Steps -----------------------------------

    def iapredStep(self):
        rois = self._getRois()
        sequences = [roi.getROISequence() for roi in rois]
        if not sequences:
            return

        # Resolved to absolute: the subprocess runs with cwd=IAPRED_HOME, so
        # a relative path here would resolve against THAT directory instead
        # of the project root.
        fastaPath = Path(self._getExtraPath('candidates.fasta')).resolve()
        with open(fastaPath, 'w') as fh:
            for i, sequence in enumerate(sequences):
                fh.write(f'>candidate_{i}\n{sequence}\n')

        rawOutputPath = Path(self._getRawOutputPath()).resolve()
        immunoPlugin.runIApred(self, f'{fastaPath} {rawOutputPath} -v',
                               cwd=immunoPlugin.getIApredDir())

    def createOutputStep(self):
        rois = self._getRois()
        if not rois:
            return

        resultDf = parseOutput(self._getRawOutputPath(), nExpected=len(rois))

        outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
        for roi, row in zip(rois, resultDf.itertuples(index=False)):
            # The 'sequence too short' case reaches here as NaN (parseOutput
            # coerces that literal text with pd.to_numeric(errors='coerce')) and
            # is stored as a null Float.
            isScored = not math.isnan(row.iapred_score)
            roi._iapredScore = Float(row.iapred_score) if isScored else Float(None)
            roi._iapredCategory = String(row.iapred_category)
            outROIs.append(roi)

        self._defineOutputs(outputROIs=outROIs)
        self._defineSourceRelation(self.inputROIs, outROIs)

    # ---------------------------------- Utils -----------------------------------

    def _getRawOutputPath(self):
        return self._getExtraPath('iapred_raw_output.csv')

    def _getRois(self):
        # Iterating a Scipion SetOfXXX reuses the same Python object per row
        # (the underlying sqlite cursor), so each item must be cloned when
        # materialized into a list or all N references end up pointing at
        # the cursor's last state.
        return [roi.clone() for roi in self.inputROIs.get()]

    # ---------------------------------- Info ------------------------------------

    def _validate(self):
        return immunoPlugin.validateIApredInstallation()

    def _summary(self):
        summary = []
        if self.isFinished():
            outROIs = getattr(self, 'outputROIs', None)
            if outROIs is not None:
                categories = [roi._iapredCategory.get() for roi in outROIs]
                summary.append(f'Categories: {", ".join(categories)}.')
        return summary
