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
This protocol is used to annotate MHC-I C-terminal antigen processing
(proteasomal cleavage) evidence on a set of peptide candidates with a
local NetCleave installation.
"""

import glob
import os
import tempfile
from pathlib import Path

import pandas as pd
from pwchem.objects import SetOfSequenceROIs
from pwem.protocols import EMProtocol
from pyworkflow.object import Boolean, Float
from pyworkflow.protocol import params

from immuno import Plugin as immunoPlugin
from ..constants import (
    NETCLEAVE_MHC_CLASS as MHC_CLASS,
    NETCLEAVE_MHC_FAMILY as MHC_FAMILY,
    NETCLEAVE_TECHNIQUE as TECHNIQUE,
)
from ..utils.netcleave_exceptions import NetCleaveExecutionError
from ..utils.netcleave_utils import findCTermMatch, parseOutput


class ProtNetCleavePrediction(EMProtocol):
    """
    AI Generated:

    Annotates (does NOT filter) every input ROI with evidence of real
    proteasomal C-terminal cleavage, using a local NetCleave installation.
    Complementary signal, not a filter: a peptide can bind MHC-I strongly
    but never actually be generated via antigen processing if the
    proteasome does not cut exactly where the accepted binding core ends.

    Meant to be run on the OUTPUT of an MHC-I promiscuity protocol (e.g.
    ``scipion-chem-netmhcpan``), whose ROIs already carry a ``_core9aa``
    attribute (the accepted binding core, narrower than the full evaluated
    window) -- if that attribute is absent, the ROI's own full sequence is
    used as its own "core" (cleavage checked right after the ROI's last
    residue instead).

    Output
    ------
    outputROIs: the same SetOfSequenceROIs as the input, annotated with
    ``_netcleaveCTermMatch`` (bool) and ``_netcleaveCTermScore`` (float,
    only meaningful if the match is True).
    """

    _label = 'netcleave c-term cleavage'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputROIs', params.PointerParam, pointerClass='SetOfSequenceROIs',
                       label='Sequence ROIs: ',
                       help='Peptide candidates (typically MHC-I promiscuity protocol output) to '
                            'annotate with C-terminal cleavage evidence.')
        form.addParam('timeoutSeconds', params.IntParam, label='Timeout (s) per candidate: ', default=300,
                       expertLevel=params.LEVEL_ADVANCED,
                       help='Maximum time NetCleave is allowed to run for a single candidate before '
                            'the step is aborted as failed. Increase on slow/CPU-only hardware.')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.netcleaveStep)
        self._insertFunctionStep(self.createOutputStep)

    # ---------------------------------- Steps -----------------------------------

    def _getRois(self):
        # Iterating a Scipion SetOfXXX reuses the same Python object per row
        # (the underlying sqlite cursor): each item must be cloned when
        # materialized into a list, or all N references end up pointing to
        # the cursor's last state.
        return [roi.clone() for roi in self.inputROIs.get()]

    def _getCore(self, roi):
        # A Scipion SetOfXXX fixes its column schema from the first
        # appended item's dynamic attributes: an upstream protocol that
        # sets '_core9aa' on some ROIs must set it (even if empty) on all
        # of them, so an empty/blank value means "no narrower core given"
        # just as much as the attribute being absent entirely.
        core = getattr(roi, '_core9aa', None)
        coreSeq = core.get() if core is not None else None
        return coreSeq if coreSeq else roi.getROISequence()

    def netcleaveStep(self):
        rois = self._getRois()
        scriptPath = Path(immunoPlugin.getNetCleaveScriptPath())

        # The subprocess below runs with cwd=<NetCleave's own install dir>
        # (required so it finds its bundled model, see the plugin's
        # constants docstring): a RELATIVE fastaPath (self._getExtraPath()
        # is relative to the project root, not an absolute path) would
        # resolve against the CHILD process's cwd instead, not this
        # protocol's. Always resolve to absolute before building the
        # subprocess command.
        with tempfile.TemporaryDirectory(prefix='netcleave_', dir=Path(self._getExtraPath()).resolve()) as tmp:
            tmpDir = Path(tmp).resolve()
            for i, roi in enumerate(rois):
                windowSeq = roi.getROISequence()
                fastaPath = tmpDir / f'seq_{i}.fasta'
                fastaPath.write_text(f'>candidate_{i}\n{windowSeq}\n', encoding='utf-8')

                args = (
                    f'--mhc_class {MHC_CLASS} --technique {TECHNIQUE} --mhc_family {MHC_FAMILY} '
                    f'--score_fasta {fastaPath}'
                )
                immunoPlugin.runNetCleave(self, args, cwd=str(scriptPath.parent))

                # NetCleave names its output '<fasta_stem>_<first_token_of_
                # fasta_header>_NetCleave.xlsx' (undocumented in --help,
                # verified by reading source): glob instead of reconstructing
                # the exact name, to be robust to naming changes.
                matches = glob.glob(str(tmpDir / f'seq_{i}_*_NetCleave.xlsx'))
                if not matches:
                    raise NetCleaveExecutionError(
                        f"NetCleave finished without error but did not generate the expected "
                        f"'seq_{i}_*_NetCleave.xlsx' file in '{tmpDir}' (candidate {i})."
                    )
                cleavageDf = parseOutput(Path(matches[0]))
                cleavageDf.to_csv(self._getExtraPath(f'candidate_{i}_cleavage.csv'), index=False)

    def createOutputStep(self):
        rois = self._getRois()

        outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
        for i, roi in enumerate(rois):
            cleavageCsv = self._getExtraPath(f'candidate_{i}_cleavage.csv')
            cleavageDf = pd.read_csv(cleavageCsv) if os.path.isfile(cleavageCsv) else None

            score = None
            if cleavageDf is not None:
                score = findCTermMatch(cleavageDf, self._getCore(roi), roi.getROISequence())

            roi._netcleaveCTermMatch = Boolean(score is not None)
            roi._netcleaveCTermScore = Float(score) if score is not None else Float(None)
            outROIs.append(roi)

        if len(outROIs) > 0:
            self._defineOutputs(outputROIs=outROIs)
            self._defineSourceRelation(self.inputROIs, outROIs)

    # ---------------------------------- Validation -------------------------------

    def _validate(self):
        return immunoPlugin.validateNetCleaveInstallation()

    def _summary(self):
        summary = []
        if self.isFinished():
            outROIs = getattr(self, 'outputROIs', None)
            if outROIs is not None:
                nMatch = sum(1 for roi in outROIs if roi._netcleaveCTermMatch.get())
                summary.append(f'{nMatch}/{len(outROIs)} candidate(s) with a confirmed C-terminal cleavage site.')
        return summary
