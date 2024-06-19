# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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
import os, webbrowser

from pyworkflow.protocol import params, Protocol
import pyworkflow.viewer as pwviewer

from pwchem.viewers.viewers_sequences import SequenceAliView

from ..protocols import ProtBepiPredPrediction

class ViewerBepiPred(pwviewer.ProtocolViewer):
  _label = 'BepiPred viewer'
  _targets = [ProtBepiPredPrediction]

  def _defineParams(self, form):
    form.addSection(label='Visualization of BepiPred epitoeps')
    group = form.addGroup('BepiPred scores view')
    group.addParam('displayBepiPred', params.LabelParam,
                   label='Display with BepiPred scores: ',
                   help='Display html with raw BepiPred scores')

    group = form.addGroup('BepiPred output epitopes view')
    group.addParam('displayROIs', params.LabelParam,
                   label='Display output BepiPred ROIs: ',
                   help='Display output ROIs as epitopes with AliView')

  def _getVisualizeDict(self):
    dispDic = {'displayBepiPred': self._showBepiPred,
               'displayROIs': self._showROIs}
    return dispDic

  def _showBepiPred(self, paramName=None):
      htmlFile = self.protocol._getExtraPath('output_interactive_figures.html')
      webbrowser.open(htmlFile)

  def _showROIs(self, paramName=None):
    obj = self.getProtocol().outputROIs
    outDir = self.getOutDir()
    outPath = os.path.join(outDir, f'viewSequences_{obj.getSequenceObj().getId()}.fasta')
    obj.exportToFile(outPath)

    seqFiles = [outPath]

    return [SequenceAliView(seqFiles, cwd=outDir)]

  def getOutDir(self):
    return os.path.abspath(self.getProtocol()._getExtraPath()
                           if self.getProtocol() else self.getProject().getTmpPath())

  def getProtocol(self):
    if hasattr(self, 'protocol') and isinstance(self.protocol, Protocol):
      return self.protocol