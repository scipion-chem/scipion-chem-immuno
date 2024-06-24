# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

import os

from pwem.protocols import EMProtocol
from pwem.objects import SetOfSequences
from pyworkflow.protocol import params

from pwchem.objects import SetOfSequenceROIs

from immuno import Plugin as immunoPlugin

class ProtVaxignMLEpitopeEvaluation(EMProtocol):
  """Run epitope evaluation on a set of epitopes (SetOfSequenceROIs)"""
  _label = 'Vaxign-ML epitope evaluation'

  _organOptions = ['Gram+', 'Gram-', 'Virus']
  _softParams = {'ABCpred': ['abcWindow', 'abcThres', 'abcFilter'],
                 'LBtope': ['lbModel', 'lbThres', 'lbLength']}

  def __init__(self, **kwargs):
    EMProtocol.__init__(self, **kwargs)

  def _defineParams(self, form):
    form.addSection(label='Input')
    iGroup = form.addGroup('Input')
    iGroup.addParam('inputSource', params.EnumParam, label='Input: ',
                    default=0, choices=['Set of Sequences', 'Set of sequence ROIs'],
                    display=params.EnumParam.DISPLAY_HLIST,
                    help="Whether to predict the protegenicity on a set of sequences or sequence ROIs")
    iGroup.addParam('inputSequences', params.PointerParam, pointerClass="SetOfSequences", allowsNull=True,
                    label='Input protein sequences: ', condition='inputSource==0',
                    help="Protein sequences to perform the screening on")
    iGroup.addParam('inputROIs', params.PointerParam, pointerClass="SetOfSequenceROIs",
                    label='Input sequence ROIs: ', condition='inputSource==1', allowsNull=True,
                    help="Set of sequence ROIs to label with the present MHC-II alleles")

    aGroup = form.addGroup('Parameters')
    aGroup.addParam('organ', params.EnumParam, choices=self._organOptions,
                    label='Choose epitope source: ', default=0,
                    help='Organism type for the pathogenic epitopes.')

    form.addParallelSection(threads=4, mpi=1)

  def _insertAllSteps(self):
    self._insertFunctionStep(self.evaluationStep)
    self._insertFunctionStep(self.defineOutputStep)

  def evaluationStep(self):
    if self.countViableSequences() > 0:
      inFasta = self.generateFasta()
      org = self.getEnumText('organ').lower()
      kwargs = {'i': inFasta, 'o': self.getOutDir(), 't': org, 'p': self.numberOfThreads.get()}

      immunoPlugin.runVaxignML(self, kwargs, cwd=self._getExtraPath())
    else:
      print('There are no sequences over 50 residues in the input. VaxignML cannot be executed')

  def defineOutputStep(self):
      outFile = self.getOutDir('inputSequences.result.tsv')
      if os.path.exists(outFile):
        resDic = self.parseResults(outFile)

        if self.inputSource.get() == 0:
          outSeqs = SetOfSequences().create(outputPath=self._getPath())
          for seq in self.inputSequences.get():
            seqId = seq.getId()
            seq._vaxignML = params.Float(float(resDic[seqId][1]))
            outSeqs.append(seq)

          self._defineOutputs(outputSequences=outSeqs)

        else:
          inROIs = self.inputROIs.get()
          outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
          for roi in inROIs:
            roiId = roi.getROIId()
            roi._vaxignML = params.Float(float(resDic[roiId][1]))
            outROIs.append(roi)

          self._defineOutputs(outputROIs=outROIs)


  ##################### UTILS #####################
  def getInputSet(self):
    return self.inputSequences.get() if self.inputSource.get() == 0 else self.inputROIs.get()

  def getInputSeqFunc(self):
    return 'getSequence' if self.inputSource.get() == 0 else 'getROISequence'

  def getInputNameFunc(self):
    return 'getSeqName' if self.inputSource.get() == 0 else 'getROIId'
  
  def countViableSequences(self):
    count = 0
    inSet = self.getInputSet()
    seqFunc, nameFunc = self.getInputSeqFunc(), self.getInputNameFunc()
    for item in inSet:
      seq, name = getattr(item, seqFunc)(), getattr(item, nameFunc)()
      if len(seq) >= 50:
        count += 1
    return count
  
  def getOutDir(self, path=''):
    return os.path.join(os.path.abspath(self._getExtraPath('vaxResults')), path)

  def generateFasta(self):
    inFasta = os.path.abspath(self._getExtraPath('inputSequences.fa'))
    if self.inputSource.get() == 0:
      self.inputSequences.get().exportToFile(inFasta)
    else:
      self.inputROIs.get().exportToFile(inFasta, mainSeq=False)
    return inFasta
  
  def parseResults(self, resFile):
    resDic = {}
    with open(resFile) as f:
      f.readline()
      for line in f:
        sline = line.strip().split()
        resDic[sline[0]] = sline[-2:]
    return resDic
  
  
  def _warnings(self):
    warns = []
    inSet = self.getInputSet()
    funcStr = self.getInputSeqFunc()

    for item in inSet:
      if len(getattr(item, funcStr)()) < 50:
        warns.append('The input contains sequences less than 50 residues long. VaxignML cannot work with them, so those '
                    'sequence ROIs will have VaxignML score of 0')
        break
    return warns
