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
from pyworkflow.protocol import params

from pwchem.objects import Sequence, SequenceROI, SetOfSequenceROIs

from immuno import Plugin as iiitdPlugin
from ..constants import SEL_PARAM_MAP

class ProtIIITDEpitopeSelection(EMProtocol):
  """Run epitope selections on a set of protein sequences (SetOfSequences)"""
  _label = 'immuno epitope selection'

  _selectorOptions = ['ABCpred', 'LBtope']
  _lbModels = ['LBtope_Fixed', 'LBtope_Fixed_non_redundant',
               'LBtope_Variable', 'LBtope_Variable_non_redundant', 'LBtop_Confirm']

  _softParams = {'ABCpred': ['abcWindow', 'abcThres', 'abcFilter'],
                 'LBtope': ['lbModel', 'lbThres', 'lbLength']}

  def __init__(self, **kwargs):
    EMProtocol.__init__(self, **kwargs)

  def _defineParams(self, form):
    form.addSection(label='Input')
    iGroup = form.addGroup('Input')
    iGroup.addParam('inputSequence', params.PointerParam, pointerClass="Sequence",
                    label='Input sequence: ',
                    help="Input protein sequence as where epitopes will be identified")

    form.addSection(label='Add selectors')
    aGroup = form.addGroup('Define selector')
    aGroup.addParam('chooseSelector', params.EnumParam, choices=self._selectorOptions,
                    label='Choose selector: ', default=0,
                    help='Epitope selection software to use.')
    aGroup.addParam('selectorName', params.StringParam, label='Selector name: ',
                    default='', expertLevel=params.LEVEL_ADVANCED,
                    help='Set the name for the defined selector.')

    aGroup.addParam('abcWindow', params.EnumParam, choices=["16", "18"], label='ABCpred window size: ', default=0,
                    display=params.EnumParam.DISPLAY_HLIST, condition='chooseSelector==0',
                    help='Window size for the prediction. Correspond to the size of the predicted epitopes.')
    aGroup.addParam('abcThres', params.FloatParam, label='ABCpred threshold: ', default=0.51,
                    condition='chooseSelector==0', expertLevel=params.LEVEL_ADVANCED,
                    help='Threshold to consider a positive predicted epitope (from 0.1 to 1).')
    aGroup.addParam('abcFilter', params.EnumParam, choices=["on", "off"], 
                    label='ABCpred overlap filter: ', default=0, condition='chooseSelector==0',
                    expertLevel=params.LEVEL_ADVANCED, display=params.EnumParam.DISPLAY_HLIST,
                    help='Window size for the prediction. Correspond to the size of the predicted epitopes.')

    aGroup.addParam('lbModel', params.EnumParam, choices=self._lbModels, label='LBtope model to use: ', default=2,
                    condition='chooseSelector==1', help='LBtope model to use for epitope prediction')
    aGroup.addParam('lbThres', params.EnumParam, choices=["20", "40", "60", "80"], condition='chooseSelector==1', 
                    label='LBtope probability threshold: ', default=2, expertLevel=params.LEVEL_ADVANCED,
                    help='Probability threshold to consider a positive predicted epitope.')
    aGroup.addParam('lbLength', params.IntParam, label='LBtope epitope lenght: ', default=15,
                    condition='chooseSelector==1 and lbModel>1',
                    help='Size of the predicted epitopes.')

    aGroup.addParam('addSel', params.LabelParam, label='Add defined selector: ',
                   help='Add defined selector to perform the epitope prediction')

    sGroup = form.addGroup('Selectors summary')
    sGroup.addParam('inSels', params.TextParam, width=70, default='',
                   label='Selectors summary: ',
                   help='Summary of the epitope selections that will be performed')

    form.addParallelSection(threads=4, mpi=1)


  def _insertAllSteps(self):
    self._insertFunctionStep(self.selectionStep)

  def selectionStep(self):
    nt = self.numberOfThreads.get()
    sDics = self.getWebSelectorDics()
    sDics = self.addInputSequences(sDics)

    epiDics = iiitdPlugin.selectEpitopes(sDics, nt, iiitdPlugin.getBrowserData())

    inpSeq = self.inputSequence.get()
    outROIs = SetOfSequenceROIs(filename=self._getPath(f'sequenceROIs.sqlite'))

    for (selKey, softName), epiDic in epiDics.items():
      seqEpDic = epiDic[f'seq1']
      if seqEpDic:
        for epSeq, epIdx, epSc in zip(seqEpDic['Sequence'], seqEpDic['Position'], seqEpDic['Score']):
          idxs = [int(epIdx), int(epIdx) + len(epSeq)]
          roiName = '{}_ROI_{}-{}'.format(selKey, *idxs)
          roiSeq = Sequence(sequence=epSeq, name=roiName, id=roiName,
                            description=f'{selKey} epitope')
          seqROI = SequenceROI(sequence=inpSeq, seqROI=roiSeq, roiIdx=idxs[0], roiIdx2=idxs[1])
          seqROI._epitopeType = params.String('B')
          seqROI._source = params.String(softName)
          setattr(seqROI, softName, params.Float(epSc))
          outROIs.append(seqROI)

    if len(outROIs) > 0:
      self._defineOutputs(**{f'outputROIs': outROIs})

  ##################### UTILS #####################
  def addInputSequences(self, sDics):
    faFile = self._getExtraPath('inputSequence.fa')
    self.inputSequence.get().exportToFile(faFile)
    for sName in sDics:
      sDics[sName].update({'i': faFile})
    return sDics

  def buildElementDic(self):
    sName, soft = self.selectorName.get(), self.getEnumText('chooseSelector')
    if not sName.strip():
      sName = self.getDefSName(soft)

    sDic = {sName: {'software': soft}}
    for paramName in self._softParams[soft]:
      sDic[sName].update({paramName: self.getParamValue(paramName)})
    return sDic

  def parseElementsDic(self):
    ''' Parse the selector dictionaries included in the input list
    :return: dic, {selName: {software: softName, paramName: paramValue}} with the chosen Scipion parameters
    '''
    sDic = {}
    for line in self.inSels.get().split('\n'):
      if line.strip():
        sd = f'{{{line.split(")")[1]}}}'
        sDic.update(eval(sd))
    return sDic

  def getWebSelectorDics(self):
    ''' Returns the selector dictionary with the parameter names expected by the web server
    :return: dic, {selName: {software: softName, paramName: paramValue}} with the webserver chosen parameters
    '''
    wsDic = {}
    sDic = self.parseElementsDic()
    for sName, curSDic in sDic.items():
      wsDic[sName] = {}
      for paramName, paramValue in curSDic.items():
        if paramName in SEL_PARAM_MAP:
          paramName = SEL_PARAM_MAP[paramName]
        wsDic[sName][paramName] = paramValue
    return wsDic


  def getParamValue(self, paramName):
    if isinstance(self.getParam(paramName), params.EnumParam):
      value = self.getEnumText(paramName)
    else:
      value = getattr(self, paramName).get()
    return value

  def getDefSName(self, soft):
    sDic, i = self.parseElementsDic(), 1
    sName = f'{soft}-{i}'
    while sName in sDic:
      i += 1
      sName = f'{soft}-{i}'
    return sName

  def _validate(self):
    vs = []
    if len(self.getWebSelectorDics()) < 1:
      vs.append('You need to add at least one selector to run the protocol')
    return vs

  def _summary(self):
    sm = []
    if self.inSels.get().strip():
      sm.append(self.inSels.get())
    return sm