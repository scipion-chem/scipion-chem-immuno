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

from pwem.protocols import EMProtocol
from pyworkflow.protocol import params

from pwchem.objects import SetOfSequenceROIs

from .. import Plugin as iiitdPlugin
from ..constants import TOXIN2WARN
from ..utils import mapEvalParamNames

class ProtIIITDEvaluations(EMProtocol):
  """Run evaluations on a set of epitopes (SetOfSequenceROIs)"""
  _label = 'immuno epitope evaluations'

  _evaluatorOptions = ['ToxinPred', 'AlgPred2', 'IL4pred', 'IL10pred', 'IFNepitope', 'ToxinPred2']

  _toxinSVMMethods = ["SVM(Swiss-Prot)", "SVM(Swiss-Prot)+Motif", "SVM(TrEMBL)", "SVM(TrEMBL)+Motif"]
  _toxinQMMethods = ["Monopeptide(Swiss-Prot)", "Monopeptide(TrEMBL)", "Dipeptide(Swiss-Prot)", "Dipeptide(TrEMBL)"]
  _algMethods = ["AAC based RF", "Hybrid (RF+BLAST+MERCI)"]
  _il4Methods = ["SVM", "MERCI", "Hybrid", "Swiss-prot"]
  _il10Methods = ["SVM", "Random Forest"]
  _ifnMethods = ["Motif", "SVM", "Hybrid"]
  _ifnModels = ["IFN-gamma versus Non IFN-gamma", "IFN-gamma versus other cytokine", "IFN-gamma versus random"]
  _toxin2Methods = ["AAC based RF", "Hybrid (RF+BLAST+MERCI)"]


  _softParams = {'ToxinPred': ['toxinMethod', 'toxinSVMMethod', 'toxinQMMethod', 'toxinEval', 'toxinThval'],
                 'AlgPred2': ['algMethod', 'algThres'],
                 'IL4pred': ['il4Method', 'il4Thval'],
                 'IL10pred': ['il10Method', 'il10Thval'],
                 'IFNepitope': ['ifnMethod'],
                 'ToxinPred2': ['toxin2Method', 'toxin2Eval']
                 }

  def __init__(self, **kwargs):
    EMProtocol.__init__(self, **kwargs)

  def _defineEvalParams(self, aGroup, allCond=True):
    '''Define the evaluation options and the parameters for each of them.
    allCond: condition to apply for all the parameters

    WARNING: This function is used by a scipion-chem metaprotocol to use and define this parameters by its own,
    modify with care'''
    aGroup.addParam('chooseIIITDEvaluator', params.EnumParam, choices=self._evaluatorOptions,
                    label='Choose evaluator: ', default=0, condition=f'{allCond}',
                    help=f'Epitope evaluation software to use.\n{TOXIN2WARN}')

    aGroup.addParam('toxinMethod', params.EnumParam, choices=['SVM', 'Quantitative Matrix'], default=0,
                    label='ToxinPred method: ', display=params.EnumParam.DISPLAY_HLIST, 
                    condition=f'{allCond} and chooseIIITDEvaluator==0',
                    help='Which kind of ToxinPred method to use')
    aGroup.addParam('toxinSVMMethod', params.EnumParam, choices=self._toxinSVMMethods, default=0, 
                    label='ToxinPred SVM method: ', 
                    condition=f'{allCond} and chooseIIITDEvaluator==0 and toxinMethod==0',
                    help='Which SVM ToxinPred method to use')
    aGroup.addParam('toxinQMMethod', params.EnumParam, choices=self._toxinQMMethods, default=0,
                    label='ToxinPred QM method: ', condition=f'{allCond} and chooseIIITDEvaluator==0 and toxinMethod==1',
                    help='Which Quantitative Matrix ToxinPred method to use')
    aGroup.addParam('toxinEval', params.FloatParam, label='E-value cutoff: ', default=10.0,
                    condition=f'{allCond} and chooseIIITDEvaluator==0 and toxinMethod==0 and toxinSVMMethod in [1, 3]',
                    help='E-value to MAST search in Motif based methods')
    aGroup.addParam('toxinThval', params.FloatParam, label='SVM threshold: ', default=0.0,
                    condition=f'{allCond} and chooseIIITDEvaluator==0 and toxinMethod==0 and toxinSVMMethod in [0, 2]',
                    help='Threshold for the SVM predictions (-1, 1)')

    aGroup.addParam('algMethod', params.EnumParam, label='AlgPred2 model: ', default=0,
                    condition=f'{allCond} and chooseIIITDEvaluator==1', choices=self._algMethods,
                    help='Machine Learning Technique used for developing model.')
    aGroup.addParam('algThres', params.FloatParam, choices=["on", "off"],
                    label='AlgPred2 threshold: ', default=0.3, condition=f'{allCond} and chooseIIITDEvaluator==1',
                    help='Threshold for the predictions (-0.5, 2)')

    aGroup.addParam('il4Method', params.EnumParam, choices=self._il4Methods, label='IL4pred method to use: ', default=2,
                    condition=f'{allCond} and chooseIIITDEvaluator==2', help='IL4pred model to use for epitope evaluation')
    aGroup.addParam('il4Thval', params.FloatParam, label='SVM threshold: ', default=0.2,
                    condition=f'{allCond} and chooseIIITDEvaluator==2 and il4Method==0',
                    help='Threshold for the SVM predictions (-1, 1)')

    aGroup.addParam('il10Method', params.EnumParam, choices=self._il10Methods, label='IL10pred method to use: ',
                    default=0, condition=f'{allCond} and chooseIIITDEvaluator==3',
                    help='IL10pred model to use for epitope evaluation')
    aGroup.addParam('il10Thval', params.FloatParam, label='IL10pred threshold: ', default=-0.3,
                    condition=f'{allCond} and chooseIIITDEvaluator==3',
                    help='Threshold for the SVM/Random Forest predictions (-1, 1)')

    aGroup.addParam('ifnMethod', params.EnumParam, choices=self._ifnMethods, label='IFNepitope approach to use: ',
                    default=2, condition=f'{allCond} and chooseIIITDEvaluator==4',
                    help='IFNepitope approach to use for epitope evaluation')

    aGroup.addParam('toxin2Method', params.EnumParam, choices=self._toxin2Methods, default=0,
                    label='ToxinPred2 method: ', condition=f'{allCond} and chooseIIITDEvaluator==5',
                    help=f'Which ToxinPred2 method to use.\n{TOXIN2WARN}')
    aGroup.addParam('toxin2Eval', params.FloatParam, label='Threshold value: ', default=0.6,
                    condition=f'{allCond} and chooseIIITDEvaluator==5',
                    help=f'Threshold for the SVM predictions (-0.5, 2).\n{TOXIN2WARN}')

    return aGroup

  def _defineParams(self, form):
    form.addSection(label='Input')
    iGroup = form.addGroup('Input')
    iGroup.addParam('inputROIs', params.PointerParam, pointerClass="SetOfSequenceROIs",
                    label='Input epitopes: ',
                    help="Input set of epitope sequences as SetOfSequenceROIs")

    form.addSection(label='Add evaluations')
    aGroup = form.addGroup('Define evaluator')
    aGroup = self._defineEvalParams(aGroup)
    aGroup.addParam('evaluatorIIITDName', params.StringParam, label='Evaluator name: ',
                    default='', expertLevel=params.LEVEL_ADVANCED,
                    help='Set the name for the defined evaluator.')
    aGroup.addParam('addEval', params.LabelParam, label='Add defined evaluator: ',
                    help='Add defined evaluator to perform the epitope prediction')
    sGroup = form.addGroup('Evaluators summary')
    sGroup.addParam('inEvals', params.TextParam, width=70, default='',
                    label='Evaluators summary: ',
                    help='Summary of the epitope evaluations that will be performed')

    form.addParallelSection(threads=4, mpi=1)


  def _insertAllSteps(self):
    self._insertFunctionStep(self.evaluationStep)

  def evaluationStep(self):
    nt = self.numberOfThreads.get()
    sDics = self.getWebEvaluatorDics()
    sequences = self.getInputSequences()

    epiDic = iiitdPlugin.performEvaluations(sequences, sDics, nt, iiitdPlugin.getBrowserData())

    outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
    for i, roi in enumerate(self.inputROIs.get()):
      for (evalKey, softName), scores in epiDic.items():
        setattr(roi, evalKey, params.Float(scores[i]))
      outROIs.append(roi)

    if len(outROIs) > 0:
      self._defineOutputs(outputROIs=outROIs)


  ##################### UTILS #####################
  def getInputSequences(self):
    seqs = {}
    for roi in self.inputROIs.get():
      seqs[roi.getROIId()] = roi.getROISequence()
    return seqs

  def buildElementDic(self):
    sName, soft = self.evaluatorIIITDName.get(), self.getEnumText('chooseIIITDEvaluator')
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
    for line in self.inEvals.get().split('\n'):
      if line.strip():
        sd = f'{{{line.split(") ")[1]}}}'
        sDic.update(eval(sd))
    return sDic

  def getDefSName(self, soft):
    sDic, i = self.parseElementsDic(), 1
    sName = f'{soft}-{i}'
    while sName in sDic:
      i += 1
      sName = f'{soft}-{i}'
    return sName

  def getParamValue(self, paramName):
    if isinstance(self.getParam(paramName), params.EnumParam):
      value = self.getEnumText(paramName)
    else:
      value = getattr(self, paramName).get()
    return value

  def getWebEvaluatorDics(self):
    ''' Returns the selector dictionary with the parameter names expected by the web server
    :return: dic, {selName: {software: softName, paramName: paramValue}} with the webserver chosen parameters
    '''
    sDic = self.parseElementsDic()
    wsDic = mapEvalParamNames(sDic)
    return wsDic

  def _validate(self):
    vs = []
    if len(self.getWebEvaluatorDics()) < 1:
      vs.append('You need to add at least one evaluator to run the protocol')
    return vs

  def _summary(self):
    sm = []
    if self.inEvals.get().strip():
      sm.append(self.inEvals.get())
    return sm
