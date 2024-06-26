# **************************************************************************
# *
# * Authors:	Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
This package contains protocols for creating and using IIITD Raghava software
"""

import multiprocessing

from scipion.install.funcs import InstallHelper

from pwchem import Plugin as pwchemPlugin

from .utils import *
from .constants import *

# Pluging variables
_logo = 'iiitd_logo.png'

class Plugin(pwchemPlugin):
	"""
	"""

	@classmethod
	def _defineVariables(cls):
		cls._defineVar(IIITD_DIC['activation'], cls.getEnvActivationCommand(IIITD_DIC))
		cls._defineVar(IIITD_DIC['browser'], 'Chrome')
		cls._defineVar(IIITD_DIC['browserPath'], '/usr/bin/google-chrome')

	@classmethod
	def defineBinaries(cls, env, default=True):
		"""This function defines the binaries for each package."""
		cls.addIIITDPackage(env)
		cls.addVaxignMLPackage(env)

	@classmethod
	def addIIITDPackage(cls, env, default=True):
		installer = InstallHelper(IIITD_DIC['name'], packageHome=cls.getVar(IIITD_DIC['home']),
															packageVersion=IIITD_DIC['version'])
		# Installing IIITD package
		installer.getCondaEnvCommand(pythonVersion='3.10', requirementsFile=False) \
			.addCondaPackages(['selenium'], channel='conda-forge') \
			.addPackage(env, ['git', 'conda'], default=default)

	@classmethod
	def addVaxignMLPackage(cls, env, default=True):
		# Installing Vaxign-ML package
		installer = InstallHelper(VAXIGNML_DIC['name'], packageHome=cls.getVar(VAXIGNML_DIC['home']),
															packageVersion=VAXIGNML_DIC['version'])
		installer.addCommand('docker pull e4ong1031/vaxign-ml:latest'). \
			addCommand('wget https://raw.githubusercontent.com/VIOLINet/Vaxign-ML-docker/master/VaxignML.sh'). \
			addCommand('chmod a+x VaxignML.sh').addCommand("sed -i 's/sudo//g' VaxignML.sh", 'VAXIGN_INSTALLED').\
			addPackage(env, dependencies=['chmod', 'wget', 'docker'], default=default)


	# ---------------------------------- Protocol functions-----------------------
	@classmethod
	def selectEpitopes(cls, selDics, jobs=1, browserData={}):
		'''Call the selectors specified in selecDic with the stored parameters using multiprocessing with n jobs.
			- selecDics : list of dictionaries as {selectorKey: {"software": softwareName, parameterName: parameterValue, }, }
			- jobs: number of jobs for multiprocessing

			Returns a Panda Dataframe with the selected epitopes with the following information columns:
			[Source, ProteinId, Position, Epitope, Score]
		'''
		# Create a pool of worker processes
		nJobs = len(selDics) if len(selDics) < jobs else jobs
		pool = multiprocessing.Pool(processes=nJobs)

		resultsDic = {}
		for selKey, selDic in selDics.items():
			softName = selDic['software']
			del selDic['software']
			resultsDic[(selKey, softName)] = pool.apply_async(runEpitopeSelection,
																												args=(softName, selDic, browserData))

		reportPoolStatus(resultsDic)

		pool.close()
		pool.join()

		epiDics = {}
		for (selKey, softName), res in resultsDic.items():
			epiDics[(selKey, softName)] = res.get()

		return epiDics

	@classmethod
	def performEvaluations(cls, sequences, evalDics, jobs=1, browserData={}, verbose=True):
		'''Generalize caller to the evaluation functions.
    - sequences: dict with sequences in the form: {seqId: sequence}
    - evalDics: dictionary as {evalKey: {parameterName: parameterValue}}
    - jobs: int, number of jobs for parallelization
    Returns a dictionary of the form: {(evalKey, softwareName): [scores]}
    '''
		funcDic = {
			'ToxinPred': callToxinPred, 'AlgPred2': callAlgPred2, 'ToxinPred2': callToxinPred2,
			'IL4pred': callIL4pred, 'IL10pred': callIL10pred, 'IFNepitope': callIFNepitope,
		}

		# Create a pool of worker processes
		nJobs = len(evalDics) if len(evalDics) < jobs else jobs
		pool = multiprocessing.Pool(processes=nJobs)

		resultsDic = {}
		for evalKey, evalDic in evalDics.items():
			softName = evalDic['software']
			smallEvalDic = evalDic.copy()
			del smallEvalDic['software']
			if softName in funcDic:
				resultsDic[(evalKey, softName)] = pool.apply_async(funcDic[softName],
																													 args=(sequences, smallEvalDic, browserData))

		if verbose:
			reportPoolStatus(resultsDic)

		pool.close()
		pool.join()

		epiDics = {}
		for (evalKey, softName), res in resultsDic.items():
			epiDics[(evalKey, softName)] = res.get()['Score']

		return epiDics

	# ---------------------------------- Utils functions-----------------------
	@classmethod
	def getBrowserData(cls):
		return {'name': cls.getVar(IIITD_DIC['browser']), 'path': cls.getVar(IIITD_DIC['browserPath'])}