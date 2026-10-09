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

import multiprocessing, shutil, subprocess

from scipion.install.funcs import InstallHelper

from pwchem.utils import insistentRun, getReplaceCommand

from .utils import *
from .constants import *

# Pluging variables
_logo = 'immuno_logo.png'


class Plugin(pwchemPlugin):
	"""
	"""

	@classmethod
	def _defineVariables(cls):
		cls._defineVar(IIITD_DIC['activation'], cls.getEnvActivationCommand(IIITD_DIC))
		cls._defineVar(IIITDW_DIC['activation'], cls.getEnvActivationCommand(IIITDW_DIC))
		cls._defineVar(IIITDW_DIC['browser'], 'Chrome')
		cls._defineVar(IIITDW_DIC['browserPath'], '/usr/bin/google-chrome')

		cls._defineEmVar(VAXIGNML_DIC['home'], VAXIGNML_DIC['name'] + '-' + VAXIGNML_DIC['version'])

		cls._defineEmVar(EPIDOPE_DIC['home'], cls.getEnvName(EPIDOPE_DIC))

		cls._defineEmVar(SCANNET_DIC['home'], cls.getEnvName(SCANNET_DIC))
		cls._defineVar(SCANNET_DIC['activation'], cls.getEnvActivationCommand(SCANNET_DIC))

		cls._defineEmVar(DISCOTOPE_DIC['home'], cls.getEnvName(DISCOTOPE_DIC))
		cls._defineVar(DISCOTOPE_DIC['activation'], cls.getEnvActivationCommand(DISCOTOPE_DIC))

		cls._defineEmVar(TMBED_DIC['home'], cls.getEnvName(TMBED_DIC))
		cls._defineVar(TMBED_DIC['activation'], cls.getEnvActivationCommand(TMBED_DIC))

		# SignalP-6.0 never auto-installs (academic license): no _defineEmVar
		cls._defineVar(SIGNALP_DIC['python_bin'], '')
		cls._defineVar(SIGNALP_DIC['binary_name'], SIGNALP_DEFAULT_BINARY_NAME)
		cls._defineVar(SIGNALP_DIC['model_dir'], '')

		cls._defineEmVar(IAPRED_DIC['home'], cls.getEnvName(IAPRED_DIC))
		cls._defineVar(IAPRED_DIC['activation'], cls.getEnvActivationCommand(IAPRED_DIC))
		cls._defineVar(IAPRED_DIC['script_name'], IAPRED_DEFAULT_SCRIPT_NAME)

		cls._defineEmVar(NETCLEAVE_DIC['home'], cls.getEnvName(NETCLEAVE_DIC))
		cls._defineVar(NETCLEAVE_DIC['activation'], cls.getEnvActivationCommand(NETCLEAVE_DIC))

		cls._defineEmVar(STACKGLYEMBED_DIC['home'], cls.getEnvName(STACKGLYEMBED_DIC))
		cls._defineVar(STACKGLYEMBED_DIC['activation'], cls.getEnvActivationCommand(STACKGLYEMBED_DIC))
		cls._defineVar(STACKGLYEMBED_DIC['esm_model_name'], STACKGLYEMBED_DEFAULT_ESM_MODEL_NAME)

	@classmethod
	def defineBinaries(cls, env, default=True):
		"""This function defines the binaries for each package."""
		cls.addIIITDPackage(env)
		cls.addIIITDWPackage(env)
		cls.addIL6PredPackage(env)
		cls.addVaxignMLPackage(env)
		cls.addEpiDopePackage(env)
		cls.addScanNetPackage(env)
		cls.addDiscoTopePackage(env)
		cls.addTMbedPackage(env)
		cls.addIApredPackage(env)
		cls.addNetCleavePackage(env)
		cls.addStackGlyEmbedPackage(env)
		# SignalP-6.0: no-op, never installed automatically (academic license,
		# not redistributable) -- see validateInstallation/README.rst.

	@classmethod
	def addIIITDPackage(cls, env, default=True):
		installer = InstallHelper(IIITD_DIC['name'], packageHome=cls.getVar(IIITD_DIC['home']),
															packageVersion=IIITD_DIC['version'])
		# Installing IIITD package
		installer.getCondaEnvCommand(pythonVersion=IIITD_DIC['python'], requirementsFile=False) \
			.addCommand(f'{cls.getEnvActivationCommand(IIITD_DIC)} && pip install {" ".join(IIITD_PACKAGES)}',
									'PIP_MODS_INSTALLED') \
			.addCommand(f'{cls.getFixScripts(IIITD_DIC, IIITD_FIXES)}', 'SCRIPS_FIXED') \
			.addPackage(env, ['conda', 'pip'], default=default)

	@classmethod
	def addIIITDWPackage(cls, env, default=True):
		installer = InstallHelper(IIITDW_DIC['name'], packageHome=cls.getVar(IIITDW_DIC['home']),
															packageVersion=IIITDW_DIC['version'])
		# Installing IIITD package
		installer.getCondaEnvCommand(pythonVersion='3.10', requirementsFile=False) \
			.addCondaPackages(['selenium'], channel='conda-forge') \
			.addPackage(env, ['conda'], default=default)


	@classmethod
	def addIL6PredPackage(cls, env, default=True):
		installer = InstallHelper(IL6PRED_DIC['name'], packageHome=cls.getVar(IL6PRED_DIC['home']),
															packageVersion=IL6PRED_DIC['version'])
		# Installing IL6PRED package
		installer.getCondaEnvCommand(pythonVersion=IL6PRED_DIC['python'], requirementsFile=False) \
			.addCondaPackages(['tqdm'], channel='conda-forge') \
			.addCommand(f'{cls.getEnvActivationCommand(IL6PRED_DIC)} && pip install il6pred', 'IL6PRED_PIP_INSTALLED') \
			.addCommand(f'{cls.getFixScripts(IL6PRED_DIC, IL6_FIXES)}', 'SCRIPS_FIXED') \
			.addPackage(env, ['conda', 'pip'], default=default)

	@classmethod
	def addVaxignMLPackage(cls, env, default=True):
		# Installing Vaxign-ML package
		# VAXIGN_INSTALLED = '%s_installed' % VAXIGNML_DIC['name']
		# vaxignML_commands = f'docker images -q e4ong1031/vaxign-ml > {VAXIGN_INSTALLED} && '
		# vaxignML_commands += f'find ./ -maxdepth 1 -size 0 -name {VAXIGN_INSTALLED} -delete'
		#
		# vaxignML_commands = [(vaxignML_commands, VAXIGN_INSTALLED)]
		# env.addPackage(VAXIGNML_DIC['name'], version=VAXIGNML_DIC['version'],
		# 							 tar='void.tgz', commands=vaxignML_commands, default=True)
		
		
		correctInstall = 'VAXIGN_INSTALLED'
		installer = InstallHelper(VAXIGNML_DIC['name'], packageHome=cls.getVar(VAXIGNML_DIC['home']),
															packageVersion=VAXIGNML_DIC['version'])
		# Checks if the docker image exists, fails if not found
		installer.addCommand('docker pull e4ong1031/vaxign-ml:latest', correctInstall). \
			addPackage(env, dependencies=['docker'], default=default)


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
	def performEvaluations(cls, sequences, evalDics, jobs=1, browserData={}, outDir='/tmp', verbose=True):
			sDics, sWebDics = {k: v for k, v in evalDics.items() if v['software'] in STAND_SOFT}, \
												{k: v for k, v in evalDics.items() if v['software'] not in STAND_SOFT}

			epiDics = {}
			if len(sDics) > 0:
				epiDics.update(cls.performStandEvaluations(sequences, sDics, outDir))
			if len(sWebDics) > 0:
				epiDics.update(cls.performWebEvaluations(sequences, sWebDics, jobs, browserData, verbose))
			return epiDics

	@classmethod
	def performStandEvaluations(cls, sequences, evalDics, outDir):
		'''Generalize caller to the standalone evaluation functions.
    - sequences: dict with sequences in the form: {seqId: sequence}
    - evalDics: dictionary as {evalKey: {parameterName: parameterValue}}
    - jobs: int, number of jobs for parallelization
    Returns a dictionary of the form: {(evalKey, softwareName): [scores]}
    '''

		resultsDic, fKeys, outDic, sevalDics = {}, {}, {}, {}
		for evalKey, evalDic in evalDics.items():
			softName = evalDic['software']
			fKeys[(evalKey, softName)] = list(sequences.keys())
			smallEvalDic = evalDic.copy()
			del smallEvalDic['software']
			sevalDics[(evalKey, softName)] = smallEvalDic
			outDic[(evalKey, softName)] = os.path.join(outDir, f'{evalKey}_output.csv')
			resultsDic[(evalKey, softName)] = callIIITD(sequences, softName, smallEvalDic, outDic[(evalKey, softName)])

		# Check Subprocesses status and restart if failed
		pDic = {(evalKey, softName): None for (evalKey, softName) in resultsDic}
		while None in pDic.values():
			time.sleep(1)
			pDic = {}
			for (evalKey, softName), p in resultsDic.items():
				pDic[(evalKey, softName)] = p.poll()
				if pDic[(evalKey, softName)] == 1:
					resultsDic[(evalKey, softName)] = callIIITD(sequences, softName,
																											sevalDics[(evalKey, softName)], outDic[(evalKey, softName)])

		# Parse output results
		epiDics = {}
		for (evalKey, softName), res in resultsDic.items():
			fScores = parseIIITD(outDic[(evalKey, softName)], softName)
			allScores, i = [], 0
			for seqId in sequences:
				if seqId in fKeys[(evalKey, softName)]:
					allScores.append(fScores[i])
					i += 1
				else:
					allScores.append(0)

			epiDics[(evalKey, softName)] = allScores
		return epiDics

	@classmethod
	def performWebEvaluations(cls, sequences, evalDics, jobs=1, browserData={}, verbose=True):
		'''Generalize caller to the web evaluation functions.
    - sequences: dict with sequences in the form: {seqId: sequence}
    - evalDics: dictionary as {evalKey: {parameterName: parameterValue}}
    - jobs: int, number of jobs for parallelization
    Returns a dictionary of the form: {(evalKey, softwareName): [scores]}
    '''
		funcDic = {
			TOXINPRED: callToxinPred, ALGPRED2: callAlgPred2, TOXINPRED2: callToxinPred2,
			IL4PRED: callIL4pred, IL10PRED: callIL10pred, IFNEPITOPE: callIFNepitope,
		}

		# Create a pool of worker processes
		nJobs = len(evalDics) if len(evalDics) < jobs else jobs
		pool = multiprocessing.Pool(processes=nJobs)

		resultsDic, fKeys = {}, {}
		for evalKey, evalDic in evalDics.items():
			softName = evalDic['software']
			fKeys[(evalKey, softName)] = list(sequences.keys())
			smallEvalDic = evalDic.copy()
			del smallEvalDic['software']
			if softName in funcDic:
				resultsDic[(evalKey, softName)] = pool.apply_async(funcDic[softName],
																													 args=(sequences, browserData, smallEvalDic))

		if verbose:
			reportPoolStatus(resultsDic)

		pool.close()
		pool.join()

		epiDics = {}
		for (evalKey, softName), res in resultsDic.items():
			fScores = res.get()['Score'] if 'Score' in res.get() else []
			allScores, i = [], 0
			for seqId in sequences:
				if seqId in fKeys[(evalKey, softName)]:
					allScores.append(fScores[i])
					i += 1
				else:
					allScores.append(0)

			epiDics[(evalKey, softName)] = allScores
		return epiDics

	@classmethod
	def runVaxignML(cls, protocol, kwargs, cwd=None):
		""" Run vaxignML command from a given protocol.
		kwargs must contain:
		{"i": inputFasta, "o": outputDir, "-t": organism}
		other optional parameters are:
		{"s": modelPath, "p": numberProcessors}
		"""
		protId = protocol.getObjId()
		tmpDir = f'/tmp/VaxignML_{protId}'
		iFile, oDir = kwargs['i'], kwargs['o']
		kwargs['o'] = tmpDir

		# Vaxign-ML aborts when the requested process count reaches the number of
		# available cores, so the value is capped below that limit.
		if 'p' in kwargs:
			kwargs['p'] = max(1, min(int(kwargs['p']), multiprocessing.cpu_count() - 1))

		program = f"docker run --rm -v {iFile}:{iFile} -v {tmpDir}:{tmpDir} " \
							f"-v {tmpDir}/_FEATURE/PSORTB:/tmp/results " \
							"e4ong1031/vaxign-ml:latest python3.6 VaxignML.py "
		args = [f'-{k} {v}' for k,v in kwargs.items()]
		args = ' '.join(args)

		insistentRun(protocol, program, args, cwd=cwd, popen=True, stdout=subprocess.DEVNULL)
		# subprocess.check_call(program + args, shell=True, cwd=cwd, stdout=subprocess.DEVNULL)

		# Copying results dir with no-root user
		shutil.copytree(tmpDir, oDir)

		# Remove root results dir
		if os.path.exists(tmpDir):
			program = f"docker run --rm -v /:/mnt e4ong1031/vaxign-ml:latest rm -rf "
			args = f'/mnt/{tmpDir}'
			insistentRun(protocol, program, args, cwd=cwd, popen=True)
			# subprocess.check_call(program + args, shell=True, cwd=cwd, stdout=subprocess.DEVNULL)


	@classmethod
	def addEpiDopePackage(cls, env, default=True):
		# The environment MUST be created from the epidope.yml EpiDope ships
		# (Python 3.6, TensorFlow 1.13, ELMo/AllenNLP), installed directly at
		# EPIDOPE_HOME via '-p' so downstream code can rely on a stable
		# '<EPIDOPE_HOME>/bin/epidope' path regardless of where conda keeps its
		# named environments.
		installer = InstallHelper(EPIDOPE_DIC['name'], packageHome=cls.getVar(EPIDOPE_DIC['home']),
															packageVersion=EPIDOPE_DIC['version'])

		home = cls.getVar(EPIDOPE_DIC['home'])
		installer.addCommand(
			f'git clone --depth 1 {EPIDOPE_UPSTREAM_URL}.git',
			'EPIDOPE_CLONED'
		).addCommand(
			f'cd EpiDope && conda env create -f epidope.yml -p {home}',
			'EPIDOPE_INSTALLED'
		).addPackage(env, dependencies=['conda', 'git'], default=default)

	@classmethod
	def addStackGlyEmbedPackage(cls, env, default=True):
		home = cls.getVar(STACKGLYEMBED_DIC['home'])
		esmModelName = cls.getVar(STACKGLYEMBED_DIC['esm_model_name'])
		t5ModelName = STACKGLYEMBED_DEFAULT_T5_MODEL_NAME
		proteinbertDir = cls.getStackGlyEmbedProteinBertDir()
		hfCacheDir = cls.getStackGlyEmbedHfCacheDir()

		# ProteinBERT's load_pretrained_model() defaults to an INTERACTIVE
		# confirmation prompt before downloading (validate_downloading=True
		# by default, confirmed by reading nadavbra/protein_bert's own
		# existing_model_loading.py) -- would hang a non-interactive
		# 'scipion3 installb' forever waiting on stdin. Explicitly disabled
		# here, install-time only (the runtime call in predict_local.py
		# keeps download_model_dump_if_not_exists=False, so a MISSING dump
		# still fails loudly instead of silently reaching the network).
		#
		# Real bug found+fixed via an actual 'scipion3 installb'
		# run: load_pretrained_model()'s own DEFAULT_REMOTE_MODEL_DUMP_URL
		# is an FTP URL (ftp://ftp.cs.huji.ac.il/...) -- confirmed a real,
		# reproducible network-level connection timeout on port 21 in this
		# environment (many networks/firewalls block outbound FTP). Fixed
		# by pre-downloading the exact same file over HTTPS from a real,
		# verified mirror (github.com/nadavbra/proteinbert_data_files,
		# same author as the upstream repo, byte-identical 192MB file
		# confirmed via a real HTTPS request) directly to the path
		# load_pretrained_model() itself checks first
		# ('{proteinbertDir}/default.pkl', confirmed reading its source:
		# DEFAULT_LOCAL_MODEL_DUMP_FILE_NAME='default.pkl') -- since that
		# exact file already exists, load_pretrained_model() skips its own
		# download entirely and just loads it, so this also works
		# unmodified for anyone whose network allows the FTP path (the curl
		# step becomes a harmless no-op once cached, and the call below is
		# kept as-is as a fallback for the (empty) case where it doesn't).
		proteinbertMirrorUrl = (
			"https://github.com/nadavbra/proteinbert_data_files/raw/master/"
			"epoch_92400_sample_23500000.pkl"
		)
		downloadProteinBertCmd = (
			f"curl -fsSL -o {proteinbertDir}/default.pkl {proteinbertMirrorUrl}"
		)
		primeProteinBertCmd = (
			f"{downloadProteinBertCmd} && "
			f"python -c \"from proteinbert import load_pretrained_model; "
			f"load_pretrained_model(local_model_dump_dir='{proteinbertDir}', "
			f"download_model_dump_if_not_exists=True, validate_downloading=False)\""
		)
		# ESM-2 and ProtT5 both resolve via transformers' HF Hub cache
		# (HF_HOME redirected to a plugin-local dir instead of the user's
		# global ~/.cache/huggingface): calling from_pretrained() once here
		# with network available caches them for the fully offline runtime
		# call in predict_local.py (which already sets
		# HF_HUB_OFFLINE=1/TRANSFORMERS_OFFLINE=1).
		primeHfCmd = (
			f"HF_HOME={hfCacheDir} python -c \""
			f"from transformers import AutoTokenizer, EsmModel, T5EncoderModel, T5Tokenizer; "
			f"AutoTokenizer.from_pretrained('{esmModelName}'); "
			f"EsmModel.from_pretrained('{esmModelName}'); "
			f"T5Tokenizer.from_pretrained('{t5ModelName}', do_lower_case=False); "
			f"T5EncoderModel.from_pretrained('{t5ModelName}')\""
		)

		installer = InstallHelper(STACKGLYEMBED_DIC['name'], packageHome=home,
								  packageVersion=STACKGLYEMBED_DIC['version'])

		# requirements.txt IS used (patched in place first), following the
		# same real precedent already established elsewhere in the
		# scipion-chem/scipion-em ecosystem for an upstream dependency file
		# that needs a small fix before it is installable as-is --
		# confirmed by reading the real source of scipion-chem-bepipred
		# ('sed -i 's/^torch==/#torch==/g' requirements.txt') and
		# scipion-chem-iedb ('sed -i 's/^torch==1.12.0/torch/g'
		# requirements.txt'), both `sed`-patch a bad line in-place and THEN
		# install from the file, rather than abandoning it for a hardcoded
		# list. Verified against the real upstream file: it is
		# StackGlyEmbed's own TRAINING-framework file (Python 3.9.19), so
		# patched the same way here:
		#   - 'pytorch==2.2.2' commented out (not a real PyPI name as-is;
		#     also too old for modern transformers regardless -- real torch
		#     install happens as its own supplementary step below, from the
		#     CPU-only wheel index)
		#   - 'pickle5==0.0.11' commented out (Python<3.8-only backport,
		#     this env is 3.10; also not imported by predict_local.py)
		#   - PyQt5/shap/IPython/imblearn commented out (real training-only
		#     extras, confirmed absent from predict_local.py's own
		#     imports/pickle usage -- PyQt5 in particular is a heavy GUI/
		#     X11 dependency, the single biggest real friction risk in a
		#     headless conda env)
		# 'xgboost==' is NOT commented out, unlike the others above (real
		# bug found+fixed via an actual 'scipion3 test' run,
		# AFTER everything else already worked): predict_local.py never
		# writes 'import xgboost' directly, but it unpickles real
		# 'XGB_base_layer_*.sav' files (predict.py::_base_layer_predictions)
		# whose stored objects are actual xgboost.XGBClassifier instances --
		# unpickling requires the defining class's module to be importable
		# even with no direct import statement in the consuming script.
		# Confirmed real: 'ModuleNotFoundError: No module named xgboost'
		# at that exact unpickling call. This is a genuine runtime
		# dependency, not a training-only extra like its neighbors.
		# numpy/pandas/scikit-learn/matplotlib/xgboost keep their real
		# upstream pins this way, so a future version bump upstream is
		# tracked automatically -- same benefit as netcleave/iapred/
		# scannet/discotope. transformers/tensorflow/protein-bert are NOT
		# in this file at all (it predates a prediction path): the
		# README's own separate "Prediction / Prerequisites" section
		# documents those instead, installed as a supplementary step, same
		# as IApred's imbalanced-learn/matplotlib/seaborn gap. PyPI package
		# name is 'protein-bert' (hyphen), not 'proteinbert' -- that is
		# only the importable module name.
		patchRequirementsCmd = (
			"sed -i "
			"-e '/^pytorch==/s/^/#/' "
			"-e '/^pickle5==/s/^/#/' "
			"-e '/^PyQt5==/s/^/#/' "
			"-e '/^shap==/s/^/#/' "
			"-e '/^IPython==/s/^/#/' "
			"-e '/^imblearn==/s/^/#/' "
			"requirements.txt"
		)

		# Clone BEFORE creating the conda env (real bug found+fixed via an
		# actual 'scipion3 installb' run -- see
		# netcleave/__init__.py::addNetCleavePackage for the full
		# explanation: 'getCondaEnvCommand' leaves its own completion
		# marker inside 'home', which then blocks a subsequent 'git clone'
		# into that same now-nonempty directory).
		installer.addCommand(
			f"git clone --depth 1 {STACKGLYEMBED_UPSTREAM_URL} {home}",
			'STACKGLYEMBED_CLONED'
		).getCondaEnvCommand(
			STACKGLYEMBED_DIC['name'], binaryVersion=STACKGLYEMBED_DIC['version'], pythonVersion='3.10'
		).addCommand(
			# 'transformers<5' pinned explicitly (real bug found+fixed via
			# an actual 'scipion3 installb' run -- same class
			# of bug as scipion-chem-tmbed's own transformers pin): an
			# unpinned 'pip install transformers' pulls the latest 5.x
			# line, which requires 'torch>=2.4' -- confirmed real:
			# 'EsmModel requires the PyTorch library but it was not found'
			# even with a real torch installed, if it's too old.
			# 'sentencepiece' also added explicitly (real, separate bug
			# found the same run): T5Tokenizer needs it and neither
			# StackGlyEmbed nor transformers declare it as a hard
			# dependency (it is an optional extra upstream).
			#
			# torch installed SEPARATELY from the PyTorch CPU-only wheel
			# index, plus an explicit purge of any nvidia-*/cuda-*/triton
			# packages afterwards (real bug found+fixed, THIRD
			# issue found installing this same run): a plain 'pip install
			# torch' on Linux resolves the CUDA-enabled build by default
			# regardless of whether a GPU/CUDA driver is actually present
			# -- on this GPU-less machine that build caused a real,
			# reproduced segfault (SIGSEGV, confirmed via 'python -X
			# faulthandler') deep inside 'torch._dynamo' the moment
			# 'transformers' imports it (masking_utils -> generation.utils
			# -> torch._dynamo). Installing from
			# 'download.pytorch.org/whl/cpu' alone was NOT enough to fully
			# fix it once stray nvidia-*/cuda-bindings/triton packages were
			# already present (confirmed real: the segfault persisted after
			# switching torch's own build to '+cpu' until those extra
			# packages -- pulled in as transitive deps of the CUDA build,
			# e.g. by 'tensorflow[and-cuda]'-style bundling -- were also
			# uninstalled) -- purging them explicitly as a real, verified
			# fix, not just switching torch's own wheel.
			#
			# 'tensorflow<2.16' pinned (real, separate bug found the same
			# run, AFTER the segfault was fixed): protein-bert's own
			# model_generation.py calls
			# 'keras.optimizers.Adam(lr=..., ...)', the legacy Keras 2 API
			# -- Keras 3 (bundled by default from TF 2.16 onwards) removed
			# the 'lr' alias entirely ('ValueError: Argument(s) not
			# recognized: {\'lr\': ...}', confirmed real running an actual
			# prediction). Pinning tensorflow keeps Keras 2 as the default
			# backend, matching what protein-bert's own code was written
			# against.
			f"cd {home} && {patchRequirementsCmd} && "
			f"{cls.getEnvActivationCommand(STACKGLYEMBED_DIC)} && "
			"pip install -r requirements.txt && "
			"pip install --index-url https://download.pytorch.org/whl/cpu torch && "
			"pip install 'tensorflow<2.16' 'transformers<5' protein-bert sentencepiece && "
			"pip uninstall -y cuda-bindings cuda-pathfinder cuda-toolkit nvidia-cublas "
			"nvidia-cuda-cupti nvidia-cuda-nvrtc nvidia-cuda-runtime nvidia-cudnn-cu13 "
			"nvidia-cufft nvidia-cufile nvidia-curand nvidia-cusolver nvidia-cusparse "
			"nvidia-cusparselt-cu13 nvidia-nccl-cu13 nvidia-nvjitlink nvidia-nvshmem-cu13 "
			"nvidia-nvtx triton || true",
			'STACKGLYEMBED_DEPS_INSTALLED'
		).addCommand(
			f"mkdir -p {proteinbertDir} {hfCacheDir} && "
			f"{cls.getEnvActivationCommand(STACKGLYEMBED_DIC)} && "
			f"{primeProteinBertCmd} && {primeHfCmd}",
			'STACKGLYEMBED_INSTALLED'
		).addPackage(env, dependencies=['conda', 'git'], default=default)

	@classmethod
	def addNetCleavePackage(cls, env, default=True):
		home = cls.getVar(NETCLEAVE_DIC['home'])

		installer = InstallHelper(NETCLEAVE_DIC['name'], packageHome=home,
															packageVersion=NETCLEAVE_DIC['version'])

		# Installed from NetCleave's own requirements.txt, a strict superset of
		# what this protocol needs, so a future upstream pin change is picked
		# up automatically.
		#
		# Clone BEFORE creating the conda env: 'InstallHelper.addCommand'
		# writes its completion marker directly inside 'packageHome', and
		# 'getCondaEnvCommand' is itself one such call, so running it first
		# leaves a marker in 'home' and the subsequent 'git clone' fails with
		# 'destination path already exists and is not an empty directory'.
		installer.addCommand(
			f"git clone --depth 1 {NETCLEAVE_UPSTREAM_URL} {home}",
			'NETCLEAVE_CLONED'
		).getCondaEnvCommand(
			NETCLEAVE_DIC['name'], binaryVersion=NETCLEAVE_DIC['version'], pythonVersion='3.10'
		).addCommand(
			f"{cls.getEnvActivationCommand(NETCLEAVE_DIC)} && "
			f"cd {home} && pip install -r requirements.txt",
			'NETCLEAVE_INSTALLED'
		).addPackage(env, dependencies=['conda', 'git'], default=default)

	@classmethod
	def addIApredPackage(cls, env, default=True):
		home = cls.getVar(IAPRED_DIC['home'])

		installer = InstallHelper(IAPRED_DIC['name'], packageHome=home,
															packageVersion=IAPRED_DIC['version'])

		# Installed from IApred's own requirements.txt so a future upstream pin
		# change is picked up automatically, but that file genuinely
		# under-declares its runtime dependencies: 'imbalanced-learn' is
		# commented out in it and matplotlib/seaborn are not listed at all,
		# even though its own functions.py imports all three. Those three are
		# installed as a supplementary step, not as a replacement.
		#
		# Clone BEFORE creating the conda env, for the same reason as ScanNet:
		# 'getCondaEnvCommand' leaves its own completion marker inside 'home',
		# which then blocks a subsequent 'git clone' into that same
		# now-nonempty directory.
		installer.addCommand(
			f"git clone --depth 1 {IAPRED_UPSTREAM_URL} {home}",
			'IAPRED_CLONED'
		).getCondaEnvCommand(
			IAPRED_DIC['name'], binaryVersion=IAPRED_DIC['version'], pythonVersion='3.10'
		).addCommand(
			f"{cls.getEnvActivationCommand(IAPRED_DIC)} && "
			f"cd {home} && pip install -r requirements.txt && "
			"pip install imbalanced-learn matplotlib seaborn",
			'IAPRED_INSTALLED'
		).addPackage(env, dependencies=['conda', 'git'], default=default)

	@classmethod
	def addScanNetPackage(cls, env, default=True):
		home = cls.getVar(SCANNET_DIC['home'])

		installer = InstallHelper(SCANNET_DIC['name'], packageHome=home,
															packageVersion=SCANNET_DIC['version'])

		# Clone BEFORE creating the conda env: 'getCondaEnvCommand' leaves its
		# own completion marker inside 'home', which then blocks a subsequent
		# 'git clone' into that same now-nonempty directory.
		installer.addCommand(
			f"git clone --depth 1 {SCANNET_UPSTREAM_URL} {home}",
			'SCANNET_CLONED'
		).getCondaEnvCommand(
			SCANNET_DIC['name'], binaryVersion=SCANNET_DIC['version'], pythonVersion='3.6.12'
		).addCommand(
			f"{cls.getEnvActivationCommand(SCANNET_DIC)} && "
			f"cd {home} && pip install -r requirements.txt",
			'SCANNET_INSTALLED'
		).addPackage(env, dependencies=['conda', 'git'], default=default)

	@classmethod
	def addDiscoTopePackage(cls, env, default=True):
		# Python is pinned to 3.14, per the upstream project's own README, not
		# the stale 'Python :: 3.9' classifier in its setup.py.
		home = cls.getVar(DISCOTOPE_DIC['home'])
		weightsCacheDir = cls.getDiscoTopeWeightsCacheDir()

		installer = InstallHelper(DISCOTOPE_DIC['name'], packageHome=home,
															packageVersion=DISCOTOPE_DIC['version'])

		installer.addCommand(
			f"git clone --depth 1 {DISCOTOPE_UPSTREAM_URL} {home}",
			'DISCOTOPE_CLONED'
		).getCondaEnvCommand(
			DISCOTOPE_DIC['name'], binaryVersion=DISCOTOPE_DIC['version'], pythonVersion='3.14'
		).addCommand(
			# 'unzip' installed as a conda package INSIDE this env, not relied
			# upon as a system binary: 'conda activate' replaces PATH entirely.
			f"{cls.getEnvActivationCommand(DISCOTOPE_DIC)} && "
			"conda install -y -c conda-forge unzip && "
			f"cd {home} && pip install -r requirements.txt && pip install . && unzip -q models.zip",
			'DISCOTOPE_DEPS_INSTALLED'
		).addCommand(
			f"mkdir -p {weightsCacheDir} && "
			f"{cls.getEnvActivationCommand(DISCOTOPE_DIC)} && "
			f"TORCH_HOME={weightsCacheDir} python -c "
			f"\"from discotope3.esm.pretrained import esm_if1_gvp4_t16_142M_UR50; "
			f"esm_if1_gvp4_t16_142M_UR50()\"",
			'DISCOTOPE_INSTALLED'
		).addPackage(env, dependencies=['conda', 'git', 'unzip'], default=default)

	@classmethod
	def addTMbedPackage(cls, env, default=True):
		home = cls.getVar(TMBED_DIC['home'])
		modelDir = cls.getTMbedModelDir()
		primeT5Cmd = (
			f"python -c \"from tmbed.embed import T5Encoder; "
			f"T5Encoder(model_path='{modelDir}', use_gpu=False)\""
		)

		installer = InstallHelper(TMBED_DIC['name'], packageHome=home,
															packageVersion=TMBED_DIC['version'])

		# TMbed is not published on PyPI: installed from its tagged GitHub
		# release. transformers is pinned <5: TMbed's own embed.py still
		# calls T5Tokenizer.batch_encode_plus, removed in transformers 5.x.
		installer.getCondaEnvCommand(
			TMBED_DIC['name'], binaryVersion=TMBED_DIC['version'], pythonVersion='3.10'
		).addCommand(
			f"{cls.getEnvActivationCommand(TMBED_DIC)} && "
			f"pip install git+{TMBED_DOWNLOAD_URL}.git@v{TMBED_DIC['version']} && "
			"pip install 'transformers<5' protobuf tiktoken",
			'TMBED_INSTALLED'
		).addCommand(
			f"mkdir -p {modelDir} && {cls.getEnvActivationCommand(TMBED_DIC)} && {primeT5Cmd}",
			'TMBED_WEIGHTS_CACHED'
		).addPackage(env, dependencies=['conda', 'git'], default=default)

	@classmethod
	def validateStackGlyEmbedInstallation(cls):
		errors = []

		modelsDir = cls.getStackGlyEmbedModelsDir()
		if not os.path.isdir(modelsDir):
			errors.append(f"Could not find the 'prediction/' folder under STACKGLYEMBED_HOME: '{modelsDir}'.")

		proteinbertDump = os.path.join(cls.getStackGlyEmbedProteinBertDir(), 'default.pkl')
		if not os.path.isfile(proteinbertDump):
			errors.append(f"ProteinBERT model dump not found: '{proteinbertDump}'.")

		hfCacheDir = cls.getStackGlyEmbedHfCacheDir()
		if not os.path.isdir(hfCacheDir) or not os.listdir(hfCacheDir):
			errors.append(f"ESM-2/ProtT5 HuggingFace cache is empty: '{hfCacheDir}'.")

		if not errors and not cls.checkCallEnv(
				STACKGLYEMBED_DIC, 'import torch, tensorflow, transformers, proteinbert'):
			errors.append("Activation of the StackGlyEmbed conda environment failed.")

		if errors:
			errors.append(STACKGLYEMBED_NOINSTALL_WARNING)
		return errors

	@classmethod
	def validateNetCleaveInstallation(cls):
		errors = []
		scriptPath = cls.getNetCleaveScriptPath()
		if not os.path.isfile(scriptPath):
			errors.append(f"Could not find 'NetCleave.py' under NETCLEAVE_HOME: '{cls.getNetCleaveDir()}'.")
		elif not cls.checkCallEnv(NETCLEAVE_DIC, 'import tensorflow'):
			errors.append("Activation of the NetCleave conda environment failed.")
		if errors:
			errors.append(NETCLEAVE_NOINSTALL_WARNING)
		return errors

	@classmethod
	def validateIApredInstallation(cls):
		errors = []
		scriptPath = cls.getIApredScriptPath()
		home = cls.getIApredDir()
		if not os.path.isfile(scriptPath):
			errors.append(f"Could not find '{cls.getVar(IAPRED_DIC['script_name'])}' under IAPRED_HOME: '{home}'.")
		elif not os.path.isdir(os.path.join(home, 'models')):
			errors.append(f"Could not find the 'models/' folder under IAPRED_HOME: '{home}'.")
		elif not cls.checkCallEnv(IAPRED_DIC, 'import sklearn'):
			errors.append("Activation of the IApred conda environment failed.")
		if errors:
			errors.append(IAPRED_NOINSTALL_WARNING)
		return errors

	@classmethod
	def validateEpiDopeInstallation(cls):
		errors = []
		binaryPath = cls.getEpiDopeBinaryPath()
		if not os.path.isfile(binaryPath):
			errors.append(f"Could not find the 'epidope' entry point under EPIDOPE_HOME: '{binaryPath}'.")
		if errors:
			errors.append(EPIDOPE_NOINSTALL_WARNING)
		return errors

	@classmethod
	def validateScanNetInstallation(cls):
		errors = []
		home = cls.getScanNetDir()
		if not os.path.isfile(os.path.join(home, 'predict_bindingsites.py')):
			errors.append(f"Could not find 'predict_bindingsites.py' under SCANNET_HOME: '{home}'.")
		elif not cls.checkCallEnv(SCANNET_DIC, 'import numpy'):
			errors.append("Activation of the ScanNet conda environment failed.")
		if errors:
			errors.append(SCANNET_NOINSTALL_WARNING)
		return errors

	@classmethod
	def validateDiscoTopeInstallation(cls):
		errors = []
		mainScript = cls.getDiscoTopeMainScriptPath()
		modelsDir = cls.getDiscoTopeModelsDir()
		if not os.path.isfile(mainScript):
			errors.append(f"Could not find 'discotope3/main.py' under DISCOTOPE_HOME: '{cls.getDiscoTopeDir()}'.")
		elif not os.path.isdir(modelsDir):
			errors.append(f"Could not find the unzipped 'models/' folder under DISCOTOPE_HOME: '{modelsDir}'.")
		elif not cls.checkCallEnv(DISCOTOPE_DIC, 'import discotope3'):
			errors.append("Activation of the DiscoTope-3.0 conda environment failed.")
		if errors:
			errors.append(DISCOTOPE_NOINSTALL_WARNING)
		return errors

	@classmethod
	def validateTMbedInstallation(cls):
		errors = []
		modelDir = cls.getTMbedModelDir()
		missing = [fn for fn in TMBED_T5_MODEL_REQUIRED_FILES if not os.path.isfile(os.path.join(modelDir, fn))]
		if not any(os.path.isfile(os.path.join(modelDir, fn)) for fn in TMBED_T5_MODEL_WEIGHT_FILE_ALTERNATIVES):
			missing.append(f"one of {TMBED_T5_MODEL_WEIGHT_FILE_ALTERNATIVES}")
		if not any(os.path.isfile(os.path.join(modelDir, fn)) for fn in TMBED_T5_MODEL_TOKENIZER_FILE_ALTERNATIVES):
			missing.append(f"one of {TMBED_T5_MODEL_TOKENIZER_FILE_ALTERNATIVES}")
		if missing:
			errors.append(f"TMBED_HOME ('{cls.getTMbedDir()}') is missing expected ProtT5 file(s) "
										f"under '{modelDir}': {missing}.")
		elif not cls.checkCallEnv(TMBED_DIC, 'tmbed --help', isModuleImport=False):
			errors.append("Activation of the TMbed conda environment failed.")
		if errors:
			errors.append(TMBED_NOINSTALL_WARNING)
		return errors

	@classmethod
	def validateSignalPInstallation(cls):
		errors = []
		pythonBin = cls.getVar(SIGNALP_DIC['python_bin'])
		binaryPath = cls.getSignalPBinaryPath()
		if not pythonBin or not os.path.isfile(pythonBin):
			errors.append(f"SIGNALP_PYTHON_BIN is not set or does not exist: '{pythonBin}'.")
		elif not binaryPath or not os.path.isfile(binaryPath):
			errors.append(f"Could not find the local SignalP-6.0 binary at '{binaryPath}'.")
		modelDir = cls.getVar(SIGNALP_DIC['model_dir'])
		if not modelDir or not os.path.isdir(os.path.join(modelDir or '', 'sequential_models_signalp6')):
			errors.append(f"Could not find 'sequential_models_signalp6/' under SIGNALP_MODEL_DIR: '{modelDir}'.")
		if errors:
			errors.append(SIGNALP_NOINSTALL_WARNING)
		return errors

	@classmethod
	def checkCallEnv(cls, packageDic, checkCmd, isModuleImport=True):
		actCommand = cls.getVar(packageDic['activation'])
		pyCmd = f'python -c "{checkCmd}"' if isModuleImport else checkCmd
		try:
			if 'conda' in actCommand and 'shell.bash hook' not in actCommand:
				actCommand = f'{cls.getCondaActivationCmd()}{actCommand}'
			subprocess.check_output(f'{actCommand} && {pyCmd}', shell=True)
			return True
		except subprocess.CalledProcessError:
			return False

	# ---------------------------------- Utils (EpiDope/ScanNet/DiscoTope/TMbed/SignalP) --

	@classmethod
	def getStackGlyEmbedDir(cls):
		return cls.getVar(STACKGLYEMBED_DIC['home'])

	@classmethod
	def getStackGlyEmbedModelsDir(cls):
		return os.path.join(cls.getStackGlyEmbedDir(), 'prediction')

	@classmethod
	def getStackGlyEmbedProteinBertDir(cls):
		return os.path.join(cls.getStackGlyEmbedDir(), '.proteinbert_models')

	@classmethod
	def getStackGlyEmbedHfCacheDir(cls):
		return os.path.join(cls.getStackGlyEmbedDir(), '.hf_cache')

	@classmethod
	def getStackGlyEmbedT5ModelName(cls):
		return STACKGLYEMBED_DEFAULT_T5_MODEL_NAME

	@classmethod
	def getNetCleaveDir(cls):
		return cls.getVar(NETCLEAVE_DIC['home'])

	@classmethod
	def getNetCleaveScriptPath(cls):
		return os.path.join(cls.getNetCleaveDir(), 'NetCleave.py')

	@classmethod
	def getIApredDir(cls):
		return cls.getVar(IAPRED_DIC['home'])

	@classmethod
	def getIApredScriptPath(cls):
		return os.path.join(cls.getIApredDir(), cls.getVar(IAPRED_DIC['script_name']))

	@classmethod
	def getEpiDopeBinaryPath(cls):
		return os.path.join(cls.getVar(EPIDOPE_DIC['home']), 'bin', 'epidope')

	@classmethod
	def getScanNetDir(cls):
		return cls.getVar(SCANNET_DIC['home'])

	@classmethod
	def getDiscoTopeDir(cls):
		return cls.getVar(DISCOTOPE_DIC['home'])

	@classmethod
	def getDiscoTopeMainScriptPath(cls):
		return os.path.join(cls.getDiscoTopeDir(), 'discotope3', 'main.py')

	@classmethod
	def getDiscoTopeModelsDir(cls):
		return os.path.join(cls.getDiscoTopeDir(), 'models')

	@classmethod
	def getDiscoTopeWeightsCacheDir(cls):
		return os.path.join(cls.getDiscoTopeDir(), '.torch_cache')

	@classmethod
	def getTMbedDir(cls):
		return cls.getVar(TMBED_DIC['home'])

	@classmethod
	def getTMbedModelDir(cls):
		return os.path.join(cls.getTMbedDir(), 'prott5_weights')

	@classmethod
	def getSignalPBinaryPath(cls):
		pythonBin = cls.getVar(SIGNALP_DIC['python_bin'])
		if not pythonBin:
			return None
		return os.path.join(os.path.dirname(pythonBin), cls.getVar(SIGNALP_DIC['binary_name']))

	# ---------------------------------- Protocol functions (EpiDope/ScanNet/DiscoTope/TMbed) --

	@classmethod
	def runStackGlyEmbed(cls, protocol, scriptPath, args, cwd=None):
		# HF_HOME is set on os.environ directly and no explicit 'env' override
		# is passed, the same pattern DiscoTope uses for TORCH_HOME: passing
		# env=cls.getEnviron() would risk silently dropping HF_HOME rather
		# than failing loudly.
		os.environ['HF_HOME'] = cls.getStackGlyEmbedHfCacheDir()
		activation = cls.getVar(STACKGLYEMBED_DIC['activation'])
		fullProgram = f'{activation} && python {scriptPath}'
		protocol.runJob(fullProgram, args, cwd=cwd)

	@classmethod
	def runNetCleave(cls, protocol, args, cwd=None):
		# NetCleave.py resolves its bundled model path relative to the process
		# cwd, not to its own location: the caller must pass
		# cwd=getNetCleaveDir().
		activation = cls.getVar(NETCLEAVE_DIC['activation'])
		fullProgram = f'{activation} && python {cls.getNetCleaveScriptPath()}'
		protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)

	@classmethod
	def runIApred(cls, protocol, args, cwd=None):
		# IApred.py resolves its own 'models_folder = "models"' relative to the
		# process cwd, not to the script's location: the caller must pass
		# cwd=getIApredDir().
		activation = cls.getVar(IAPRED_DIC['activation'])
		fullProgram = f'{activation} && python {cls.getIApredScriptPath()}'
		protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)

	@classmethod
	def runEpiDope(cls, protocol, args, cwd=None):
		# EpiDope's 'epidope' entry point is a self-contained shim whose shebang
		# already points at its own environment's interpreter: it must NOT be
		# wrapped in a conda activation command, which produced intermittent
		# spurious failures with piped stdout/stderr.
		protocol.runJob(cls.getEpiDopeBinaryPath(), args, env=cls.getEnviron(), cwd=cwd)

	@classmethod
	def runScanNet(cls, protocol, args, cwd=None):
		# ScanNet resolves its own 'models/' path relative to the process cwd,
		# not the script's own location -- caller must pass cwd=getScanNetDir().
		activation = cls.getVar(SCANNET_DIC['activation'])
		fullProgram = f'{activation} && python predict_bindingsites.py'
		protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)

	@classmethod
	def runDiscoTope(cls, protocol, args, cwd=None):
		# runJob's 'env' kwarg expects a pyworkflow Environ object, not a plain
		# dict -- TORCH_HOME is set on os.environ directly instead, which a
		# subprocess launched with no explicit 'env' override inherits.
		os.environ['TORCH_HOME'] = cls.getDiscoTopeWeightsCacheDir()
		activation = cls.getVar(DISCOTOPE_DIC['activation'])
		fullProgram = f'{activation} && python {cls.getDiscoTopeMainScriptPath()}'
		protocol.runJob(fullProgram, args, cwd=cwd)

	@classmethod
	def runTMbed(cls, protocol, args, cwd=None):
		""" Run TMbed's 'predict' subcommand through the dedicated conda env
		(TMbed requires an exact torch/transformers/sentencepiece stack). """
		activation = cls.getVar(TMBED_DIC['activation'])
		fullProgram = f'{activation} && tmbed'
		protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)

	# ---------------------------------- Utils functions-----------------------
	@classmethod
	def getBrowserData(cls):
		return {'name': cls.getVar(IIITDW_DIC['browser']), 'path': cls.getVar(IIITDW_DIC['browserPath'])}

	@classmethod
	def getEnvScriptsPath(cls, envDic, software, scriptName):
		return pwchemPlugin.getEnvPath(envDic, f'lib/python{envDic["python"]}/site-packages/{software}'
																					 f'/python_scripts/{scriptName}.py')

	@classmethod
	def getFixScripts(cls, envDic, fixDic):
		cmds = []
		for software, repDic in fixDic.items():
			for fileName, repList in repDic.items():
				scriptFile = cls.getEnvScriptsPath(envDic, software.lower(), fileName)
				for repPair in repList:
					cmds.append(getReplaceCommand(scriptFile, repPair[0], repPair[1]))
		return ' && '.join(cmds)





