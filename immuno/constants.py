# coding: latin-1
# **************************************************************************
# *
# * Authors:  Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
# *
# * Biocomputing Unit, CNB-CSIC
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

# Common constants
DEFAULT_VERSION = '1.0'

# Package dictionaries
IIITDW_DIC = {'name': 'IIITD_WEB',    'version': '3.0',
              'home': 'IIITD_WEB_HOME', 'activation': 'IIITD_WEB_ACTIVATION_CMD',
              'browser': 'IIITD_BROWSER', 'browserPath': 'IIITD_BROWSER_PATH'}

IIITD_DIC = {'name': 'IIITD',    'version': '3.0', 'python': '3.7',
             'home': 'IIITD_HOME', 'activation': 'IIITD_ACTIVATION_CMD'}

IL6PRED_DIC = {'name': 'IL6PRED',    'version': '1.1', 'python': '3.7',
               'home': 'IL6PRED_HOME', 'activation': 'IL6PRED_ACTIVATION_CMD'}

# toxinpred2/algpred2 pinned to the exact versions this project already
# validates against in its own standalone plugins (scipion-chem-toxinpred2,
# scipion-chem-algpred): an unpinned 'pip install' installs whatever is
# latest on PyPI at install time, which can silently drift from what was
# actually tested here. The rest are left unpinned (out of scope for this
# change -- no verified version pin sourced for them yet).
IIITD_PACKAGES = ["toxinpred3", "toxinpred2==1.1", "algpred2==1.4", "ifnepitope2", "il5pred", "il13pred", "clbtope"]

# Standalone packages
TOXINPRED, TOXINPRED2, IFNEPITOPE, ALGPRED2, IL5PRED, IL6PRED, IL13PRED = 'ToxinPred3', 'ToxinPred2', 'IFNepitope2', \
                                                                          'AlgPred2',  'IL5pred', 'IL6pred', 'IL13pred'
STAND_SOFT = [TOXINPRED, TOXINPRED2, IFNEPITOPE, ALGPRED2, IL5PRED, IL6PRED, IL13PRED]


OLD_AAC = '''def aac_comp(file,out):
    std = list("ACDEFGHIKLMNPQRSTVWY")
    f = open(out, "w")
    sys.stdout = f
    df = pd.DataFrame(file)
    zz = df.iloc[:,0]
    for j in zz:
        for i in std:
            count = 0
            for k in j:
                temp1 = k
                if temp1 == i:
                    count += 1
                composition = (count/len(j))*100
            print("%.2f"%composition, end = ",")
        print("")
    f.truncate()'''

NEW_AAC = '''def aac_comp(file,out):
    std = list("ACDEFGHIKLMNPQRSTVWY")
    df1 = pd.DataFrame(file, columns=["Seq"])
    dd = []
    for j in df1["Seq"]:
        cc = []
        for i in std:
            count = 0
            for k in j:
                temp1 = k
                if temp1 == i:
                    count += 1
                composition = (count/len(j))*100
            cc.append(composition)
        dd.append(cc)
    df2 = pd.DataFrame(dd)
    head = []
    for mm in std:
        head.append("AAC_"+mm)
    df2.columns = head
    df2.to_csv(out, index=None, header=False)'''

IIITD_FIXES = {TOXINPRED: {'toxinpred3': [("np.loadtxt(file_name, delimiter=',')", "np.loadtxt(file_name, delimiter=',', ndmin=2)"),
                                          ("np.loadtxt(file_name3, delimiter=',')", "np.loadtxt(file_name3, delimiter=',', ndmin=2)")]},
               # toxinpred2==1.1's own Hybrid-mode BLAST_processor(blast_result, blast_processed,
               # name1) references the OLD pre-refactor variable name 'seqid' instead of its own
               # 'name1' parameter in one branch (NameError: name 'seqid' is not defined) -- confirmed
               # real via a direct 'toxinpred2 -m 2' run, not assumed. Genuine upstream bug in the
               # PyPI package itself, patched the same way as the other known IIITD bugs.
               TOXINPRED2: {'toxinpred2': [("np.loadtxt(file_name, delimiter=',')", "np.loadtxt(file_name, delimiter=',', ndmin=2)"),
                                           ("for j in seqid:", "for j in name1:")]},
               IFNEPITOPE: {'ifnepitope2': [('blastdb1 = nf_path + "\/\.\.\/blast_db\/human_db"', ''),
                            ('blastdb2 = nf_path + "\/\.\.\/blast_db\/mouse_db"', 'blastdb = nf_path + "\/\.\.\/blast_db"')]},
               # NOTE: getReplaceCommand (pwchem.utils) builds a plain 'sed -z s/.../.../ ' command
               # with NO regex-metacharacter escaping. sed's basic regex treats '[' ']' as a bracket
               # expression (character CLASS), not literal brackets -- an unescaped '[df2.Subject==i]'
               # matches none of the real file's lines, so the substitution silently no-ops (sed
               # doesn't error on a non-matching pattern). Confirmed real: a direct 'algpred2 -m 2' run
               # crashed with "'DataFrame' object has no attribute 'concat'" -- the exact bug this fix
               # is supposed to patch -- because it was never actually applied. '[' and ']' are escaped
               # below as '\[' '\[' so sed matches them literally; '.' and '(' ')' are left as-is
               # (harmless here: '.' as wildcard still matches the literal char it stands for, and '('
               # ')' are literal, not special, in sed's basic regex).
               ALGPRED2: {'algpred2': [("np.loadtxt(file_name, delimiter=',')", "np.loadtxt(file_name, delimiter=',', ndmin=2)"),
                          ("df3 = df3.concat(df2.loc\\[df2.Subject==i\\]\\[0:5\\], axis=1).reset_index(drop=True)",
                           "df3 = pd.concat([df3, df2.loc[df2.Subject==i][0:5]], axis=1).reset_index(drop=True)")]}
               }
IL6_FIXES = {IL6PRED: {'il6': [("clf = load_model(model_path)", ""),
                       ("model_path = os.path.join(script_directory, \'..\', \'Models\', \'RF_model\')", ""),
                       ("Prediction Module start from here =====================",
                        """Prediction Module start from here =====================\n"""
                        """    script_directory = os.path.dirname(os.path.abspath(__file__))\n"""
                        """    model_path = os.path.join(script_directory, \'..\', \'Models\', \'RF_model\')\n"""
                        """    clf = load_model(model_path)""")]},
               }


# Only web software
IL4PRED,  IL10PRED = 'IL4pred', 'IL10pred'


VAXIGNML_DIC =     {'name': 'vaxign-ML', 'version': DEFAULT_VERSION, 'home': 'VAXIGNML_HOME'}

bepiPattern = 'immuno'

READ_URL = 'https://github.com/scipion-chem/scipion-chem-IIITD'


SEL_PARAM_MAP = {'abcWindow': 'window', 'abcThres': 'Threshold', 'abcFilter': 'filter',
                 'lbModel': "for", 'lbLength': 'expect', 'lbThres': 'val'}
SEL_PARAM_VALUE_MAP = {'LBtope_Fixed': 'fix', 'LBtope_Fixed_non_redundant': 'fixnr',
                       'LBtope_Variable': 'flx', 'LBtope_Variable_non_redundant': "flxnr", 'LBtop_Confirm': 'flx2'}

EVAL_PARAM_MAP = {
  TOXINPRED: {'method': {'SVM (Swiss-Prot)': 1, 'SVM (Swiss-Prot) + Motif': 2, 'SVM (TrEMBL)': 3, 'SVM (TrEMBL) + Motif': 4,
                'QM Monopeptide(Swiss-Prot)': 5, 'QM Monopeptide (TrEMBL)': 6, 'QM Dipeptide(Swiss-Prot)': 7,
                'QM Dipeptide (TrEMBL)': 8}
  },
  IL4PRED: {'method': {'SVM': 0, 'Merci motif': 1, 'Hybrid (SVM + motif)': 2, 'SwissProt': 3}
  },
  IL10PRED: {'method': {'SVM': 0, 'Random Forest': 1}
  },
  ALGPRED2: {'terminus': {'AAC based RF': 0, 'Hybrid (RF+BLAST+MERCI)': 4}
  },
  TOXINPRED2: {'terminus': {'AAC based RF': 0, 'Hybrid (RF+BLAST+MERCI)': 4}
  }
}

SEQ_LIMITS = {TOXINPRED: 50, IL4PRED: 50, IFNEPITOPE: 30}

TOXIN2WARN = '''ToxinPred2 is developed for predicting toxicity of proteins. In case user is interested in predicting 
toxicity of peptides then users should use ToxinPred3, which is specifically designed for peptides'''

ONLINE_WARN = '''This program does not run in local but through a website. This can significantly increase the protocol 
computation time.'''

SELSUM = '''1) "ABCpred-1": {'software': 'ABCpred', 'abcWindow': '16', 'abcThres': 0.51, 'abcFilter': 'on'}
2) "LBtope-1": {'software': 'LBtope', 'lbModel': 'LBtope_Variable', 'lbThres': '60', 'lbLength': 15}
'''

EVALSUM = '''1) "ToxinPred3-1": {'software': 'ToxinPred3', 'toxinMethod': 'Machine Learning (ML)', 'toxinThval': 0.38}
2) "IFNepitope2-1": {'software': 'IFNepitope2', 'ifnHost': 'Human', 'ifnThval': 0.49, 'ifnWindow': 8}
3) "AlgPred2-1": {'software': 'AlgPred2', 'algMethod': 'Hybrid (RF+BLAST+MERCI)', 'algThval': 0.3}
4) "IL4pred-1": {'software': 'IL4pred', 'il4Method': 'Hybrid', 'il4Thval': 0.2}
5) "IL5pred-1": {'software': 'IL5pred', 'il5Thval': 0.21, 'il5Window': 9}
6) "IL6pred-1": {'software': 'IL6pred', 'il6Thval': 0.11, 'il6Window': 10}
7) "IL10pred-1": {'software': 'IL10pred', 'il10Method': 'SVM', 'il10Thval': -0.3}
8) "IL13pred-1": {'software': 'IL13pred', 'il13Thval': 0.06, 'il13Window': 9}
9) "ToxinPred2-1": {'software': 'ToxinPred2', 'toxin2Method': 'Hybrid (RF+BLAST+MERCI)', 'toxin2Thval': 0.6}
'''

# ============================================================================
# Conformational/structure-based and sequence-based B-cell epitope predictors
# folded in here from their own standalone plugins (scipion-chem-scannet,
# scipion-chem-discotope, scipion-chem-tmbed, scipion-chem-signalp) and, in
# EpiDope's case, from scipion-chem core: small, easily-installable programs
# belong in scipion-chem-immuno rather than in their own dedicated plugin,
# and a predictor that ships its own binary does not belong in the core.
# ============================================================================

# EpiDope (Collatz et al. 2021) predicts LINEAR (sequence-based) B-cell
# epitopes from per-residue scores. Installed automatically by cloning the
# upstream repo and creating its conda environment exactly from the
# epidope.yml it ships: it pins an old, fragile stack (Python 3.6,
# TensorFlow 1.13, ELMo/AllenNLP), so the environment is built directly at
# EPIDOPE_HOME via '-p' to give downstream code a stable
# '<EPIDOPE_HOME>/bin/epidope' path regardless of where conda keeps its
# named environments.
EPIDOPE_DIC = {
    'name': 'epidope',
    'version': '1.0',
    'home': 'EPIDOPE_HOME',
}

EPIDOPE_UPSTREAM_URL = 'https://github.com/rnajena/EpiDope'

# Score cut-off carried over from this project's own validated EpiDope runs.
EPIDOPE_DEFAULT_THRESHOLD = 0.818

EPIDOPE_NOINSTALL_WARNING = (
    'Installation could not be completed because the local EpiDope '
    'installation has not been found. Run \'scipion3 installb epidope\' to '
    'install it automatically. Please check the plugin README file for more '
    f'details: {EPIDOPE_UPSTREAM_URL}'
)

# Raw output columns, one row per residue, in the per-accession CSV EpiDope
# writes to '<outdir>/epidope/<accession>.csv' (confirmed reading EpiDope's
# own epidope2.py::output_results).
EPIDOPE_RAW_COLUMNS = ('position', 'aminoacid', 'score')
EPIDOPE_ACCESSION_COLUMN = 'Accession'
EPIDOPE_RESIDUE_COLUMN = 'Residue'
EPIDOPE_SCORE_COLUMN = 'EpiDope score'

# ScanNet (Tubiana, Schneidman-Duhovny & Wolfson 2022, Apache-2.0) predicts
# CONFORMATIONAL (structure-based) B-cell epitopes. Installed automatically
# by cloning the upstream repo and building a dedicated conda env with its
# legacy dependency stack (Python 3.6.12, old TensorFlow/Keras/numba, taken
# as-is from ScanNet's own requirements.txt). No Docker runtime.
SCANNET_DIC = {
    'name': 'ScanNet',
    'version': DEFAULT_VERSION,
    'home': 'SCANNET_HOME',
    'activation': 'SCANNET_ACTIVATION_CMD',
}

IMMUNO_README_URL = 'https://github.com/scipion-chem/scipion-chem-immuno'
SCANNET_UPSTREAM_URL = 'https://github.com/jertubiana/ScanNet'

SCANNET_NOINSTALL_WARNING = (
    'Installation could not be completed because the local ScanNet '
    "installation has not been found or its conda environment could not be "
    "activated. Run 'scipion3 installb ScanNet' to install it automatically. "
    f"Please check this plugin's README file for more details: {IMMUNO_README_URL}"
)

# Prediction mode is ALWAYS '--mode epitope --noMSA', not configurable: the
# MSA mode requires a local sequence-database (UniRef30) + HH-blits
# installation, a heavy dependency this project deliberately avoids.
SCANNET_PREDICTION_MODE = 'epitope'

# Raw output columns (confirmed reading ScanNet's own
# predict_bindingsites.py::write_predictions): 'Model,Chain,Residue Index,
# Sequence,Binding site probability'.
SCANNET_RAW_CHAIN_COLUMN = 'Chain'
SCANNET_RAW_RESIDUE_COLUMN = 'Sequence'
SCANNET_RAW_SCORE_COLUMN = 'Binding site probability'

# DiscoTope-3.0 (DTU Health Tech, CC BY-NC 4.0 free academic use): installed
# automatically via git+pip (no academic-request form). Setup: clone, pip
# install, unzip the bundled XGBoost ensemble weights (models.zip), and
# pre-warm the ESM-IF1 weight cache.
DISCOTOPE_DIC = {
    'name': 'DiscoTope',
    'version': '3.0',
    'home': 'DISCOTOPE_HOME',
    'activation': 'DISCOTOPE_ACTIVATION_CMD',
}

DISCOTOPE_UPSTREAM_URL = 'https://github.com/Magnushhoie/DiscoTope-3.0'

DISCOTOPE_NOINSTALL_WARNING = (
    'Installation could not be completed because the local DiscoTope-3.0 '
    "installation has not been found or its conda environment could not be "
    "activated. Run 'scipion3 installb DiscoTope' to install it "
    f"automatically. Please check this plugin's README file for more "
    f'details: {IMMUNO_README_URL}'
)

# 'calibrated_score', not the raw 'DiscoTope-3.0_score': the authors publish
# reference thresholds for this column (Frontiers in Immunology 2024, Hoie
# et al.), not for the raw one.
DISCOTOPE_RAW_RESIDUE_COLUMN = 'residue'
DISCOTOPE_RAW_SCORE_COLUMN = 'calibrated_score'
DISCOTOPE_DEFAULT_THRESHOLD = 0.90

# TMbed (Bernhofer & Rost 2022, Apache-2.0): installed automatically in its
# own conda env. Its ProtT5 encoder weights (public, unauthenticated
# HuggingFace download) are pre-warmed into a plugin-local directory during
# installation.
TMBED_DIC = {
    'name': 'TMbed',
    'version': '1.0.2',
    'home': 'TMBED_HOME',
    'activation': 'TMBED_ACTIVATION_CMD',
}

# Files expected inside the TMbed model dir after T5Encoder's own
# save_pretrained() call. Accepts EITHER serialization format (older
# pytorch_model.bin/spiece.model, or newer safetensors/fast-tokenizer json)
# since save_pretrained() re-serializes using whatever the installed
# transformers version defaults to.
TMBED_T5_MODEL_REQUIRED_FILES = ('config.json',)
TMBED_T5_MODEL_WEIGHT_FILE_ALTERNATIVES = ('model.safetensors', 'pytorch_model.bin')
TMBED_T5_MODEL_TOKENIZER_FILE_ALTERNATIVES = ('tokenizer.json', 'spiece.model')

TMBED_DOWNLOAD_URL = 'https://github.com/BernhoferM/TMbed'

TMBED_NOINSTALL_WARNING = (
    'Installation could not be completed because the local TMbed conda '
    "environment and/or its cached ProtT5 encoder weights have not been "
    "found. Run 'scipion3 installb TMbed' to install it automatically. "
    f"Please check this plugin's README file for more details: {IMMUNO_README_URL}"
)

# TMbed '--out-format 1': merges strand/helix confidence tiers into a single
# upper-case letter per class (B/H), reports the signal peptide as 'S', and
# reports non-membrane residues as 'i'/'o' instead of format 0's '.'.
TMBED_OUT_FORMAT = '1'

# SignalP-6.0 is academic-use only software (DTU Health Tech), not
# redistributable: same class of constraint as BepiPred/NetMHCpan/
# NetMHCIIpan -- never installed automatically. The user downloads it
# manually (institutional email required) and points to it via
# scipion.conf.
SIGNALP_DIC = {
    'name': 'SignalP',
    'version': '6.0',
    'python_bin': 'SIGNALP_PYTHON_BIN',
    'binary_name': 'SIGNALP_BINARY_NAME',
    'model_dir': 'SIGNALP_MODEL_DIR',
}

SIGNALP_DEFAULT_BINARY_NAME = 'signalp6'
SIGNALP_DEFAULT_ORGANISM = 'other'

SIGNALP_DOWNLOAD_URL = 'https://services.healthtech.dtu.dk/services/SignalP-6.0/'

SIGNALP_NOINSTALL_WARNING = (
    'Installation could not be completed because the local SignalP-6.0 '
    'installation has not been found. Due to academic license restrictions, '
    f'DTU Health Tech does not allow redistributing this package: download it '
    f'manually from {SIGNALP_DOWNLOAD_URL} (requires an academic account), build a '
    'dedicated venv (Python 3.10, torch>1.7,<2, numpy<2) and set '
    'SIGNALP_PYTHON_BIN/SIGNALP_MODEL_DIR in scipion.conf. Please check '
    f"this plugin's README file for more details: {IMMUNO_README_URL}"
)

# Multi-epitope construct assembly, folded in from scipion-chem-epitope-construct.
# pwchem core keeps its own genetic-algorithm-based multi-epitope protocol
# (ProtOptimizeMultiEpitope); this one uses a different, non-GA chimeric
# assembly method and stays a separate protocol here in immuno rather than
# merging the two. Wraps NO external tool: pure selection/assembly logic.

# Linker sequences (standard-of-field convention for multi-epitope vaccine
# design, multiple independent sources agree -- not a fixed biological
# rule):
CONSTRUCT_LINKER_BCELL = 'KK'          # intra-B-cell: preserves individual epitope specificity
CONSTRUCT_LINKER_HTL = 'GPGPG'         # intra-HTL: universal spacer (Livingston et al. 2002)
CONSTRUCT_LINKER_CTL = 'AAY'           # intra-CTL: mammalian proteasomal cleavage site
CONSTRUCT_LINKER_INTERBLOCK = 'GPGPG'  # between distinct-class blocks: same universal spacer
CONSTRUCT_LINKER_ADJUVANT = 'EAAAK'   # optional adjuvant linker (Arai et al. 2001, rigid)

# Default number of epitopes kept per class -- necessary in practice: a
# real GP120 run produced 18 valid HTL candidates alone, too many for a
# manageable construct.
CONSTRUCT_DEFAULT_TOP_N_PER_CLASS = 3

# ============================================================================
# Whole-sequence annotators folded in from their own standalone plugins
# (scipion-chem-iapred, scipion-chem-netcleave, scipion-chem-stackglyembed).
# All three annotate, they never filter: the decision of what to discard
# belongs to the workflow, not to an individual evaluation protocol.
# ============================================================================

# IApred (Miles et al. 2025) scores whole-protein intrinsic antigenicity
# with an SVM over aggregated physicochemical features. It replaces
# VaxiJen, the historical reference for this task, which is not
# open-source and offers neither a downloadable standalone binary nor a
# documented public API. Installed automatically by cloning the upstream
# repository, whose pretrained SVM model comes along with the clone, into
# a dedicated conda environment.
IAPRED_DIC = {
    'name': 'IApred',
    'version': DEFAULT_VERSION,
    'home': 'IAPRED_HOME',
    'activation': 'IAPRED_ACTIVATION_CMD',
    'script_name': 'IAPRED_SCRIPT_NAME',
}

IAPRED_UPSTREAM_URL = 'https://github.com/sebamiles/IApred'
IAPRED_DEFAULT_SCRIPT_NAME = 'IApred.py'

IAPRED_NOINSTALL_WARNING = (
    'Installation could not be completed because the local IApred '
    'installation has not been found or its conda environment could not be '
    "activated. Run 'scipion3 installb IApred' to install it automatically. "
    'Please check the plugin README file for more details: '
    f'{IAPRED_UPSTREAM_URL}'
)

# IApred.py has a real, undocumented limitation, confirmed by reading the
# upstream script rather than assumed: for sequences shorter than 20 aa it
# writes the literal text 'Sequence too short' in the score column instead
# of a number, and 'N/A' in the category column. An assembled
# multi-epitope construct is normally well over 20 aa, but a single short
# B-cell candidate can hit it.
IAPRED_MIN_SEQUENCE_LENGTH = 20
IAPRED_TOO_SHORT_CATEGORY = f'Not evaluated (sequence < {IAPRED_MIN_SEQUENCE_LENGTH} aa)'

# NetCleave (Amengual-Rigo et al. 2021) predicts C-terminal proteasomal
# cleavage, the antigen-processing step that decides whether an MHC-I
# binder is ever actually generated. Installed automatically by cloning
# the upstream repository, whose bundled pretrained model
# (data/models/I_mass-spectrometry_HLA/) comes along with the clone, into
# a dedicated conda environment for its tensorflow/keras/scikit-learn
# stack.
#
# Licensing caveat: the upstream repository carries NO LICENSE file at
# all, confirmed against the repository itself rather than assumed. It is
# installed anyway, as is common for academic research code published
# alongside a paper, but this is a judgment call and not a
# license-cleared certainty.
NETCLEAVE_DIC = {
    'name': 'NetCleave',
    'version': DEFAULT_VERSION,
    'home': 'NETCLEAVE_HOME',
    'activation': 'NETCLEAVE_ACTIVATION_CMD',
}

NETCLEAVE_UPSTREAM_URL = 'https://github.com/pepamengual/NetCleave'

NETCLEAVE_NOINSTALL_WARNING = (
    'Installation could not be completed because the local NetCleave '
    'installation has not been found or its conda environment could not be '
    "activated. Run 'scipion3 installb NetCleave' to install it "
    'automatically. Please check the plugin README file for more details: '
    f'{NETCLEAVE_UPSTREAM_URL}'
)

# The bundled pretrained model (mass spectrometry, MHC-I, generic HLA
# family) ships inside the cloned repository. NetCleave.py resolves it as
# 'data/models/{mhc_class}_{technique}_{mhc_family}' RELATIVE TO THE
# PROCESS CWD, not to its own directory, so these three values must match
# that folder name exactly and the caller must run it with
# cwd=getNetCleaveDir(). The model is only ever scored against via
# --score_fasta, never retrained at runtime.
NETCLEAVE_MHC_CLASS = 'I'
NETCLEAVE_TECHNIQUE = 'mass_spectrometry'
NETCLEAVE_MHC_FAMILY = 'HLA'

# CORRECTED: the previous 'MIT license' claim here was WRONG,
# never actually verified -- the real upstream repo (nafcoder/StackGlyEmbed)
# has NO LICENSE file at all (confirmed via 'gh api': the repo's 'license'
# field is null, no LICENSE-named file in its root). No explicit open-source
# license means no formal redistribution/reuse permission is granted by
# default under copyright law, though this is a common (if technically
# ambiguous) situation for academic research code shared on GitHub
# alongside a paper. Auto-cloning it is a judgment call, not a
# license-cleared certainty like the other five auto-installed tools in
# this project (all confirmed real OSS licenses) -- flagged for the user
# to decide, not resolved unilaterally here. The heavy/heterogeneous
# dependency stack (torch + tensorflow + transformers + protein-bert)
# itself is still fine for Scipion's conda installer to handle; disk usage
# alone is not a reason to require a manual install.
# defineBinaries clones the upstream repo (for the pretrained classifier
# pickles under prediction/) and pre-warms all three embedding models'
# weight caches (ProteinBERT, ESM-2, ProtT5) so the protocol runs fully
# offline afterwards (predict_local.py itself already sets
# HF_HUB_OFFLINE=1/TRANSFORMERS_OFFLINE=1).
STACKGLYEMBED_DIC = {
    'name': 'StackGlyEmbed',
    'version': DEFAULT_VERSION,
    'home': 'STACKGLYEMBED_HOME',
    'activation': 'STACKGLYEMBED_ACTIVATION_CMD',
    'esm_model_name': 'STACKGLYEMBED_ESM_MODEL_NAME',
}

STACKGLYEMBED_DEFAULT_ESM_MODEL_NAME = 'facebook/esm2_t33_650M_UR50D'
# Same encoder as the original repo's 'Rostlab/prot_t5_xl_uniref50', just
# fp16/no decoder (see predict_local.py's module docstring, point 2).
STACKGLYEMBED_DEFAULT_T5_MODEL_NAME = 'Rostlab/prot_t5_xl_half_uniref50-enc'

# CORRECTED: the real repo is nafcoder/StackGlyEmbed (verified
# via 'gh api' -- 'GaryChan-lab/StackGlyEmbed' returns 404, does not
# exist). 'nafcoder' matches the paper's first author (Nafi).
STACKGLYEMBED_UPSTREAM_URL = 'https://github.com/nafcoder/StackGlyEmbed'

STACKGLYEMBED_NOINSTALL_WARNING = (
    'Installation could not be completed because the local StackGlyEmbed '
    "installation has not been found or its conda environment could not be "
    "activated. Run 'scipion3 installb StackGlyEmbed' to install it "
    'automatically. Please check the plugin README file for more details: '
    f'{STACKGLYEMBED_UPSTREAM_URL}'
)

# N-glycosylation sequon: Asn-Xaa-[Ser/Thr], Xaa != Pro. Lookahead-based
# (not a plain finditer) so overlapping sequons are all found (e.g.
# 'NPNSTPNST' reports N at position 1 AND 7). Returns 1-indexed Asn
# positions relative to whatever sequence string it is run against.
STACKGLYEMBED_SEQUON_PATTERN = r'(?=(N[^P][ST]))'
