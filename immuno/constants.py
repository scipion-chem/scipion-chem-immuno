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
IIITD_DIC = {'name': 'IIITD',    'version': '3.0',
             'home': 'IIITD_HOME', 'activation': 'IIITD_ACTIVATION_CMD',
             'browser': 'IIITD_BROWSER', 'browserPath': 'IIITD_BROWSER_PATH'}

VAXIGNML_DIC =     {'name': 'vaxign-ML', 'version': DEFAULT_VERSION, 'home': 'VAXIGNML_HOME'}

bepiPattern = 'immuno'

READ_URL = 'https://github.com/scipion-chem/scipion-chem-IIITD'


SEL_PARAM_MAP = {'abcWindow': 'window', 'abcThres': 'Threshold', 'abcFilter': 'filter'}

EVAL_PARAM_MAP = {
  'ToxinPred': {'method': {'SVM (Swiss-Prot)': 1, 'SVM (Swiss-Prot) + Motif': 2, 'SVM (TrEMBL)': 3, 'SVM (TrEMBL) + Motif': 4,
                'QM Monopeptide(Swiss-Prot)': 5, 'QM Monopeptide (TrEMBL)': 6, 'QM Dipeptide(Swiss-Prot)': 7,
                'QM Dipeptide (TrEMBL)': 8}
  },
  'IL4pred': {'method': {'SVM': 0, 'Merci motif': 1, 'Hybrid (SVM + motif)': 2, 'SwissProt': 3}
  },
  'IL10pred': {'method': {'SVM': 0, 'Random Forest': 1}
  },
  'AlgPred2': {'terminus': {'AAC based RF': 0, 'Hybrid (RF+BLAST+MERCI)': 4}
  },
  'ToxinPred2': {'terminus': {'AAC based RF': 0, 'Hybrid (RF+BLAST+MERCI)': 4}
  }
}

TOXIN2WARN = '''ToxinPred2 is developed for predicting toxicity of proteins. In case user is interested in predicting 
toxicity of peptides then users should use our old server ToxinPred, which is specifically designed for peptides'''

SELSUM = '''1) "ABCpred-1": {'software': 'ABCpred', 'abcWindow': '16', 'abcThres': 0.51, 'abcFilter': 'on'}
2) "LBtope-1": {'software': 'LBtope', 'lbModel': 'LBtope_Variable', 'lbThres': '60', 'lbLength': 15}
'''

EVALSUM = '''1) "ToxinPred-1": {'software': 'ToxinPred', 'toxinMethod': 'SVM', 'toxinSVMMethod': 'SVM(Swiss-Prot)', 'toxinQMMethod': 'Monopeptide(Swiss-Prot)', 'toxinEval': 10.0, 'toxinThval': 0.0}
2) "AlgPred2-1": {'software': 'AlgPred2', 'algMethod': 'AAC based RF', 'algThres': 0.3}
3) "IL4pred-1": {'software': 'IL4pred', 'il4Method': 'Hybrid', 'il4Thval': 0.2}
4) "IL10pred-1": {'software': 'IL10pred', 'il10Method': 'SVM', 'il10Thval': -0.3}
5) "IFNepitope-1": {'software': 'IFNepitope', 'ifnMethod': 'Hybrid'}
'''