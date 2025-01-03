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

IIITD_PACKAGES = ["toxinpred3", "toxinpred2", "algpred2", "ifnepitope2", "il5pred", "il13pred", "clbtope"]

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
               TOXINPRED2: {'toxinpred2': [("np.loadtxt(file_name, delimiter=',')", "np.loadtxt(file_name, delimiter=',', ndmin=2)")]},
               IFNEPITOPE: {'ifnepitope2': [('blastdb1 = nf_path + "\/\.\.\/blast_db\/human_db"', ''),
                            ('blastdb2 = nf_path + "\/\.\.\/blast_db\/mouse_db"', 'blastdb = nf_path + "\/\.\.\/blast_db"')]},
               ALGPRED2: {'algpred2': [("np.loadtxt(file_name, delimiter=',')", "np.loadtxt(file_name, delimiter=',', ndmin=2)"),
                          ("df3 = df3.concat(df2.loc[df2.Subject==i][0:5], axis=1).reset_index(drop=True)",
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

EVALSUM = f'''1) "ToxinPred-1": {{'software': '{TOXINPRED}', 'toxinMethod': 'SVM', 'toxinSVMMethod': 'SVM(Swiss-Prot)', 'toxinQMMethod': 'Monopeptide(Swiss-Prot)', 'toxinEval': 10.0, 'toxinThval': 0.0}}
2) "AlgPred2-1": {{'software': '{ALGPRED2}', 'algMethod': 'AAC based RF', 'algThval': 0.3}}
3) "IL4pred-1": {{'software': '{IL4PRED}', 'il4Method': 'Hybrid', 'il4Thval': 0.2}}
4) "IL10pred-1": {{'software': '{IL10PRED}', 'il10Method': 'SVM', 'il10Thval': -0.3}}
5) "IFNepitope-1": {{'software': '{IFNEPITOPE}', 'ifnMethod': 'Hybrid'}}
'''
