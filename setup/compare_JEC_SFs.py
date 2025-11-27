# %%
import os
import sys
import argparse
from collections import OrderedDict
from collections import OrderedDict as OD
import csv
import json

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib
import hist
import mplhep as hep

from xgboost import XGBRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from hyperopt import fmin, tpe, hp, STATUS_OK, Trials, space_eval # hyperparameter optimization
import pickle


from CommonTools import (
    convert_Et_to_logEt, 
    convert_logEt_to_Et, 
    convert_GenEt_to_GenEtByL1Et, 
    convert_GenEtByL1Et_to_GenEt, 
    convert_GenEt_to_logGenEtByL1Et, 
    convert_logGenEtByL1Et_to_GenEt,
    convert_CaloToolMPEta_to_IEta,
    prepareDataframeForSFs,
    GaussianFunction,
    calculate_errorOfRatio
)


# JEC2024Ev0_DCSOnlyJSON_14_0_7_ZSHF3p5GeV
JEC_SFs_files_dict = {   
    "2023_MC": {
        "HBEF": {
            "PtAll": {
                "fileName" : '../data/L1T_Jet_SFs_tmp_vL1JetEt_PUS_PhiRing_PFJetEtCorr_MLTarget_GenEt_versionTmp_forPU33.csv',
            }
        }
    },

    "JEC2024Ev0_DCSOnlyJSON_14_0_7": {
        "HBEF": {
            "PtAll": {
                "fileName" : '../data/L1T_Jet_SFs_tmp_vL1JetEt_PUS_PhiRing_PFJetEtCorr_MLTarget_GenEt_versionTmp_forPU33.csv',
            }
        }
    },
    
    "JEC2024Ev0_DCSOnlyJSON_14_0_7_ZSHF3p5GeV": {
        "saveAs": '../data/L1T_Jet_SFs_JEC2024Ev0_DCSOnlyJSON_14_0_7_ZSHF3p5GeV_L1JetEt_PUS_PhiRing_HBE_logGenEtByL1Et_atPU33_HF_GenEtByL1Et_atPU33.csv',
        "HBE": {
            "PtAll": {
                "fileName" : '../data/BDTModel_vL1JetEt_PUS_PhiRing_PFJetEtCorr_MLTarget_logGenEtByL1Et_versionTmp_HBEF_PtAll.pkl',
                "MLTarget": 'logGenEtByL1Et',
                "PUForSFComputation": 33
            },
        },
        "HF": {
            "PtAll": {
                "fileName" : '../data/BDTModel_vL1JetEt_PUS_PhiRing_PFJetEtCorr_MLTarget_GenEtByL1Et_versionTmp_HBEF_PtAll.pkl',
                "MLTarget": 'GenEtByL1Et',
                "PUForSFComputation": 33
            },
        },        
    },
}
sOutDir = "plots_CompareJEC_SFs_JEC2024Ev0_DCSOnlyJSON_14_0_7_ZSHF3p5GeV"






printLevel = PrintLevel = 5
iEtaBins = [i for i in range(1, 42) if i!=29]
sL1JetEt_PUS_ChunkyDonut = 'L1JetEt_PUS_ChunkyDonut'
sL1JetEt_PUS_PhiRing     = 'L1JetEt_PUS_PhiRing'
sOfflineJetEt            = 'PFJetEtCorr'
sGenJetEt                = 'GenJetEt'
sL1JetTowerIEtaAbs       = 'L1JetTowerIEtaAbs'
L1JetPtThrsh             = 10.0 # GeV
L1JetPtMax               = 255.0 # GeV
RefJetPtThrsh            = 10.0 # GeV
snVtx                    = 'nVertexReco'


sSF                    = "SF"


NCompPtBins = 16 # 16 # No. of compressed pT bins
calibSF_L1JetPtRange = [15., 255., 1.] # [<lowest pT>,  <hightest pT>,  <pT bin width>] # pT range for SFs to read from Syed's SF.csv file
LUT_PtRange = [0., 255., 1.] # pT range for SFs for LUT
SF_forZeroPt = 1.0


IEta_Cat_forML = OD()
#IEta_Cat_forML['HB'] = [ 1, 16]
#IEta_Cat_forML['HE12a'] = [17, 26]
#IEta_Cat_forML['HE2b'] = [27, 28]
#IEta_Cat_forML['HF30to32'] = [30, 32]
#IEta_Cat_forML['HF33to36'] = [33, 36]
#IEta_Cat_forML['HF37to41'] = [37, 41]
#IEta_Cat['HBEF'] = [ 1, 41]
IEta_Cat_forML['HBEF'] = [ 1, 41]
IEta_Cat_forML['HBE']  = [ 1, 28]
IEta_Cat_forML['HF']   = [30, 41]

Pt_Cat_forML = OD()
#Pt_Cat_forML['Ptlt60'] = [ 0, 60]
#Pt_Cat_forML['Ptgt60'] = [60, L1JetPtMax]
#Pt_Cat_forML['Ptlt25']   = [ 0, 25]
#Pt_Cat_forML['Pt25to35'] = [25, 35]
#Pt_Cat_forML['Pt35to60'] = [35, 60]
#Pt_Cat_forML['Pt60to90'] = [60, 90]
#Pt_Cat_forML['Ptgt90']   = [90, L1JetPtMax]
Pt_Cat_forML['PtAll'] = [10, L1JetPtMax]



print(f"JEC_SFs_files_dict: {json.dumps(JEC_SFs_files_dict, indent=4) }")

# %%
JEC_SFs_dict = {}

for JECVersionName in JEC_SFs_files_dict:
    sL1JetEt                 = sL1JetEt_PUS_ChunkyDonut if 'Donut' in JECVersionName else sL1JetEt_PUS_PhiRing
    
    print(f"{JECVersionName = } ")
    print(f"{sL1JetEt = }")
    
    JEC_SFs_toUse = None
    
    for iEtaCat in JEC_SFs_files_dict[JECVersionName]:
        if iEtaCat not in IEta_Cat_forML: 
            print(f"iEtaCat {iEtaCat} not in IEta_Cat_forML {IEta_Cat_forML}")
            continue
            
        iEtaBinRange = list(range(IEta_Cat_forML[iEtaCat][0], IEta_Cat_forML[iEtaCat][1]+1))
        print(f"{iEtaCat = }, {iEtaBinRange = }")
        
        for PtCat in JEC_SFs_files_dict[JECVersionName][iEtaCat]:
            PtRangeMin   = Pt_Cat_forML[PtCat][0]
            PtRangeMax   = Pt_Cat_forML[PtCat][1]   
            fileName_i           = JEC_SFs_files_dict[JECVersionName][iEtaCat][PtCat]["fileName"]
            MLTarget_i           = JEC_SFs_files_dict[JECVersionName][iEtaCat][PtCat]["MLTarget"] if "MLTarget" in JEC_SFs_files_dict[JECVersionName][iEtaCat][PtCat] else ""
            PUForSFComputation_i = JEC_SFs_files_dict[JECVersionName][iEtaCat][PtCat]["PUForSFComputation"] if "PUForSFComputation" in JEC_SFs_files_dict[JECVersionName][iEtaCat][PtCat] else 0
            
            print(f"{PtCat = }: {PtRangeMin = }, {PtRangeMax = }")
            print(f"{fileName_i = }, {MLTarget_i = }, {PUForSFComputation_i = }")
            data_SFs_i = None
            
            if '.pkl' in fileName_i:
                xgb_i = pickle.load(open(fileName_i, "rb"))                  
                train_vars = xgb_i.get_booster().feature_names
                #print(f"xgb_i: {xgb_i}")                
                print(f"xgb_i.get_booster().feature_names: {xgb_i.get_booster().feature_names}")
                
                # prepareDataframeForSFs -------
                snVtx_toUse = snVtx if snVtx in train_vars else ''
                data_SFs_i = prepareDataframeForSFs(
                    sL1JetTowerIEtaAbs=sL1JetTowerIEtaAbs, iEtaBinRange=iEtaBinRange, 
                    sL1JetEt=sL1JetEt, PtRangeMin=PtRangeMin, PtRangeMax=PtRangeMax, 
                    snVtx=snVtx_toUse, nVtx=PUForSFComputation_i
                )                                
                
                sL1JetEt_forML  = None
                if   MLTarget_i == 'GenEt':
                    sL1JetEt_forML  = sL1JetEt    
                    data_SFs_i[sL1JetEt_forML] = data_SFs_i[sL1JetEt]

                elif MLTarget_i == 'logGenEt':    
                    sL1JetEt_forML  = 'log%s' % (sL1JetEt)
                    data_SFs_i[sL1JetEt_forML] = convert_Et_to_logEt( data_SFs_i[sL1JetEt] )

                elif MLTarget_i == 'GenEtByL1Et':    
                    sL1JetEt_forML  = sL1JetEt
                    data_SFs_i[sL1JetEt_forML] = data_SFs_i[sL1JetEt]

                elif MLTarget_i == 'logGenEtByL1Et':    
                    sL1JetEt_forML  = 'log%s' % (sL1JetEt)
                    data_SFs_i[sL1JetEt_forML] = convert_Et_to_logEt( data_SFs_i[sL1JetEt] )
                
                if sL1JetEt_forML not in train_vars:
                    print(f" {sL1JetEt_forML} not in train_vars {train_vars} \t\t **** ERROR **** ")
                    exit(0)
                    
                # xgb_i.predict() -----
                sL1JetEt_forML_predict = "%s_predict" % (sL1JetEt_forML)
                sL1JetEt_predict       = "%s_predict" % (sL1JetEt)

                data_SFs_i[sL1JetEt_forML_predict] = xgb_i.predict(data_SFs_i[train_vars])
                
                if   MLTarget_i == 'GenEt':
                    data_SFs_i[sL1JetEt_predict] = data_SFs_i[sL1JetEt_forML_predict]

                elif MLTarget_i == 'logGenEt':
                    data_SFs_i[sL1JetEt_predict] = convert_logEt_to_Et( data_SFs_i[sL1JetEt_forML_predict] )

                elif MLTarget_i == 'GenEtByL1Et':
                    data_SFs_i[sL1JetEt_predict] = convert_GenEtByL1Et_to_GenEt( data_SFs_i[sL1JetEt_forML_predict], data_SFs_i[sL1JetEt] )

                elif MLTarget_i == 'logGenEtByL1Et':
                    data_SFs_i[sL1JetEt_predict] = convert_logGenEtByL1Et_to_GenEt( data_SFs_i[sL1JetEt_forML_predict], data_SFs_i[sL1JetEt] )

                
                # SF --------------------
                data_SFs_i[sSF]                    = data_SFs_i[sL1JetEt_predict] / data_SFs_i[sL1JetEt]
                if printLevel >= 11:
                    print("iEtaBins_i: {}".format(iEtaBins_i))
                    print("data_SFs_i: {}".format(data_SFs_i.describe()))
                    
                    
            elif '.csv' in fileName_i:
                # prepareDataframeForSFs -------
                snVtx_toUse = ''
                data_SFs_i = prepareDataframeForSFs(
                    sL1JetTowerIEtaAbs=sL1JetTowerIEtaAbs, iEtaBinRange=iEtaBinRange, 
                    sL1JetEt=sL1JetEt, PtRangeMin=PtRangeMin, PtRangeMax=PtRangeMax, 
                    snVtx=snVtx_toUse, nVtx=PUForSFComputation_i
                )                                
                data_SFs_i[sSF] = 1.0                

                # read JEC SF ---------------------
                with open(fileName_i, mode='r') as fipFileCalibSF_toRead:
                    calibSF_csv_reader = csv.DictReader(fipFileCalibSF_toRead)
                    
                    for calibSF_csv_row in calibSF_csv_reader:
                        iEta_tmp = int(calibSF_csv_row[sL1JetTowerIEtaAbs]) # str
                        l1JetPt_tmp = float(calibSF_csv_row[ sL1JetEt ])
                        SF_tmp = float(calibSF_csv_row[ sSF ]) 
                        
                        data_SFs_i.loc[ (data_SFs_i[sL1JetTowerIEtaAbs] == iEta_tmp) & (abs(data_SFs_i[sL1JetEt] - l1JetPt_tmp) < 1e-3), sSF] = SF_tmp
                        
                #print(f"After {data_SFs_i.describe() = }")
                #print(f"{data_SFs_i.head() = }")
                
            if JEC_SFs_toUse is None:
                JEC_SFs_toUse = data_SFs_i
            else:
                JEC_SFs_toUse = pd.concat([JEC_SFs_toUse, data_SFs_i])    

                
    print(f"JEC_SFs_toUse: {JEC_SFs_toUse.describe()}"  )
    JEC_SFs_dict[JECVersionName] = JEC_SFs_toUse
    
    if "saveAs" in JEC_SFs_files_dict[JECVersionName]:
        sOpFileName_SFs = JEC_SFs_files_dict[JECVersionName]["saveAs"]
        JEC_SFs_toUse.to_csv(sOpFileName_SFs, index=False)
        print("Wrote {} \n".format(sOpFileName_SFs)) 
                
            
                
    
    



# %%
marker_color_list = ['r', 'b', 'darkviolet', 'c', 'orange', 'green']
marker_style_list = ["o", "o", "o",          "o", "o",      "o", "o", "X", '>', '^', 'v', "s", "+", 'x', '*']
marker_size_list  = [3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2]

sOutDir_toUse = sOutDir
os.makedirs(sOutDir_toUse, exist_ok=True)

for iEtaBin in iEtaBins:
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(5,4), layout='constrained')
    yMin_ =  9999
    yMax_ = -9999
    i_ = 0
    
    for JECVersionName in JEC_SFs_files_dict:
        sL1JetEt                 = sL1JetEt_PUS_ChunkyDonut if 'Donut' in JECVersionName else sL1JetEt_PUS_PhiRing
        #JEC_SFs_dict[JECVersionName]
        
        JEC_SFs_toUse_1 = JEC_SFs_dict[JECVersionName]
        JEC_SFs_toUse = JEC_SFs_toUse_1[(
            (JEC_SFs_toUse_1[sL1JetTowerIEtaAbs] == iEtaBin)
        )]
        
        print(f"{JEC_SFs_toUse[sL1JetEt] = }")
        print(f"{JEC_SFs_toUse[sSF] = }")
        
        axs.plot(JEC_SFs_toUse[sL1JetEt], JEC_SFs_toUse[sSF], label=JECVersionName, color='None', marker=marker_style_list[i_], markersize=marker_size_list[i_], markerfacecolor=marker_color_list[i_], markeredgecolor='None')        
        i_ += 1
        
        yMin1_ = np.min(JEC_SFs_toUse[sSF])
        yMax1_ = np.max(JEC_SFs_toUse[sSF])        
        if yMin1_ < yMin_: yMin_ = yMin1_
        if yMax1_ > yMax_: yMax_ = yMax1_

        
    axs.set_xlabel('%s [GeV]' % ('L1T pT'))
    axs.set_ylabel(r'SF = GEN pT/L1T pT')
    axs.set_title('iEta %d ' % (iEtaBin))
    if yMin_ < 0.0 or yMax_ > 4.0:
        axs.set_ylim(0., 4)
    else:
        axs.set_ylim(yMin_, yMax_)
    axs.set_xlim(15, 255)
    #axs.legend(bbox_to_anchor=(0.1, 1.05), loc='upper left', borderaxespad=0.9, ncol=2)
    axs.legend(bbox_to_anchor=(0.05, 1.05), loc='upper left', borderaxespad=0.9, ncol=1)
    axs.margins(y=0.3)
    axs.grid()
    fig.savefig('%s/JECSFs_vs_pT_iEta_%d.pdf' % (sOutDir_toUse, iEtaBin))   
    plt.close(fig)    
        
        

    

# %%



