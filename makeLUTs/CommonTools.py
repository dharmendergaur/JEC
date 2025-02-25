from collections import OrderedDict as OD

import numpy as np
import pandas as pd

map_CaloIEta_to_CaloTool_mpEta = OD([
    (1, 1),
    (2, 2),
    (3, 3),
    (4, 4),
    (5, 5),
    (6, 6),
    (7, 7),
    (8, 8),
    (9, 9),
    (10, 10),
    
    (11, 11),
    (12, 12),
    (13, 13),
    (14, 14),
    (15, 15),
    (16, 16),
    (17, 17),
    (18, 18),
    (19, 19),
    (20, 20),
    
    (21, 21),
    (22, 22),
    (23, 23),
    (24, 24),
    (25, 25),
    (26, 26),
    (27, 27),
    (28, 28),
    (29, 29),
    (30, 29),
    
    (31, 30),
    (32, 31),
    (33, 32),
    (34, 33),
    (35, 34),
    (36, 35),
    (37, 36),
    (38, 37),
    (39, 38),
    (40, 39),

    (41, 40),
])


# log(Et) -----------------------------------------------
def convert_Et_to_logEt(Et):
    return np.log(Et)

def convert_logEt_to_Et(logEt):
    return np.exp(logEt)


# GenEt/L1Et --------------------------------------------
def convert_GenEt_to_GenEtByL1Et(GenEt, L1Et):
    return GenEt / L1Et

def convert_GenEtByL1Et_to_GenEt(GenEtByL1Et, L1Et):
    return GenEtByL1Et * L1Et


# log(GenEt/L1Et) ---------------------------------------
def convert_GenEt_to_logGenEtByL1Et(GenEt, L1Et):
    return np.log( GenEt / L1Et )

def convert_logGenEtByL1Et_to_GenEt(logGenEtByL1Et, L1Et):
    return L1Et * np.exp(logGenEtByL1Et)
#--------------------------------------------------------


def prepareDataframeForSFs(sL1JetTowerIEtaAbs, iEtaBinRange, sL1JetEt, PtRangeMin=10.0, PtRangeMax=255.0, snVtx='', nVtx=48):
    dict_iEta_Et = OD([ (sL1JetTowerIEtaAbs, []), (sL1JetEt, []) ])
    if snVtx:
        dict_iEta_Et[snVtx] = []
        
    for iEta in iEtaBinRange:
        list_pt      = np.arange(PtRangeMin, PtRangeMax+1.0)
        list_ietabin = [iEta] * len(list_pt)
        dict_iEta_Et[sL1JetTowerIEtaAbs].extend(list_ietabin) 
        dict_iEta_Et[sL1JetEt].extend(list_pt) 
        if snVtx:
            list_nVtx = [nVtx] * len(list_pt)
            dict_iEta_Et[snVtx].extend(list_nVtx) 
            
    data_SFs = pd.DataFrame(dict_iEta_Et)
    return data_SFs
#--------------------------------------------------------


def convert_CaloToolMPEta_to_IEta(CaloToolMPEta):
    IEta = None
    for IEta_tmp in map_CaloIEta_to_CaloTool_mpEta.keys():
        if map_CaloIEta_to_CaloTool_mpEta[ IEta_tmp ] == CaloToolMPEta:
            IEta = IEta_tmp

    return IEta



def GaussianFunction(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def calculate_errorOfRatio(N, D, eN, eD):
    #R = N / D
    return np.sqrt( ( ((1/D)**2) * (eN**2) ) +  ( ((N/(D**2))**2) * (eD**2) ) )
