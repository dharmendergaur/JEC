#!/bin/env python

'''
Update Pt-Eta bins for SFs to be compatible with firware requirements
'''


import os
import sys
import math 
from array import array
#import ROOT as R
#R.gROOT.SetBatch(False)  ## Don't print histograms to screen while processing

from collections import OrderedDict
import csv
import glob
import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
from operator import xor
import array as array
import statistics
import copy

PrintLevel       = 0 #2
JetShapes        = ['Default' ]
PUSAlgosAll      = ['RawPUS_phiDefault'] # ['RawPUS', 'RawPUS_phiDefault']

nLinesToRead = -1 # Max line to read from SF.csv file 

# Read L1Jet CalibLayer2 SFs from csv files provided by Syed
useAdditionalSFForLUT = False # 8/7 factor is needed for PhiRing as different PU estimation considered in CMSSW and Andrew's computation
icalibSF = 0 # 0, 1
calibSFLable = ['SF'][icalibSF]  
sipFileCalibSF = {
    'Default': {
       #'RawPUS': { # Chunky donut
       #    'fileName': '../data/L1T_Jet_SFs_2024v0_20240209_L1JetEt_PUS_ChunkyDonut_v0_HBE_logGenEtByL1Et_atPU33_HF_GenEtByL1Et_atPU33.csv', 
       #    'SFLabel': ['SF'][icalibSF],
       #    'L1JetPtVarName':'L1JetEt_PUS_ChunkyDonut',
       #    'additionalCorrForLUT': 1.0,
       #},
       
       'RawPUS_phiDefault': {
           'fileName': '', #'../data/L1T_Jet_SFs_2024v0_20240209_L1JetEt_PUS_PhiRing_HBE_logGenEtByL1Et_atPU33_HF_GenEtByL1Et_atPU33.csv',
           'SFLabel': ['SF'][icalibSF],
           'L1JetPtVarName':'L1JetEt_PUS_PhiRing',
           'additionalCorrForLUT': 1, # 8/7 factor is needed for PhiRing as different PU estimation considered in CMSSW and Andrew's computation
       },
    }, 


}

## 2024v0
#sipFileCalibSF['Default']['RawPUS_phiDefault']['fileName'] = '../data/L1T_Jet_SFs_2024v0_20240209_L1JetEt_PUS_PhiRing_HBE_logGenEtByL1Et_atPU33_HF_GenEtByL1Et_atPU33.csv'
## 2024v0p1_ZSHF3p5GeV
#sipFileCalibSF['Default']['RawPUS_phiDefault']['fileName'] = '../data/L1T_Jet_SFs_2024v0p1_13_3_0_ZSHF3p5GeV_20240219_L1JetEt_PUS_PhiRing_HBE_logGenEtByL1Et_atPU33_HF_GenEtByL1Et_atPU33.csv'
## 2024v0p2_13_3_0_L1SFvOlivier20240219
#sipFileCalibSF['Default']['RawPUS_phiDefault']['fileName'] = '../data/L1T_Jet_SFs_2024v0p2_13_3_0_L1SFvOlivier20240219_20240219_L1JetEt_PUS_PhiRing_HBE_logGenEtByL1Et_atPU33_HF_GenEtByL1Et_atPU33.csv'
## 2024v0p4_L1SFLLR20240311woZSHF
#sipFileCalibSF['Default']['RawPUS_phiDefault']['fileName'] = '../data/L1T_Jet_SFs_2024v0p4_13_3_0_L1SFLLR20240311woZSHF_L1JetEt_PUS_PhiRing_HBE_logGenEtByL1Et_atPU33_HF_GenEtByL1Et_atPU33.csv'
## 2024v0p5_13_3_0_L1SFLLR20240311wZSHF4p5GeV
#sipFileCalibSF['Default']['RawPUS_phiDefault']['fileName'] = '../data/L1T_Jet_SFs_2024v0p5_13_3_0_L1SFLLR20240311wZSHF4p5GeV_L1JetEt_PUS_PhiRing_HBE_logGenEtByL1Et_atPU33_HF_GenEtByL1Et_atPU33.csv'
## 2024v0p6_13_3_0_JECMultiplierMax4
#sipFileCalibSF['Default']['RawPUS_phiDefault']['fileName'] = '../data/L1T_Jet_SFs_2024v0_20240209_L1JetEt_PUS_PhiRing_HBE_logGenEtByL1Et_atPU33_HF_GenEtByL1Et_atPU33.csv'
## 2024v0p6_13_3_0_PUEtCapAt200GeV
#sipFileCalibSF['Default']['RawPUS_phiDefault']['fileName'] = '../data/L1T_Jet_SFs_2024v0p6_13_3_0_PUEtCapAt200GeV_L1JetEt_PUS_PhiRing_HBE_logGenEtByL1Et_atPU33_HF_GenEtByL1Et_atPU33.csv'
## JEC2024v0p7_13_3_0_L1SFLLR20240311wZSHF4p5GeV_w2023Ddata
#sipFileCalibSF['Default']['RawPUS_phiDefault']['fileName'] = '../data/L1T_Jet_SFs_JEC2024v0p7_13_3_0_L1SFLLR20240311wZSHF4p5GeV_w2023Ddata_L1JetEt_PUS_PhiRing_HBE_logGenEtByL1Et_atPU33_HF_GenEtByL1Et_atPU33.csv'
## JEC2024v0p0_13_3_0_w2023Ddata
#sipFileCalibSF['Default']['RawPUS_phiDefault']['fileName'] = '../data/L1T_Jet_SFs_JEC2024v0p0_13_3_0_w2023Ddata_L1JetEt_PUS_PhiRing_HBE_logGenEtByL1Et_atPU33_HF_GenEtByL1Et_atPU33.csv'
## JEC2024v0p1_13_3_0_ZSHF4p5GeV_w2023Ddata
#sipFileCalibSF['Default']['RawPUS_phiDefault']['fileName'] = '../data/L1T_Jet_SFs_JEC2024v0p1_13_3_0_ZSHF4p5GeV_w2023Ddata_L1JetEt_PUS_PhiRing_HBE_logGenEtByL1Et_atPU33_HF_GenEtByL1Et_atPU33.csv'
## 2024Ev0_DCSOnlyJSON_14_0_7
#sipFileCalibSF['Default']['RawPUS_phiDefault']['fileName'] = '../data/L1T_Jet_SFs_JEC2024Ev0_DCSOnlyJSON_14_0_7_L1JetEt_PUS_PhiRing_HBE_logGenEtByL1Et_atPU33_HF_GenEtByL1Et_atPU33.csv'
## 2024Ev0_DCSOnlyJSON_14_0_7_ZSHF3p5GeV
sipFileCalibSF['Default']['RawPUS_phiDefault']['fileName'] = '../data/L1T_Jet_SFs_JEC2024Ev0_DCSOnlyJSON_14_0_7_ZSHF3p5GeV_L1JetEt_PUS_PhiRing_HBE_logGenEtByL1Et_atPU33_HF_GenEtByL1Et_atPU33.csv'

sFilePtCompressedLUT_Ref = '' #'/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v7_20240209/JEC_2024_round1/CMSSW_13_3_0/src/L1Trigger/L1TCalorimeter/data/lut_pt_compress_2017v1.txt' # '' or 'lut_pt_compress_2017v1.txt file'

PtCompressedLUTVersion = '6Bits' # 'v2018' # 'v2018',  'v2022', '6Bits'
EtaCompressedLUT = False; #True, False;
EtaCompressedLUTVersion = '' # 'v2018' # 'v2018', '', 'v2022ChunkyDonut', 'v2022PhiRing', 'v2022Merged', ''


sLUTVersion = '2024Ev0_DCSOnlyJSON_14_0_7_ZSHF3p5GeV' # '2024_v0p1_13_3_0_ZSHF3p5GeV_HighGranularity' '2024_v0p6_13_3_0_PUEtCapAt200GeV_HighGranularity' '2024_v0p0_SFGt0' '2024_v0p0_SFGt0_HighGranularity' '2024_v0p0_SFGt1Lt2' '2024_v0p0_SFGt1Lt2p0_HighGranularity' '2024_v0p1_ZSHF3p5GeV_SFGt1Lt2_HighGranularity'
JECSF_boundary = [0.0, 9999.0] # [<lower bound>, <upper bound>] [0.0, 9999.0]
JECLUTSF_boundary = [0.0, 9999.0] # [<lower bound>, <upper bound>] [0.0, 9999.0] [1.0, 2.]
nBinsMaxForEtaCompressionLUT = 64 # no. of lines in eta compression LUT
makeLUTForIEta29 = [False, 1.0]
makeLUTForIEta41 = [False] # SFs for iEta=41 are missing in SFv6. Copy SFs from IEta=40
separatePtBinForLowPt = True # True: Use separate pT bin quantile for 0 < pT <= 15 GeV where SFs are truckated. NOT APPLICABLE when using firmware-fixed pT binning
IEtaBinOffsetForEtaCompressedLUT = 0 # 0: IEtaBin_forLUT = IEtaBin = [1, 41];  -1: IEtaBin_forLUT = IEtaBin - 1 = [0, 40]. # 0 give correct calibration.
sFOut_LUT_pt_compress        = 'lut_pt_compress_%s.txt' % (sLUTVersion)
sFOut_LUT_eta_compress       = 'lut_eta_compress_%s.txt' % (sLUTVersion)
sFOut_LUT_calib_compress_HEX = 'lut_calib_%s_ECALZS.txt' % (sLUTVersion)
sFOut_LUT_calib_compress     = 'lut_calib_%s_ECALZS_decimal.txt' % (sLUTVersion)

makeLUTsInUncompressedBins = False # True # make LUTs without compressing Pt and Eta bins
sFOut_LUT_pt_uncompress  = 'lut_pt_uncompress_%s.txt' % (sLUTVersion)
sFOut_LUT_eta_uncompress = 'lut_eta_uncompress_%s.txt' % (sLUTVersion)
sFOut_LUT_calib_uncompress = 'lut_calib_%s_ECALZS_decimal_uncompress.txt' % (sLUTVersion)

nBitsForPtComp = 6 if PtCompressedLUTVersion == '6Bits' else 4   # 4 # 4: for v2018 pT compression,  6: for finer granularity pT compression
NCompPtBins = int(2**nBitsForPtComp) # 16 # No. of compressed pT bins
calibSF_L1JetPtRange = [15., 255., 1.] # [<lowest pT>,  <hightest pT>,  <pT bin width>] # pT range for SFs to read from Syed's SF.csv file
LUT_PtRange = [0., 255., 1.] # pT range for SFs for LUT
SF_forZeroPt = 1.0


## Settings for JEC LUT --------
MakeJECLUT = True
nBitsCompBin = 12 if nBitsForPtComp==6 else 11 # 11: for 16 pT compression bins,    12: for 64 pT compression bins
nBitsPtComp  = nBitsForPtComp # pT compression bits:: 4: for 16 pT compression bins, 6: for 64 pT compression bins

UnitySFInBits = 512;
# 1./512 = 0.0019531250 = least count of SF
nBitsJECSF = 18
nBitsJECMultiplier = 10 # 10 (default): multiplier_max=2, 11: multiplier_max:4 ### <<<<<<<<<<<<<<<<<<< IMPORTANT setting <<<<<<<<<<<<<<<<<<
nBitsJECAdded      = nBitsJECSF - nBitsJECMultiplier # Default: 8
MaxAdded      = pow(2, nBitsJECAdded);      # 8: 256, 7: 128
MaxMultiplier = pow(2, nBitsJECMultiplier); # 10: 1024, 11: 2048
int8_t_Min = -128;
int8_t_Max =  127;

nCompBinMax = pow(2, nBitsCompBin);  # pow(2, 11) = 2048;  pow(2, 12) = 4096
## ------------------------------


if PtCompressedLUTVersion == 'v2018':
    separatePtBinForLowPt = False
    # '2018PtCompression has 1st bins from 1 <= pT <= 19, so set SF(pT <= 19) = SF(pT(19))
    #calibSF_L1JetPtRange = [15., 255., 1.] # [<lowest pT>,  <hightest pT>,  <pT bin width>] # pT range for SFs to read from Syed's SF.csv file
    calibSF_L1JetPtRange = [15., 255., 1.] # [<lowest pT>,  <hightest pT>,  <pT bin width>] # pT range for SFs to read from Syed's SF.csv file

if PtCompressedLUTVersion == '6Bits':
    separatePtBinForLowPt = False
    if nBitsForPtComp != 6:
        print(f"PtCompressedLUTVersion: {PtCompressedLUTVersion}, but nBitsForPtComp: {nBitsForPtComp}, != 6 \t\t **** ERROR ****")
        exit(0)

    
map_CaloIEta_to_CaloTool_mpEta = OrderedDict([
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

CaloTool_mpEtas = [1, 40]

CaloToolMPEtaBinsMerge_forEtaCompressedLUT_2018 = OrderedDict([
    ( 0, [*range(1, 5+1)]),  # 0
    ( 1, [*range(6, 9+1)]),  # 1
    ( 2, [*range(10,13+1)]), # 2
    ( 3, [*range(14,15+1)]), # 3
    ( 4, [*range(16,17+1)]), # 4
    ( 5, [*range(18,19+1)]), # 5
    ( 6, [*range(20,21+1)]), # 6
    ( 7, [22]), # 7
    ( 8, [23]), # 8
    ( 9, [24]), # 9
    (10, [25]), # 10
    (11, [26]), # 11
    (12, [*range(27,28+1)]), # 12
    (13, [*range(29,31+1)]), # 13
    (14, [*range(32,35+1)]), # 14
    (15, [*range(36,40+1)]), # 15
])
CaloToolMPEtaBinsMerge_forEtaUncompressedLUT = OrderedDict()
for CaloToolMPEta in range(1,41): # CaloToolMPEta=iEta for iEta=(1, 28),  CaloToolMPEta=iEta-1 for iEta=(30, 41)
    CaloToolMPEtaBinsMerge_forEtaUncompressedLUT[CaloToolMPEta - 1] = [ CaloToolMPEta ]
    

CaloToolMPEtaBinsMerge_forEtaCompressedLUT_2022_ChunkyDonut = OrderedDict([
    ( 0, [*range(1, 5+1)]),  # 0
    ( 1, [*range(6, 9+1)]),  # 1
    ( 2, [*range(10,13+1)]), # 2
    ( 3, [*range(14,15+1)]), # 3
    ( 4, [*range(16,17+1)]), # 4
    ( 5, [*range(18,19+1)]), # 5
    ( 6, [*range(20,21+1)]), # 6
    ( 7, [22]), # 7
    ( 8, [23]), # 8
    ( 9, [24]), # 9
    (10, [25]), # 10
    (11, [26]), # 11
    (12, [*range(27,28+1)]), # 12
    (13, [*range(29,31+1)]), # 13
    (14, [*range(32,35+1)]), # 14
    (15, [*range(36,40+1)]), # 15
])
IEtaBinsMerge_forPlots_ChunkyDonut = [
    #[1],
    [*range(1,6)],
    [*range(5,10)],
    [*range(9,14)],
    [*range(13,18)],
    [*range(17,22)],
    [*range(21,26)],
    [*range(25,31)],
    [*range(30,35)],
    [*range(34,41)],
]
IEtaBinsMerge_forPlots_1_ChunkyDonut = [
    #[1],
    [*range(1,8)],
    [*range(4,12)],
    [*range(8,16)],
    [*range(12,20)],
    [*range(16,24)],
    [*range(20,28)],
    [*range(24,32)],
    [*range(28,36)],
    [*range(32,38)],
    [*range(36,41)],
]
IEtaBinsMerge_forPlots_2_ChunkyDonut = [
    #[1],
    [*range(1,4+1)],
    [*range(5,11+1)],
    [12],
    [*range(13,15+1)],
    [*range(16,17+1)],
    [*range(18,21+1)],
    [*range(22,23+1)],
    [*range(24,25+1)],
    [*range(26,27+1)],
    [28],
    [30],
    [*range(31,32+1)],
    [*range(33,36+1)],
    [37],
    [*range(38,39+1)],
    [*range(40,41+1)],
]



CaloToolMPEtaBinsMerge_forEtaCompressedLUT_2022_PhiRing = OrderedDict([
    ( 0, [*range(1, 5+1)]),  # 0
    ( 1, [*range(6, 9+1)]),  # 1
    ( 2, [*range(10,13+1)]), # 2
    ( 3, [*range(14,15+1)]), # 3
    ( 4, [*range(16,17+1)]), # 4
    ( 5, [*range(18,19+1)]), # 5
    ( 6, [*range(20,21+1)]), # 6
    ( 7, [22]), # 7
    ( 8, [23]), # 8
    ( 9, [24]), # 9
    (10, [25]), # 10
    (11, [26]), # 11
    (12, [*range(27,28+1)]), # 12
    (13, [*range(29,31+1)]), # 13
    (14, [*range(32,35+1)]), # 14
    (15, [*range(36,40+1)]), # 15
])

IEtaBinsMerge_forPlots_PhiRing = [
    #[1],
    [*range(1,6)],
    [*range(5,10)],
    [*range(9,14)],
    [*range(13,18)],
    [*range(17,22)],
    [*range(21,26)],
    [*range(25,31)],
    [*range(30,35)],
    [*range(34,41)],
]
IEtaBinsMerge_forPlots_1_PhiRing = [
    #[1],
    [*range(1,8)],
    [*range(4,12)],
    [*range(8,16)],
    [*range(12,20)],
    [*range(16,24)],
    [*range(20,28)],
    [*range(24,32)],
    [*range(28,36)],
    [*range(32,38)],
    [*range(36,41)],
]
IEtaBinsMerge_forPlots_2_PhiRing = [
    #[1],
    [*range(1,4+1)],
    [*range(5,10+1)],
    [*range(11,12+1)],
    [*range(13,15+1)],
    [*range(16,17+1)],
    [*range(18,22+1)],
    [*range(23,26+1)],
    [27],
    [28],
    [30],
    [*range(31,32+1)],
    [33],
    [*range(34,36+1)],
    [*range(37,38+1)],
    [39],
    [*range(40,41+1)],
]


CaloToolMPEtaBinsMerge_forEtaCompressedLUT_2022_ChunkyDonutAndPhiRingMerged = OrderedDict([
    ( 0, [*range(1, 4+1)]),  # 0
    ( 1, [*range(5, 11+1)]),  # 1
    ( 2, [*range(12,13+1)]), # 2
    ( 3, [*range(14,15+1)]), # 3
    ( 4, [*range(16,17+1)]), # 4
    ( 5, [*range(18,21+1)]), # 5
    ( 6, [22]), # 6
    ( 7, [23]), # 7
    ( 8, [*range(24,25+1)]), # 8
    ( 9, [26]), # 9
    (10, [27]), # 10
    (11, [28]), # 11
    (12, [29]), # 12
    (13, [*range(30,31+1)]), # 13
    (14, [*range(32,38+1)]), # 14
    (15, [*range(39,40+1)]), # 15
])
IEtaBinsMerge_forPlots_2_2022_ChunkyDonutAndPhiRingMerged = [
    #[1],
    [*range(1,4+1)],
    [*range(5,11+1)],
    [*range(12,13+1)],
    [*range(14,15+1)],
    [*range(16,17+1)],
    [*range(18,21+1)],
    [22],
    [23],
    [*range(24,25+1)],
    [26],
    [27],
    [28],
    [30],
    [*range(31,32+1)],
    [*range(33,39+1)],
    [*range(40,41+1)],
]


CaloToolMPEtaBinsMerge_forEtaCompressedLUT = None
IEtaBinsMerge_forPlots   = None
IEtaBinsMerge_forPlots_1 = None
IEtaBinsMerge_forPlots_2 = None
if 'RawPUS' in PUSAlgosAll:
    CaloToolMPEtaBinsMerge_forEtaCompressedLUT = CaloToolMPEtaBinsMerge_forEtaCompressedLUT_2022_ChunkyDonut 
    IEtaBinsMerge_forPlots   = IEtaBinsMerge_forPlots_ChunkyDonut
    IEtaBinsMerge_forPlots_1 = IEtaBinsMerge_forPlots_1_ChunkyDonut
    IEtaBinsMerge_forPlots_2 = IEtaBinsMerge_forPlots_2_ChunkyDonut
elif 'RawPUS_phiDefault' in PUSAlgosAll:
    CaloToolMPEtaBinsMerge_forEtaCompressedLUT = CaloToolMPEtaBinsMerge_forEtaCompressedLUT_2022_PhiRing
    IEtaBinsMerge_forPlots   = IEtaBinsMerge_forPlots_PhiRing
    IEtaBinsMerge_forPlots_1 = IEtaBinsMerge_forPlots_1_PhiRing
    IEtaBinsMerge_forPlots_2 = IEtaBinsMerge_forPlots_2_PhiRing
IEtaBinsMerge_forPlots_2 = IEtaBinsMerge_forPlots_2_2022_ChunkyDonutAndPhiRingMerged    



pTBinsMerge = OrderedDict([ # https://github.com/cms-l1t-offline/L1Trigger-L1TCalorimeter/blob/master/lut_pt_compress_2017v1.txt
    (0, [0]),
    (1, [1, 19]),
    (2, [20]),
    (3, [21]),
    (4, [22, 25]),
    (5, [26, 29]),
    (6, [30, 34]),
    (7, [35, 40]),
    (8, [41, 48]),
    (9, [49, 58]),
    (10, [59, 71]),
    (11, [72, 88]),
    (12, [89, 113]),
    (13, [114, 155]),
    (14, [156, 254]),
    (15, [255 ]),
])


'''
Aggreed pT compression binning with 6-bits:
SF=1 for pT = 0  GeV (1 bin)
1 SF for 1 <= pT <= 15 GeV (1 bin)
1 GeV steps from 16 - 31 GeV (16 bins)
2 GeV steps from 32 - 61 GeV (15 bins)
3 GeV steps from 62 - 103 GeV (14 bins)
6 GeV steps from 104 - 193 GeV (15 bins)
1 SF for  194 <= pT <= 254 GeV (1 bin)
1 SF for  pT = 255 GeV (1 bin)
'''
pTCompressionBinEdges_6Bits = [
    # 64 bins, 65 bin edges,
    # [0, 1, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 65, 68, 71, 74, 77, 80, 83, 86, 89, 92, 95, 98, 101, 104, 110, 116, 122, 128, 134, 140, 146, 152, 158, 164, 170, 176, 182, 188, 194, 255, 256]    
    0,
    *range(  1,  15+1, 15),
    *range( 16,  31+1,  1),
    *range( 32,  61+1,  2),
    *range( 62, 103+1,  3),
    *range(104, 193+1,  6),
    *range(194, 254+1, 61),
    255,
    256
]




useAbsEtaBins = True
ETA_Bins = []
for iEta in range(-41,42):
    if iEta in [-29, 0, 29]:        continue;
    if useAbsEtaBins and iEta < 0:  continue;
    ETA_Bins.append(str(iEta))


colors_matplotlib = ['r', 'b', 'c', 'm', 'g', 'y', 'gray', 'olive', 'teal', 'orange']

    

def splitInQuantile(arrayOriginal, nQuantiles):
    '''
    Split numpy.array into sub-lists such that each sub-list have equal sum over it's elements.
    Return type:
        idx_quantiles: List of indices in the original array where quantiles are made. 0th and last index element of the original array are not included in it.

    '''
    
    arrayCumSum = np.cumsum(arrayOriginal)
    if PrintLevel >= 2:
        print("arrayOriginal ({}): {}, \nnQuantiles:{}, \narrayCumSum ({}): {}".format(len(arrayOriginal), arrayOriginal, nQuantiles, len(arrayCumSum), arrayCumSum))

    idx_quantiles = []
    for iQuant in range(1, nQuantiles):
        for idx in range(0, len(arrayCumSum)):
            if arrayCumSum[idx] >= arrayCumSum[-1]/nQuantiles*iQuant:
                idx_quantiles.append(idx)
                break
        
    arraySplitted = [ sublist.tolist()  for sublist in np.split(arrayOriginal, idx_quantiles)]
    split_sum = [sum(sublist) for sublist in arraySplitted]

    if PrintLevel >= 2:
        print("idx_quantiles ({}): {}, \narraySplitted ({}): {}, \nsplit_sum ({}): {}, \narrayOriginal ({}): {}".format(
            len(idx_quantiles), idx_quantiles,
            len(arraySplitted), arraySplitted,
            len(split_sum), split_sum,
            len(arrayOriginal), arrayOriginal))

    return idx_quantiles, arraySplitted
    
#idx_quantile, SFs_AvgOverEtaBins_splitInQuantile



def splitInQuantile_ConstDeltaSF(Pt_Bins_Read, SFsOriginal, nQuantiles):
    '''
    Spli numpy.array into sublists such that each sub-list has constant (SF_max - SF_min)
    Return type:
        idx_quantiles: List of indices in the original array where quantiles are made. 0th and last index element of the original array are not included in it.
    '''

    ## ignore arrayOriginal[0] = SF(pT=0) = 1
    SFs_truncated = SFsOriginal[1:]
    
    idx_pTBinQuantiles = [] 
    if separatePtBinForLowPt:
        # np.isclose(Pt_Bins_Read, calibSF_L1JetPtRange[0], atol=1.e-5).nonzero() : (array([15]),)
        idx_tmp = np.isclose(Pt_Bins_Read, calibSF_L1JetPtRange[0], atol=1.e-5).nonzero()[0][0] # idx in Pt_Bins where pT=15
        #idx_pTBinQuantiles.append( idx_tmp+1 )
        
        if PrintLevel >= 15:
            print("Pt_Bins_Read ({}): {}, \ncalibSF_L1JetPtRange[0]: {}, \t idx_tmp: {} \n idx_pTBinQuantiles ({}): {} \nnp.isclose(Pt_Bins_Read, calibSF_L1JetPtRange[0]): {} \n -->: {}".format(
                len(Pt_Bins_Read), Pt_Bins_Read, calibSF_L1JetPtRange[0], idx_tmp,
                len(idx_pTBinQuantiles), idx_pTBinQuantiles,
                np.isclose(Pt_Bins_Read, calibSF_L1JetPtRange[0], atol=1.e-5).nonzero(),
                np.flatnonzero( np.isclose(Pt_Bins_Read, calibSF_L1JetPtRange[0]) ) 
            ))
            

            
    nQuantiles_toUse = (nQuantiles - 1) if separatePtBinForLowPt else nQuantiles
    SFMin = np.amin(SFs_truncated)
    SFMax = np.amax(SFs_truncated)
    dSF   = (SFMax - SFMin) / nQuantiles_toUse
    SFQuants = [ SFMin + (dSF * iQuant) for iQuant in range(nQuantiles_toUse, -1, -1) ]

    if PrintLevel >= 2:
        print("\nsplitInQuantile_ConstDeltaSF():: \nSFsOriginal {}, \nSFs_truncated {}, \nSFMin {}, SFMax {}, nQuantiles {}, nQuantiles_toUse {}, SFQuants ({}) {}, \n idx_pTBinQuantiles ({}): {}".format(
            SFsOriginal, SFs_truncated,
            SFMin, SFMax, nQuantiles, nQuantiles_toUse,
            len(SFQuants), SFQuants,
            len(idx_pTBinQuantiles), idx_pTBinQuantiles
        ))

    iSFQuant = 0 #1 if separatePtBinForLowPt else 0
    for idx_SFsOriginal, SF in enumerate(SFsOriginal):
        # ignore SF(pT=0) = 1
        if idx_SFsOriginal == 0: continue

        if PrintLevel >= 2:
            print("\t idx_SFsOriginal {}, SF {}, iSFQuant {}, SFQuants[iSFQuant] {}, idx_pTBinQuantiles ({}): {}".format(
                idx_SFsOriginal, SF,
                iSFQuant, SFQuants[iSFQuant],
                len(idx_pTBinQuantiles), idx_pTBinQuantiles
            ))
        
        if SF < SFQuants[iSFQuant]:
            idx_pTBinQuantiles.append( idx_SFsOriginal )
            iSFQuant += 1


    SFsSplittedInPtBinQuantile = [ sublist.tolist()  for sublist in np.split(SFsOriginal, idx_pTBinQuantiles)]
    split_dSF = [np.amax(sublist) - np.amin(sublist)  for sublist in SFsSplittedInPtBinQuantile]
    
    if PrintLevel >= 2:
        print("splitInQuantile_ConstDeltaSF():: idx_pTBinQuantiles final ({}): {}".format(
            len(idx_pTBinQuantiles), idx_pTBinQuantiles
        ))
        print("SFsSplittedInPtBinQuantile: {}, \nsplit_dSF: {}".format(
            SFsSplittedInPtBinQuantile, split_dSF
        ))

    return idx_pTBinQuantiles, SFsSplittedInPtBinQuantile





def splitInQuantile_ConstDeltaSF_old(Pt_Bins_Read, SFsOriginal, nQuantiles):
    '''
    Spli numpy.array into sublists such that each sub-list has constant (SF_max - SF_min)
    Return type:
        idx_quantiles: List of indices in the original array where quantiles are made. 0th and last index element of the original array are not included in it.
    '''

    ## ignore arrayOriginal[0] = SF(pT=0) = 1
    SFs_truncated = SFsOriginal[1:]
    
    idx_pTBinEdges = [] # [0]
    if separatePtBinForLowPt:
        # np.isclose(Pt_Bins_Read, calibSF_L1JetPtRange[0], atol=1.e-5).nonzero() : (array([15]),)
        idx_tmp = np.isclose(Pt_Bins_Read, calibSF_L1JetPtRange[0], atol=1.e-5).nonzero()[0][0] # idx in Pt_Bins where pT=15
        #idx_pTBinEdges.append( idx_tmp+1 )
        
        if PrintLevel >= 15:
            print("Pt_Bins_Read ({}): {}, \ncalibSF_L1JetPtRange[0]: {}, \t idx_tmp: {} \n idx_pTBinEdges ({}): {} \nnp.isclose(Pt_Bins_Read, calibSF_L1JetPtRange[0]): {} \n -->: {}".format(
                len(Pt_Bins_Read), Pt_Bins_Read, calibSF_L1JetPtRange[0], idx_tmp,
                len(idx_pTBinEdges), idx_pTBinEdges,
                np.isclose(Pt_Bins_Read, calibSF_L1JetPtRange[0], atol=1.e-5).nonzero(),
                np.flatnonzero( np.isclose(Pt_Bins_Read, calibSF_L1JetPtRange[0]) ) 
            ))
            

            
    nQuantiles_toUse = (nQuantiles - 1) if separatePtBinForLowPt else nQuantiles
    SFMin = np.amin(SFs_truncated)
    SFMax = np.amax(SFs_truncated)
    dSF   = (SFMax - SFMin) / nQuantiles_toUse
    SFQuants = [ SFMin + (dSF * iQuant) for iQuant in range(nQuantiles_toUse, -1, -1) ]

    if PrintLevel >= 2:
        print("\nsplitInQuantile_ConstDeltaSF():: \nSFsOriginal {}, \nSFs_truncated {}, \nSFMin {}, SFMax {}, nQuantiles {}, nQuantiles_toUse {}, SFQuants ({}) {}, \n idx_pTBinEdges ({}): {}".format(
            SFsOriginal, SFs_truncated,
            SFMin, SFMax, nQuantiles, nQuantiles_toUse,
            len(SFQuants), SFQuants,
            len(idx_pTBinEdges), idx_pTBinEdges
        ))

    iSFQuant = 0 #1 if separatePtBinForLowPt else 0
    for idx_SFsOriginal, SF in enumerate(SFsOriginal):
        # ignore SF(pT=0) = 1
        if idx_SFsOriginal == 0: continue

        if PrintLevel >= 2:
            print("\t idx_SFsOriginal {}, SF {}, iSFQuant {}, SFQuants[iSFQuant] {}, idx_pTBinEdges ({}): {}".format(
                idx_SFsOriginal, SF,
                iSFQuant, SFQuants[iSFQuant],
                len(idx_pTBinEdges), idx_pTBinEdges
            ))
        
        if SF < SFQuants[iSFQuant]:
            idx_pTBinEdges.append( idx_SFsOriginal )
            iSFQuant += 1


    #idx_pTBinEdges.append( len(SFsOriginal) )
    if PrintLevel >= 2:
        print("idx_pTBinEdges final ({}): {}".format(
            len(idx_pTBinEdges), idx_pTBinEdges
        ))

    return



def readPtCompressedBinning(PtCompressionBinning_dict):
    idx_pTBinQuantiles = []

    Pt_list = (list(PtCompressionBinning_dict.keys()))
    for idx_Pt, Pt in enumerate(Pt_list):
        if idx_Pt == 0: continue
        
        PtCompressedBin = PtCompressionBinning_dict[Pt]
        if PrintLevel >= 20:
            print("readPtCompressedBinning():: Pt {}, PtCompressedBin {}".format(Pt, PtCompressedBin))

        Pt_lastBin = Pt_list[idx_Pt - 1]
        PtCompressedBin_lastPtBin = PtCompressionBinning_dict[Pt_lastBin]
        if PtCompressedBin != PtCompressedBin_lastPtBin:
            idx_pTBinQuantiles.append(idx_Pt)

    if PrintLevel >= 20:
        print("idx_pTBinQuantiles ({}): {}".format(len(idx_pTBinQuantiles), idx_pTBinQuantiles))
        
    return idx_pTBinQuantiles

        

                            

def convert_CaloToolMPEta_to_IEta(CaloToolMPEta):
    IEta = None
    for IEta_tmp in map_CaloIEta_to_CaloTool_mpEta.keys():
        if map_CaloIEta_to_CaloTool_mpEta[ IEta_tmp ] == CaloToolMPEta:
            IEta = IEta_tmp

    return IEta
        
    


def calculateSFInBits_addendZero(jet_hwPt, SFInDecimal):
    if (PrintLevel >= 20):
        print("calculateSFInBits_addendZero():: \n");
  
    jet_calibPt_target = int( ( float(jet_hwPt) * SFInDecimal ) + 0.5); # 0.5 to round off fractional part of the number
    if (PrintLevel >= 4):
        print("calculateSFInBits_addendZero():: {}, {}, PtCorr: {} / {} ".format(float(jet_hwPt), SFInDecimal, ( float(jet_hwPt) * SFInDecimal ), jet_calibPt_target))
  
  
    multiplier = int(SFInDecimal * float(UnitySFInBits));
    # multiplier is 10-bit, which is Power(2,10) = 1024 = MaxMultiplier
    if (multiplier >= MaxMultiplier ): # { # truncate multiplier to 1023
        multiplier = MaxMultiplier - 1;
    
    multiplier_0 = multiplier;

    if (PrintLevel >= 2): # {
        print(f"{multiplier_0 = }")
    

    if (jet_hwPt == 0 and abs(SFInDecimal - 1) < 0.000001):
        addPlusMult = multiplier_0; # multiplier=512 for SFInDecimal==1
        return addPlusMult;
  
  
    multiplier_wAddendZero = []
    #for (m_i = 0; m_i < MaxMultiplier; m_i++) {
    for m_i in range(MaxMultiplier):
        m_tmp = multiplier + m_i;
        if (m_tmp >= MaxMultiplier): break;

        a_i = calculateAddend(jet_hwPt, m_tmp, jet_calibPt_target);
        if (a_i == 0): multiplier_wAddendZero.append( m_tmp );

        if (a_i < 0): break;
    

    if (PrintLevel >= 2):
        print("multiplier_wAddendZero:: {}".format(multiplier_wAddendZero))
 
    #addPlusMult;    
    if ( len(multiplier_wAddendZero) > 0 ): # {
        idx_multiplier = int( len(multiplier_wAddendZero) / 2 );
        multiplier = multiplier_wAddendZero[idx_multiplier];
        addend = 0;
        
        jetPtCorr = ((jet_hwPt * multiplier) >> 9) + addend;
        addPlusMult = (addend << nBitsJECMultiplier) + multiplier;  # (addend << 10) + multiplier;
  
    else: 
  
        multiplier = multiplier_0;
  
        absDPt_Min = 99999;
        addend = 0;
        #for (int8_t addend_i = 0; addend_i < 127; addend_i++) { // int8_t range [-128, 127]
        for addend_i in range(int8_t_Max):
            jetPtCorr_i = ((jet_hwPt * multiplier) >> 9) + addend_i;
            dPt = int(jetPtCorr_i) - int(jet_calibPt_target);
            absDPt = abs(dPt);
            if (absDPt < absDPt_Min): # {
                addend = addend_i;
                absDPt_Min = absDPt;
      

            if (dPt == 0): # {
	            break;


        addPlusMult = (addend << nBitsJECMultiplier) + multiplier;  # (addend << 10) + multiplier;
  
    
    return addPlusMult;  
  

def calculateAddend(jet_hwPt, multiplier, jet_calibPt_target):
    addend = 0;
    jetPtCorr_i = ((jet_hwPt * multiplier) >> 9) + addend;
    dPt = int(jetPtCorr_i) - int(jet_calibPt_target);
    if (dPt == 0): # {
        return addend;
    elif (dPt < 0): # {
        absDPt_Min = 99999;
        addend = 0;
        for addend_i in range(int8_t_Max): # int8_t_Max: 127
            jetPtCorr_i = ((jet_hwPt * multiplier) >> 9) + addend_i;
            dPt = int(jetPtCorr_i) - int(jet_calibPt_target);
            absDPt = abs(dPt);
            if (absDPt < absDPt_Min):
                addend = addend_i;
                absDPt_Min = absDPt;
        return addend;

    elif (dPt > 0): # {
        absDPt_Min = 99999;
        addend = 0;
        for addend_i in range(-1, int8_t_Min, -1): # int8_t_Min: -128
            jetPtCorr_i = ((jet_hwPt * multiplier) >> 9) + addend_i;
            dPt = int(jetPtCorr_i) - int(jet_calibPt_target);
            absDPt = abs(dPt);
            if (absDPt < absDPt_Min): # {
                addend = addend_i;
                absDPt_Min = absDPt;             
        return addend;
    
    return addend;


def calibrateJet(jet_hwPt, addPlusMult, printLevel_i=0):
  if (printLevel_i >= 2):
    print("calibrateJet():: \n");
  
  # multiplier = addPlusMult & 0x3ff; #  0x3ff = 1023 = 11 1111 1111
  # addend = (addPlusMult >> 10);
  multiplier = addPlusMult & 0x7ff if nBitsJECMultiplier == 11 else addPlusMult & 0x3ff # hex(2**10-1): 0x3ff,  hex(2**11-1): 0x7ff
  addend = (addPlusMult >> nBitsJECMultiplier);
  jetPtCorr = ((jet_hwPt * multiplier) >> 9) + addend;

  return jetPtCorr;

if __name__ == '__main__':


    print(f"{sLUTVersion = }, {PUSAlgosAll = }, \n{EtaCompressedLUT = }, {EtaCompressedLUTVersion = }, \n{nBitsForPtComp = }, {NCompPtBins = }, {PtCompressedLUTVersion = }, {LUT_PtRange = }, {calibSF_L1JetPtRange = }, \n\n")

    calibSF_L1JetPtRangeMin = calibSF_L1JetPtRange[0]
    calibSF_L1JetPtRangeMax = calibSF_L1JetPtRange[1]
    calibSF_L1JetPtBinWidth = calibSF_L1JetPtRange[2]

    mplColors_list = ['orange', 'cyan', 'y', 'royalblue', 'darkviolet', 'red', 'olive', 'teal', ]



    ## Set PtCompression Binning ---------------------------------------
    PtCompressionBinning_dict = OrderedDict()
    if PtCompressedLUTVersion == 'v2018':
        print("PtCompressedLUT read from {}".format(sFilePtCompressedLUT_Ref))
        with open(sFilePtCompressedLUT_Ref, mode='r') as fPtCompressedLUT_Ref:
            lines = fPtCompressedLUT_Ref.readlines()

            for line0 in lines:
                line = line0.replace('\n', '')
                columns = line.split(' ')
                if '#' in columns[0]: continue

                PtCompressionBinning_dict[int(columns[0])] = [int(columns[1])]
                if PrintLevel >= 11:
                    print("line: {}!!! \t split: {}".format(line, line.split(' ')))

    if PtCompressedLUTVersion == '6Bits':
        print(f"PtCompressedLUTVersion: {PtCompressedLUTVersion}")
        for PtCompBin_ in range(len(pTCompressionBinEdges_6Bits) - 1):
            PtMin_ = pTCompressionBinEdges_6Bits[PtCompBin_    ]
            PtMax_ = pTCompressionBinEdges_6Bits[PtCompBin_ + 1] - 1
            for Pt_ in range(PtMin_, PtMax_+1, 1):
                PtCompressionBinning_dict[int(Pt_)] = PtCompBin_
    
    if PrintLevel >= 2:
        print("PtCompressionBinning_dict ({}): \n{}".format(len(PtCompressionBinning_dict), json.dumps(PtCompressionBinning_dict, indent=4)))
    #-------------------------------------------------------------------------

    
    
    print(list(sipFileCalibSF.keys()))
    for jetShape in JetShapes:
        if jetShape not in list(sipFileCalibSF.keys()):
            print("%s not in sipFileCalibSF " % (jetShape))
            continue

        print(sipFileCalibSF[jetShape].keys())
        for PUSAlgo in PUSAlgosAll:

            if PUSAlgo not in list(sipFileCalibSF[jetShape].keys()):
                print("%s not in sipFileCalibSF " % (PUSAlgo))
                continue


            sipFileCalibSF_toRead = sipFileCalibSF[jetShape][PUSAlgo]['fileName']
            print("sipFileCalibSF[{}][{}]: {} ".format(jetShape, PUSAlgo, sipFileCalibSF_toRead))
            print("ETA_Bins ({}): {} ".format(len(ETA_Bins), ETA_Bins))
            

            data_calibSFs =  OrderedDict()
            for Eta_bin in ETA_Bins:
                data_calibSFs[int(Eta_bin)] = OrderedDict()
                for Pt in np.arange(LUT_PtRange[0], LUT_PtRange[1]+1, LUT_PtRange[2]):
                    data_calibSFs[int(Eta_bin)][Pt] = 0.0
                    

            with open(sipFileCalibSF_toRead, mode='r') as fipFileCalibSF_toRead:
                calibSF_csv_reader = csv.DictReader(fipFileCalibSF_toRead)

                nLines = 0
                
                #print("\n\ncalibSF_csv_reader: {} {}\n\n\n".format(type(calibSF_csv_reader), calibSF_csv_reader))
                for calibSF_csv_row in calibSF_csv_reader:
                    iEta_tmp = int(calibSF_csv_row['L1JetTowerIEtaAbs']) # str
                    #l1JetPt_tmp = float(calibSF_csv_row['PhiRingEnergy'])
                    l1JetPt_tmp = float(calibSF_csv_row[ sipFileCalibSF[jetShape][PUSAlgo]['L1JetPtVarName'] ])
                    SF_tmp = float(calibSF_csv_row[ sipFileCalibSF[jetShape][PUSAlgo]['SFLabel'] ]) * sipFileCalibSF[jetShape][PUSAlgo]['additionalCorrForLUT']
                    if SF_tmp < JECSF_boundary[0]: SF_tmp = JECSF_boundary[0]
                    if SF_tmp > JECSF_boundary[1]: SF_tmp = JECSF_boundary[1]
                    
                    if PrintLevel >= 11:
                        print("calibSF_csv_row: {}, iEta {} {}, l1JetPt {} {}, SF {} {}={} {} * additional factor {}".format(
                            calibSF_csv_row, type(iEta_tmp),iEta_tmp, type(l1JetPt_tmp),l1JetPt_tmp, type(SF_tmp),SF_tmp, sipFileCalibSF[jetShape][PUSAlgo]['SFLabel'],float(calibSF_csv_row[ sipFileCalibSF[jetShape][PUSAlgo]['SFLabel'] ]), sipFileCalibSF[jetShape][PUSAlgo]['additionalCorrForLUT']))

                        
                    #if l1JetPt_tmp < calibSF_L1JetPtRangeMin or l1JetPt_tmp > calibSF_L1JetPtRangeMax: continue
                    

                    data_calibSFs[iEta_tmp][l1JetPt_tmp] = SF_tmp

                    if nLinesToRead > 0 and nLines > nLinesToRead: break
                    nLines += 1
                    


            
            # Updates SFs for lowpT and hightPt bins --------------------------------------------------------------
            for Eta_bin in list(data_calibSFs.keys()):
                if len(list(data_calibSFs[Eta_bin].keys())) == 0: # Not all Eta_Bins were used to extract SFs
                    continue

                #print("data_calibSFs[{}]: {},  {}".format(Eta_bin, data_calibSFs[Eta_bin].keys(), calibSF_L1JetPtRange[0]))
                #if float(calibSF_L1JetPtRange[0]) not in data_calibSFs[Eta_bin].keys():
                if calibSF_L1JetPtRange[0] not in data_calibSFs[Eta_bin].keys():
                    print("L1 pT {} not in data_calibSFs[ieta={}] \t\t **** ERROR ****".format(calibSF_L1JetPtRange[0], Eta_bin))
                sTmp = "Eta %d: \nLowPt: " % (Eta_bin)

                for Pt in np.arange(LUT_PtRange[0], calibSF_L1JetPtRange[0]):
                    if abs(Pt - 0) < 0.0001:
                        data_calibSFs[Eta_bin][Pt] = SF_forZeroPt
                        # Add Pt=0 at the beining in OrderedDict
                        #data_calibSFs[Eta_bin].update( {Pt : SF_forZeroPt } )
                        #data_calibSFs[Eta_bin].move_to_end(Pt, last=False)
                    else:
                        data_calibSFs[Eta_bin][Pt] = data_calibSFs[Eta_bin][calibSF_L1JetPtRange[0]]
                    sTmp += " ({}, {}),".format(Pt, data_calibSFs[Eta_bin][Pt])
                sTmp += "\nHighPt: "
                for Pt in np.arange(calibSF_L1JetPtRange[1]+1., LUT_PtRange[1]+1.):
                    data_calibSFs[Eta_bin][Pt] = data_calibSFs[Eta_bin][calibSF_L1JetPtRange[1]]
                    sTmp += " ({}, {}),".format(Pt, data_calibSFs[Eta_bin][Pt])
                if PrintLevel >= 12:
                    print(sTmp)
                if PrintLevel >= 12:
                    print("data_calibSFs[{}] after: {}".format(Eta_bin, data_calibSFs[Eta_bin]))
            # ------------------------------------------------------------------------------------------------------        


            # LUT for iEta=29 --------------------------------------------------------------------------------------
            if makeLUTForIEta29[0]:
                #print("LUT_PtRange ({}): {}".format(len(list(np.arange(LUT_PtRange[0], LUT_PtRange[1]+1, LUT_PtRange[2]))), list(np.arange(LUT_PtRange[0], LUT_PtRange[1]+1, LUT_PtRange[2]))))
                data_calibSFs[29] = OrderedDict()
                for Pt in list(np.arange(LUT_PtRange[0], LUT_PtRange[1]+1, LUT_PtRange[2])):
                    data_calibSFs[29][Pt] = makeLUTForIEta29[1]

                if PrintLevel >= 12:
                    print("data_calibSFs[29]: ({}) {}".format(len(list(data_calibSFs[29].keys())), data_calibSFs[29]))
            # ------------------------------------------------------------------------------------------------------
            
            
            # LUT for iEta=41 --------------------------------------------------------------------------------------
            if makeLUTForIEta41[0] and len(list(data_calibSFs[41].keys())) == 0:
                IEta_current   = 41
                IEta_reference = 40 # SFs from IEta_reference will be copied for IEta_current
                for Pt in list(np.arange(LUT_PtRange[0], LUT_PtRange[1]+1, LUT_PtRange[2])):
                    data_calibSFs[IEta_current][Pt] = data_calibSFs[IEta_reference][Pt]

                if PrintLevel >= 12:
                    print("data_calibSFs[{}]: ({}) {}".format(IEta_reference, len(list(data_calibSFs[IEta_reference].keys())), data_calibSFs[IEta_reference]))
                    print("data_calibSFs[{}]: ({}) {}".format(IEta_current, len(list(data_calibSFs[IEta_current].keys())), data_calibSFs[IEta_current]))
            # ------------------------------------------------------------------------------------------------------


            if PrintLevel >= 1:
                print("\n\ndata_calibSFs:: {}".format('*'*10))
                for Eta_bin in list(data_calibSFs.keys()):
                    print("\t iEta {}: {}".format(Eta_bin, data_calibSFs[Eta_bin]))


                    

            EtaBins_sorted = list(data_calibSFs.keys())
            EtaBins_sorted.sort()
            PtBins_forLUT  = list(np.arange(LUT_PtRange[0], LUT_PtRange[1]+1, LUT_PtRange[2]))
            nPtBins_forLUT = len(PtBins_forLUT)
            
            # write LUTs with uncompressed Pt and eta -------------------------------------------------------------
            if makeLUTsInUncompressedBins:
                print("\nmakeLUTsInUncompressedBins:: \nEtaBins_sorted ({}): {} \nPtBins_forLUT ({}): {}".format(
                    len(EtaBins_sorted), EtaBins_sorted,
                    nPtBins_forLUT, PtBins_forLUT
                ))
                
                sDir_LUT = "LUTs_%s/%s_%s_%s" % (sLUTVersion, jetShape, PUSAlgo, calibSFLable)
                sDir_LUT = sDir_LUT.replace("(","_")
                sDir_LUT = sDir_LUT.replace(")","_")
                sDir_LUT = sDir_LUT.replace("__","_")
                if not os.path.exists(sDir_LUT):
                    os.makedirs(sDir_LUT)

                fOut_LUT_pt_uncompress    = open("%s/%s" % (sDir_LUT, sFOut_LUT_pt_uncompress),    'w')
                fOut_LUT_eta_uncompress   = open("%s/%s" % (sDir_LUT, sFOut_LUT_eta_uncompress),   'w')
                fOut_LUT_calib_uncompress = open("%s/%s" % (sDir_LUT, sFOut_LUT_calib_uncompress), 'w')

                # write headers
                fOut_LUT_eta_uncompress.write('# MP ieta compression LUT\n')
                fOut_LUT_eta_uncompress.write('# Converts abs(MP ieta) (6 bits) into 6 bit index\n')
                fOut_LUT_eta_uncompress.write('# This is NOT calo ieta\n')
                fOut_LUT_eta_uncompress.write('# anything after # is ignored with the exception of the header\n')
                fOut_LUT_eta_uncompress.write('# the header is first valid line starting with #<header> versionStr nrBitsAddress nrBitsData </header>\n')
                fOut_LUT_eta_uncompress.write('#<header> v1 6 6 </header>\n')

                fOut_LUT_pt_uncompress.write('# PT compression LUT\n')
                fOut_LUT_pt_uncompress.write('# maps 8 bits to 8 bits\n')
                fOut_LUT_pt_uncompress.write('# the 1st column is the integer value after selecting bits 1:8\n')
                fOut_LUT_pt_uncompress.write('# anything after # is ignored with the exception of the header\n')
                fOut_LUT_pt_uncompress.write('# the header is first valid line starting with #<header> versionStr nrBitsAddress nrBitsData </header>\n')
                fOut_LUT_pt_uncompress.write('#<header> v1 8 8 </header>\n')
                
                
                IEtaBin_forLUT_col0 = None
                
                for IEta in EtaBins_sorted:
                    CaloTool_mpEta      = map_CaloIEta_to_CaloTool_mpEta[IEta]
                    IEtaBin_forLUT_col0 = CaloTool_mpEta
                    IEtaBin_forLUT_col1 = CaloTool_mpEta - 1
                    IEtaBin_forLUT      = IEtaBin_forLUT_col1
                    
                    if IEta == EtaBins_sorted[0]:
                        fOut_LUT_eta_uncompress.write('0 0  # dummy ieta_bin\n')                        
                    fOut_LUT_eta_uncompress.write('%d %d  # ieta %d\n' % (int(IEtaBin_forLUT_col0), int(IEtaBin_forLUT_col1), IEta))
                    
                    
                    for Pt in PtBins_forLUT:
                        PtBin_forLUT = Pt

                        
                        if IEta == EtaBins_sorted[0]:
                            fOut_LUT_pt_uncompress.write('%d %d\n' % (int(Pt), int(PtBin_forLUT)))
                            
                        sComments = ""
                        if Pt == PtBins_forLUT[0]:
                            sComments = "   # ieta %d, pt %d" % (IEta, PtBin_forLUT)

                        bin_eta_pt_forLUT = (nPtBins_forLUT * IEtaBin_forLUT) + PtBin_forLUT
                        #fOut_LUT_calib_uncompress.write('%d %f%s\n' % (int(bin_eta_pt_forLUT), data_calibSFs[IEta][Pt], sComments ))
                        #fOut_LUT_calib_uncompress.write('%d %d %f%s\n' % (int(bin_eta_pt_forLUT), int(Pt), data_calibSFs[IEta][Pt], sComments ))
                        fOut_LUT_calib_uncompress.write('%d %d %d %f%s\n' % (int(IEtaBin_forLUT), int(PtBin_forLUT), int(Pt), data_calibSFs[IEta][Pt], sComments ))

                        if PrintLevel >= 2:
                            print("(%d, %d):  (%d, %d) -> %d: \t %f" % (
                                IEta, Pt,
                                IEtaBin_forLUT, PtBin_forLUT, bin_eta_pt_forLUT,
                                data_calibSFs[IEta][Pt]                                
                            ))
            

                IEtaBin_forLUT_col0 += 1 # the next bin = 42
                while IEtaBin_forLUT_col0 < nBinsMaxForEtaCompressionLUT:
                    IEtaBin_forLUT_col1 = IEtaBin_forLUT_col0 - 1 
                    #fOut_LUT_eta_compress.write('%d %d\n' % (int(IEtaBin_forLUT_col0), int(IEtaBin_forLUT_col1)))
                    fOut_LUT_eta_uncompress.write('%d %d\n' % (int(IEtaBin_forLUT_col0), int(0)))
                    IEtaBin_forLUT_col0 += 1


                            
                fOut_LUT_pt_uncompress.close()
                fOut_LUT_eta_uncompress.close()
                fOut_LUT_calib_uncompress.close()

                print("Wrote %s/%s, %s/%s, %s/%s" % (
                    sDir_LUT, sFOut_LUT_pt_uncompress,
                    sDir_LUT, sFOut_LUT_eta_uncompress,
                    sDir_LUT, sFOut_LUT_calib_uncompress
                ))
                


            
            # list of Eta_Bins to consider. Not all Eta_Bins were used to extract SFs
            Eta_Bins_Read  = [EtaBin  for EtaBin in list(data_calibSFs.keys())  if len(list(data_calibSFs[EtaBin].keys())) > 0 ]
            nEta_Bins_Read = len(Eta_Bins_Read)
            Pt_Bins_Read   = list(data_calibSFs[ Eta_Bins_Read[0] ].keys())
            nPt_Bins_Read  = len(Pt_Bins_Read)
            if PrintLevel >= 0:
                #print("data_calibSFs.keys() ({}): {}".format(len(list(data_calibSFs.keys())), data_calibSFs.keys()))
                print("Eta_Bins_Read ({}): {}".format(nEta_Bins_Read, Eta_Bins_Read))
                print("Pt_Bins_Read ({}): {}".format(nPt_Bins_Read, Pt_Bins_Read))


            

            ### Store SFs into (nEtaBins x nPtBins) matrix-array ---------------------------------------------------------
            #SFs_array = np.zeros(shape=(int(ETA_Bins[-1])+1, len(list(data_calibSFs[int(ETA_Bins[0])].values())) ))
            SFs_array = np.zeros(shape=(nEta_Bins_Read, len(list(data_calibSFs[Eta_Bins_Read[0]].values())) ))
            if PrintLevel >= 5:
                print("\n\nSF_all initially ({}): {}".format(SFs_array.shape,SFs_array))

            #for Eta_bin0 in ETA_Bins:
            for Eta_bin0 in Eta_Bins_Read:
                Eta_bin = int(Eta_bin0)
                if Eta_bin not in list(data_calibSFs.keys()):
                    print("Eta {} not in data_calibSFs ".format(Eta_bin))
                    continue

                if PrintLevel >= 5:
                    print("\nEta_bin {} \ndata_calibSFs[Eta_bin]: {} ".format(Eta_bin, data_calibSFs[Eta_bin]))
                    print("data_calibSFs[{}] ({}): {} ".format(Eta_bin, len(list(data_calibSFs[Eta_bin].keys())), data_calibSFs[Eta_bin].keys()))
                if len(list(data_calibSFs[Eta_bin].keys())) == 0:
                    print("Eta_bin {}: len(list(data_calibSFs[Eta_bin].keys())) == 0".format(Eta_bin))
                    continue

                if PrintLevel >= 5:
                    print("Eta_Bins_Read.index(Eta_bin)): {} ".format(Eta_Bins_Read.index(Eta_bin)))
                    
                idx_Eta_Bin_Read = Eta_Bins_Read.index(Eta_bin)

                #SFs_inEtaBin = np.array( list(data_calibSFs[Eta_bin].values()) )
                SFs_inEtaBin_list = [data_calibSFs[Eta_bin][Pt] for Pt in np.arange(LUT_PtRange[0], LUT_PtRange[1]+1, LUT_PtRange[2]) ]
                SFs_inEtaBin = np.array( SFs_inEtaBin_list )
                #SFs_array[Eta_bin] = SFs_inEtaBin
                SFs_array[idx_Eta_Bin_Read] = SFs_inEtaBin

                if PrintLevel >= 5:
                    print("\nEta_bin {}: \nSF ({}): {} \nSF_all ({}): {}".format(
                        Eta_bin,
                        SFs_inEtaBin.shape, SFs_inEtaBin,
                        SFs_array.shape,SFs_array
                    ))

            #print("\n\nSF_all final ({}): \n{}".format(SFs_array.shape, np.array2string(SFs_array, max_line_width=2000, threshold=3000, precision=6) ))
            #print("\n\nSF_all final ({}): \n{}".format(SFs_array.shape, np.array_str(SFs_array, max_line_width=2000, precision=6) ))
            #print("\n\nSF_all final ({}): \n{}".format(SFs_array.shape, SFs_array.tolist() ))
            print("\n\nSF_all final ({}): ".format(SFs_array.shape ))
            for iEta in range(SFs_array.shape[0]):
                sTmp = [ (Pt, SFs_array[iEta][Pt]) for Pt in range(SFs_array.shape[1])]
                print("\tSFs_array[{}]: {}".format(iEta, sTmp))
            # ----------------------------------------------------------------------------------------------------------
            
            ## average SF over all eta bins
            SFs_AvgOverEtaBins = np.average(SFs_array, axis=0)
            print("\nSFs_AvgOverEtaBins ({}): {}".format(len(SFs_AvgOverEtaBins), SFs_AvgOverEtaBins))


                
            if PtCompressedLUTVersion in ['v2018', '6Bits']: # ----------------------
                # use 2018 PtCompression binning
                idx_PtQuantiles                    = readPtCompressedBinning(PtCompressionBinning_dict)
                SFs_AvgOverEtaBins_splitInQuantile = [ sublist.tolist()  for sublist in np.split(SFs_AvgOverEtaBins, idx_PtQuantiles) ]

            else:
                ### calculate pT compression bins ---------------------

                ### Divide pT in quantiles with constant DeltaSF over the quatiles  ---------------------------
                idx_PtQuantiles, SFs_AvgOverEtaBins_splitInQuantile = splitInQuantile_ConstDeltaSF(Pt_Bins_Read, SFs_AvgOverEtaBins, nQuantiles=NCompPtBins)


            Pt_quantiles = []
            for idx in idx_PtQuantiles:
                Pt_quantiles.append( (Pt_Bins_Read[idx-1] + Pt_Bins_Read[idx]) / 2. )

            if PrintLevel >= 0:
                print("idx_PtQuantiles ({}): {}, \nPt_quantiles ({}): {} ".format(
                    len(idx_PtQuantiles), idx_PtQuantiles,
                    len(Pt_quantiles), Pt_quantiles ))
            if PrintLevel >= 2:
                print("SFs_AvgOverEtaBins_splitInQuantile({}): {}".format(
                    len(SFs_AvgOverEtaBins_splitInQuantile), SFs_AvgOverEtaBins_splitInQuantile ))


            SFs_AvgInPtQuant_AvgOverEtaBins = SFs_AvgOverEtaBins.copy()
            Pt_AvgInPtQuant = []
            SFs_AvgInPtQuant_AvgOverEtaBins_Coarse = []
            for iQuant in range(NCompPtBins):
                idxStart                  = 0             if iQuant == 0               else idx_PtQuantiles[(iQuant-1)]
                idxEnd                    = nPt_Bins_Read if iQuant == (NCompPtBins-1) else idx_PtQuantiles[(iQuant)]
                Pt_Bins_Read_iQuant       = Pt_Bins_Read[idxStart : idxEnd]
                SFs_AvgOverEtaBins_iQuant = SFs_AvgOverEtaBins[idxStart : idxEnd]
                avgPt_iQuant              = statistics.mean(Pt_Bins_Read_iQuant)
                avgSF_iQuant              = statistics.mean(SFs_AvgOverEtaBins_iQuant)
                if separatePtBinForLowPt and iQuant == 0:
                    avgSF_iQuant          = SFs_AvgOverEtaBins_iQuant[1]
                Pt_AvgInPtQuant.append( avgPt_iQuant )
                SFs_AvgInPtQuant_AvgOverEtaBins_Coarse.append( avgSF_iQuant )
                if PrintLevel >= 4:
                     print("\n\niQuant {}, idxStart {}, idxEnd {}, \nPt_Bins_Read_iQuant ({}): {}, SFs_AvgOverEtaBins_iQuant ({}): {}, avgPt_iQuant {}, avgSF_iQuant {}, \nPt_AvgInPtQuant ({}): {}, \nSFs_AvgInPtQuant_AvgOverEtaBins_Coarse ({}): {}, \nSFs_AvgInPtQuant_AvgOverEtaBins before ({}): {}".format(
                         iQuant, idxStart, idxEnd,
                         len(Pt_Bins_Read_iQuant),Pt_Bins_Read_iQuant, len(SFs_AvgOverEtaBins_iQuant),SFs_AvgOverEtaBins_iQuant,
                         avgPt_iQuant, avgSF_iQuant,
                         len(Pt_AvgInPtQuant), Pt_AvgInPtQuant,
                         len(SFs_AvgInPtQuant_AvgOverEtaBins_Coarse), SFs_AvgInPtQuant_AvgOverEtaBins_Coarse,
                         len(SFs_AvgInPtQuant_AvgOverEtaBins), SFs_AvgInPtQuant_AvgOverEtaBins
                     ))

                SFs_AvgInPtQuant_AvgOverEtaBins[idxStart : idxEnd] = [avgSF_iQuant] * (idxEnd - idxStart)
                if PrintLevel >= 10:
                    print("SFs_AvgInPtQuant_AvgOverEtaBins after ({}): {}".format(
                        len(SFs_AvgInPtQuant_AvgOverEtaBins), SFs_AvgInPtQuant_AvgOverEtaBins
                    ))
            # ----------------------------------------------------------------------------------------------------------

                    
            sDir1 = "plots_makeLUTs_%s_%s_%s_%s" % (sLUTVersion, jetShape, PUSAlgo, calibSFLable)
            sDir1 = sDir1.replace("(", "_")
            sDir1 = sDir1.replace(")", "_")
            sDir1 = sDir1.replace("__", "_")
            sDirPlots_JECLUT = '%s/JECLUT' % (sDir1)
            if not os.path.exists(sDir1):
                os.makedirs(sDir1)
            if not os.path.exists(sDirPlots_JECLUT):
                os.makedirs(sDirPlots_JECLUT)            


            ## plot SF_avgEtaBins in pT quantiles --------------------------------------------------    
            if 1==1:
                fig, ax = plt.subplots()  # Create a figure containing a single axes.

                ax.plot(Pt_Bins_Read, SFs_AvgOverEtaBins, label='Avg. over all iEta', color='k') # 'o'

                ax.plot(Pt_Bins_Read, SFs_AvgInPtQuant_AvgOverEtaBins, label='Avg. over all iEta, pT quantile', color='r') # 'o'
                ax.plot(Pt_AvgInPtQuant, SFs_AvgInPtQuant_AvgOverEtaBins_Coarse, label='Avg. over all iEta, pT quantile - 2', color='None', marker='o', markerfacecolor='b', markeredgecolor='b') # 'o'

                for iPt_quantile in range(len(Pt_quantiles)):
                    label1 = "Compressed pT bin edge" if iPt_quantile == 0 else None
                    ax.axvline(Pt_quantiles[iPt_quantile], color='dimgray', linestyle='dotted', label=label1)
                
                ax.set_xlabel('L1 jet pT [GeV]')
                ax.set_ylabel('SF')
                #ax.set_title("Simple Plot")  # Add a title to the axes.
                ax.legend()  # Add a legend.

                ax.set_xscale('log', base=10)
                #plt.xscale('log', base=2)

                # Set axes limit
                #plt.xlim(10, 60)
                
                #plt.show()
                plt.savefig('%s/ietaAll_SFsOriginal_2.png' % (sDir1))
            # -------------------------------------------------------------------------------------




            sDir_LUT = "LUTs_%s/%s_%s_%s" % (sLUTVersion, jetShape, PUSAlgo, calibSFLable)
            sDir_LUT = sDir_LUT.replace("(","_")
            sDir_LUT = sDir_LUT.replace(")","_")
            sDir_LUT = sDir_LUT.replace("__","_")
            if not os.path.exists(sDir_LUT):
                os.makedirs(sDir_LUT)

            fOut_LUT_pt_compress        = open("%s/%s" % (sDir_LUT, sFOut_LUT_pt_compress),        'w')
            fOut_LUT_eta_compress       = open("%s/%s" % (sDir_LUT, sFOut_LUT_eta_compress),       'w')
            fOut_LUT_calib_compress_HEX = open("%s/%s" % (sDir_LUT, sFOut_LUT_calib_compress_HEX), 'w')
            fOut_LUT_calib_compress     = open("%s/%s" % (sDir_LUT, sFOut_LUT_calib_compress),     'w')

            # write headers
            fOut_LUT_eta_compress.write('# MP ieta compression LUT\n')
            fOut_LUT_eta_compress.write('# Converts abs(MP ieta) (6 bits) into 4 bit index\n')
            fOut_LUT_eta_compress.write('# This is NOT calo ieta\n')
            fOut_LUT_eta_compress.write('# anything after # is ignored with the exception of the header\n')
            fOut_LUT_eta_compress.write('# the header is first valid line starting with #<header> versionStr nrBitsAddress nrBitsData </header>\n')
            if EtaCompressedLUT:
                fOut_LUT_eta_compress.write('#<header> v1 6 4 </header>\n')
            else:
                fOut_LUT_eta_compress.write('#<header> v1 6 6 </header>\n') # eta uncompresses LUT

            fOut_LUT_pt_compress.write('# PT compression LUT\n')
            #fOut_LUT_pt_compress.write('# maps 8 bits to 4 bits\n')
            fOut_LUT_pt_compress.write('# maps 8 bits to %d bits\n' % (nBitsForPtComp))
            fOut_LUT_pt_compress.write('# the 1st column is the integer value after selecting bits 1:8\n')
            fOut_LUT_pt_compress.write('# anything after # is ignored with the exception of the header\n')
            fOut_LUT_pt_compress.write('# the header is first valid line starting with #<header> versionStr nrBitsAddress nrBitsData </header>\n')
            #fOut_LUT_pt_compress.write('#<header> v1 8 4 </header>\n')
            fOut_LUT_pt_compress.write('#<header> v1 8 %d </header>\n' % (nBitsForPtComp))

            if MakeJECLUT:
                fOut_LUT_calib_compress_HEX.write("# address to addend+multiplicative factor LUT\n")
                fOut_LUT_calib_compress_HEX.write("# maps " + str(nBitsCompBin) + " bits to 18 bits\n")
                fOut_LUT_calib_compress_HEX.write("# 18 bits = (addend<<%d) + multiplier)\n" % (nBitsJECMultiplier))
                fOut_LUT_calib_compress_HEX.write("# addend is signed 8 bits, multiplier is 10 bits\n")
                fOut_LUT_calib_compress_HEX.write("# anything after # is ignored with the exception of the header\n")
                fOut_LUT_calib_compress_HEX.write("# the header is first valid line starting with #<header> versionStr nrBitsAddress nrBitsData </header>\n")
                fOut_LUT_calib_compress_HEX.write("#<header> v1 " + str(nBitsCompBin) + " 18 </header>\n")



            ## Update SFs in each IEta bin to have in pT_quantiles --------------------------------
            if PrintLevel >= 2:
                print("Update SFs in each IEta bin to have in pT_quantiles".format())
                #print("\n\n data_calibSFs @10 : {} \n\n".format(data_calibSFs))
                
            data_calibSFs_compressed       = OrderedDict()
            #data_calibSFs_compressed_check = data_calibSFs.copy()
            data_calibSFs_compressed_check = copy.deepcopy(data_calibSFs)

            data_calibSFs_compressedInEtaAndPt   = OrderedDict()
            data_calibSFs_compressedInPt         = OrderedDict()
            data_calibSFs_uncompress             = copy.deepcopy(data_calibSFs) 
            
            IEtaBin_forLUT_col0 = None

            CaloToolMPEtaBinsMerge = None
            if EtaCompressedLUT:
                if EtaCompressedLUTVersion == 'v2018':
                    CaloToolMPEtaBinsMerge = CaloToolMPEtaBinsMerge_forEtaCompressedLUT_2018 # use 2018 Eta compressed bins
                if PrintLevel >= 0:
                    print("EtaCompressedLUTVersion: {}, CaloToolMPEtaBinsMerge_forEtaCompressedLUT_2018".format(EtaCompressedLUTVersion))
                    
                elif EtaCompressedLUTVersion == 'v2022ChunkyDonut':
                    CaloToolMPEtaBinsMerge = CaloToolMPEtaBinsMerge_forEtaCompressedLUT_2022_ChunkyDonut
                    if 'RawPUS_phiDefault' in PUSAlgosAll:
                        print("Mismatch in {}  and {} \t\t\t **** ERROR **** \n".format(PUSAlgosAll, EtaCompressedLUTVersion))
                        exit(0)
                    if PrintLevel >= 0:
                        print("EtaCompressedLUTVersion: {}, CaloToolMPEtaBinsMerge_forEtaCompressedLUT_2022_ChunkyDonut".format(EtaCompressedLUTVersion))
                        
                elif EtaCompressedLUTVersion == 'v2022PhiRing':
                    CaloToolMPEtaBinsMerge = CaloToolMPEtaBinsMerge_forEtaCompressedLUT_2022_PhiRing
                    if 'RawPUS' in PUSAlgosAll:
                        print("Mismatch in {}  and {} \t\t\t **** ERROR **** \n".format(PUSAlgosAll, EtaCompressedLUTVersion))
                        exit(0)
                    if PrintLevel >= 0:
                        print("EtaCompressedLUTVersion: {}, CaloToolMPEtaBinsMerge_forEtaCompressedLUT_2022_PhiRing".format(EtaCompressedLUTVersion))
                        
                elif EtaCompressedLUTVersion == 'v2022Merged':
                    CaloToolMPEtaBinsMerge = CaloToolMPEtaBinsMerge_forEtaCompressedLUT_2022_ChunkyDonutAndPhiRingMerged
                    if PrintLevel >= 0:
                        print("EtaCompressedLUTVersion: {}, CaloToolMPEtaBinsMerge_forEtaCompressedLUT_2022_ChunkyDonutAndPhiRingMerged".format(EtaCompressedLUTVersion))
            else:
                CaloToolMPEtaBinsMerge = CaloToolMPEtaBinsMerge_forEtaUncompressedLUT
                if PrintLevel >= 0:
                    print("EtaCompressedLUTVersion: {}, CaloToolMPEtaBinsMerge_forEtaUncompressedLUT".format(EtaCompressedLUTVersion))
                    
            if PrintLevel >= 0:
                print("CaloToolMPEtaBinsMerge: {}".format(CaloToolMPEtaBinsMerge))


                

            if PrintLevel >= 0:
                print("\n\nCalculate pT-Eta compressed LUTS:")
                
            for CaloToolMPEtaBin_forLUT, CaloToolMPEtaRange in CaloToolMPEtaBinsMerge.items(): # run on each iEtaBinRange

                ## Plots SFs used in making LUTs
                fig, ax = plt.subplots()  # Create a figure containing a single axes.
                
                data_calibSFs_compressedInEtaAndPt[CaloToolMPEtaBin_forLUT]   = OrderedDict()
                data_calibSFs_compressedInPt[CaloToolMPEtaBin_forLUT]         = OrderedDict()
                SFs_array = np.zeros(shape=(len(CaloToolMPEtaRange), NCompPtBins ))

                PtQuantiles_info = OrderedDict()
                if PrintLevel >= 0:
                    print(f"{' '*4}CaloToolMPEtaBin_forLUT: {CaloToolMPEtaBin_forLUT}")
                 
                IEtaBin_forLUT_col1 = CaloToolMPEtaBin_forLUT
                i_ = 0
                for CaloTool_mpEta in CaloToolMPEtaRange: # runs each iEta in iEtaBinRange
                    idx_CaloTool_mpEta_inRange = CaloToolMPEtaRange.index(CaloTool_mpEta)
                    IEtaBin = convert_CaloToolMPEta_to_IEta(CaloTool_mpEta)
                    IEtaBin_forLUT_col0 = CaloTool_mpEta
                    data_calibSFs_compressedInPt[CaloToolMPEtaBin_forLUT][IEtaBin]         = OrderedDict()
                    data_calibSFs_compressed[IEtaBin]       = OrderedDict()
                    
                    if CaloToolMPEtaBin_forLUT == 0 and idx_CaloTool_mpEta_inRange == 0:
                        fOut_LUT_eta_compress.write('0 0  # dummy ieta_bin\n')
                    fOut_LUT_eta_compress.write('%d %d  # ieta %d\n' % (int(IEtaBin_forLUT_col0), int(IEtaBin_forLUT_col1), IEtaBin))

                    if PrintLevel >= 0:
                        print(f"{' '*8}CaloTool_mpEta: {CaloTool_mpEta}, IEtaBin: {IEtaBin},  idx_CaloTool_mpEta_inRange: {idx_CaloTool_mpEta_inRange} ")
                    
                    for iPtQuant in range(NCompPtBins):    # run on iPtQuant in iEta               
                        idxPtStart                  = 0             if iPtQuant == 0               else idx_PtQuantiles[(iPtQuant-1)]
                        idxPtEnd                    = nPt_Bins_Read if iPtQuant == (NCompPtBins-1) else idx_PtQuantiles[(iPtQuant)]
                        Pt_Bins_Read_iPtQuant_Full  = Pt_Bins_Read[idxPtStart : idxPtEnd] # for lut_pt_compression                        
                        if iPtQuant > 0 and Pt_Bins_Read[idxPtStart] < calibSF_L1JetPtRange[0]: 
                            # except pT=0 bin, use pT range as set in calibSF_L1JetPtRange. Note: Only for average SF calculation
                            idxPtStart = Pt_Bins_Read.index( calibSF_L1JetPtRange[0] )
                        Pt_Bins_Read_iPtQuant     = Pt_Bins_Read[idxPtStart : idxPtEnd]

                        #SFs_iPtQuant              = list(data_calibSFs[IEtaBin].values())[idxPtStart : idxPtEnd]
                        SFs_iPtQuant              = list(data_calibSFs[IEtaBin].values()).copy()[idxPtStart : idxPtEnd]

                                                
                        PtMin_iPtQuant            = Pt_Bins_Read_iPtQuant[0]
                        PtMax_iPtQuant            = Pt_Bins_Read_iPtQuant[-1]
                        avgPt_iPtQuant            = statistics.mean( Pt_Bins_Read_iPtQuant )
                        avgSF_iPtQuant            = statistics.mean( SFs_iPtQuant )
                        if separatePtBinForLowPt and iPtQuant == 0:  
                            # ignore SF(pT=0)=1 and set SF(pT<=15 GeV) = SF(pT=15GeV)
                            avgSF_iPtQuant        = SFs_iPtQuant[1]
                            avgPt_iPtQuant        = Pt_Bins_Read_iPtQuant[-1] # calculate SF_inBits at bin higher edge (pT=15 GeV)

                        if PrintLevel >= 4:
                            print(f"{' '*12}{iPtQuant=}, {idxPtStart=}, {idxPtEnd=}, {Pt_Bins_Read_iPtQuant=}, {SFs_iPtQuant=},   {avgSF_iPtQuant=}  {avgPt_iPtQuant=}")
                            
                        #data_calibSFs_compressed[IEtaBin][iPtQuant] = avgSF_iPtQuant
                        data_calibSFs_compressed[IEtaBin][avgPt_iPtQuant] = avgSF_iPtQuant
                        data_calibSFs_compressedInPt[CaloToolMPEtaBin_forLUT][IEtaBin][avgPt_iPtQuant] = avgSF_iPtQuant

                        #data_calibSFs_compressed_check[IEtaBin][idxPtStart : idxPtEnd] = [avgSF_iPtQuant] * (idxPtEnd - idxPtStart)
                        for idx_PtUncompressesBin in range(idxPtStart, idxPtEnd):
                            data_calibSFs_compressed_check[IEtaBin][idx_PtUncompressesBin] = avgSF_iPtQuant


                        # fOut_LUT_pt_compress ---------------------------
                        PtBin_forLUT = iPtQuant                    
                        if CaloToolMPEtaBin_forLUT == 0 and idx_CaloTool_mpEta_inRange == 0:
                            for Pt in Pt_Bins_Read_iPtQuant_Full:
                                fOut_LUT_pt_compress.write('%d %d\n' % (int(Pt), int(PtBin_forLUT)))
                        # ------------------------------------------------


                        if idx_CaloTool_mpEta_inRange == 0:
                            PtQuantiles_info[iPtQuant] = OrderedDict([
                                ('PtMin_iPtQuant', int(PtMin_iPtQuant + 0.5)),
                                ('PtMax_iPtQuant', int(PtMax_iPtQuant + 0.5)),
                                ('avgPt_iPtQuant', int(avgPt_iPtQuant + 0.5)),
                                ('PtBin_forLUT',   int(PtBin_forLUT))
                            ])


                        '''
                        sComments = ""
                        if iPtQuant == 0:
                            sComments = "   # ieta %d, pt %d" % (IEtaBin, PtBin_forLUT)
                        fOut_LUT_calib_compress.write('%d %d %d %d %d %d %f%s\n' % (int(IEtaBin), int(PtMin_iPtQuant + 0.5), int(PtMax_iPtQuant + 0.5), int(avgPt_iPtQuant + 0.5), int(IEtaBin_forLUT_col1), int(PtBin_forLUT), avgSF_iPtQuant, sComments ))
                        '''

                        
                    if PrintLevel >= 10:
                        print("\ndata_calibSFs_compressed[{}]: ({}) {} \ndata_calibSFs_compressed_check[{}]: ({}) {}".format(
                            IEtaBin, len(data_calibSFs_compressed[IEtaBin].keys()), data_calibSFs_compressed[IEtaBin],
                            IEtaBin, len(data_calibSFs_compressed_check[IEtaBin].keys()), data_calibSFs_compressed_check[IEtaBin]
                        ))

                    if PrintLevel >= 20:
                        print("\n\n data_calibSFs[{}] @11  : {} \n\n".format(IEtaBin, data_calibSFs[IEtaBin]))

                        
                    SFs_array[idx_CaloTool_mpEta_inRange] = np.array( list(data_calibSFs_compressedInPt[CaloToolMPEtaBin_forLUT][IEtaBin].values()) )

                    ax.plot(data_calibSFs[IEtaBin].keys(), data_calibSFs[IEtaBin].values(), label='uncompressed SFs, iEta=%s' % (IEtaBin), color=colors_matplotlib[i_])
                    ax.plot(data_calibSFs_compressed[IEtaBin].keys(), data_calibSFs_compressed[IEtaBin].values(), color=colors_matplotlib[i_], marker='o', linestyle='')
                    i_ = i_ + 1
                    



                    
                ## average SF over all eta bins
                SFs_AvgOverEtaBinsRange = np.average(SFs_array, axis=0)
                if PrintLevel >= 10:
                    print("\nCaloToolMPEtaBin_forLUT: {}, CaloToolMPEtaRange: {}, SFs_AvgOverEtaBinsRange ({}): {}".format(CaloToolMPEtaBin_forLUT, CaloToolMPEtaRange, len(SFs_AvgOverEtaBinsRange), SFs_AvgOverEtaBinsRange))

                if PrintLevel >= 0: #2:
                    print("\nPrint SFs in CaloToolMPEtaBin_forLUT {},  CaloToolMPEtaRange {}".format(CaloToolMPEtaBin_forLUT, CaloToolMPEtaRange))
                    if EtaCompressedLUT:
                        for CaloTool_mpEta in CaloToolMPEtaRange: # runs each iEta in iEtaBinRange
                            IEtaBin = convert_CaloToolMPEta_to_IEta(CaloTool_mpEta)
                            print("CaloTool_mpEta {}, IEtaBin {} SFs_0 ({}): {}".format(CaloTool_mpEta, IEtaBin, len(list(data_calibSFs_compressed[IEtaBin].values())), data_calibSFs_compressed[IEtaBin].values()))
                            #print("CaloTool_mpEta {}, IEtaBin {} SFs_1 ({}): {}".format(CaloTool_mpEta, IEtaBin, len(list(data_calibSFs_compressedInPt[CaloToolMPEtaBin_forLUT][IEtaBin].values())), data_calibSFs_compressedInPt[CaloToolMPEtaBin_forLUT][IEtaBin].values()))
                    print("SFs_AvgOverEtaBinsRange ({}): {}".format(len(SFs_AvgOverEtaBinsRange), SFs_AvgOverEtaBinsRange))

                
                for iPtQuant, PtQuantile_info in PtQuantiles_info.items():
                    PtMin_iPtQuant = PtQuantile_info['PtMin_iPtQuant']
                    PtMax_iPtQuant = PtQuantile_info['PtMax_iPtQuant']
                    avgPt_iPtQuant = PtQuantile_info['avgPt_iPtQuant']
                    avgSF_iPtQuant = SFs_AvgOverEtaBinsRange[iPtQuant]

                    if avgSF_iPtQuant < JECLUTSF_boundary[0]: avgSF_iPtQuant = JECLUTSF_boundary[0]
                    if avgSF_iPtQuant > JECLUTSF_boundary[1]: avgSF_iPtQuant = JECLUTSF_boundary[1]

                    data_calibSFs_compressedInEtaAndPt[CaloToolMPEtaBin_forLUT][avgPt_iPtQuant] = avgSF_iPtQuant

                    PtBin_forLUT = iPtQuant 
                    sComments  = ""
                    sSFClosure = ""
                    if iPtQuant == 0:
                        #sComments = "   # ieta %d, pt %d" % (IEtaBin, PtBin_forLUT)
                        sComments = "   # ieta {}, pt {}".format(CaloToolMPEtaRange, PtBin_forLUT)
                        if int(CaloToolMPEtaBin_forLUT) == 0:
                            fOut_LUT_calib_compress.write('# columns \n')
                            fOut_LUT_calib_compress.write('# %s %s %s %s %s %s %s\n' % ('CaloToolMPEtaBin', 'PtMin', 'PtMax', 'avgPt', 'IEtaBin_forLUT_col1', 'PtBin_forLUT', 'SF' ))
                    #fOut_LUT_calib_compress.write('%d %d %d %d %d %d %f%s\n' % (int(CaloToolMPEtaBin_forLUT), int(PtMin_iPtQuant + 0.5), int(PtMax_iPtQuant + 0.5), int(avgPt_iPtQuant + 0.5), int(IEtaBin_forLUT_col1), int(PtBin_forLUT), avgSF_iPtQuant, sComments ))

                    ## Calculate JEC SFs in HEX
                    if MakeJECLUT:
                        iEta_ = int(CaloToolMPEtaBin_forLUT)
                        PtMin_iQuant_ = int(PtMin_iPtQuant + 0.5)
                        PtMax_iQuant_ = int(PtMax_iPtQuant + 0.5)
                        Pt_ = int(avgPt_iPtQuant + 0.5)

                        etaBin_ = int(IEtaBin_forLUT_col1) # 4
                        ptBin_  = int(PtBin_forLUT)  # 5
                        compBin_       = (etaBin_ << nBitsPtComp) | ptBin_;

                        SFInDecimal_ = avgSF_iPtQuant

                        jetHwPt_       = Pt_ * 2; # jetHwPt = 2*jetPtInGeV https://github.com/cms-sw/cmssw/blob/master/L1Trigger/L1TCalorimeter/src/firmware/Stage2Layer2JetAlgorithmFirmwareImp1.cc#L717 . Applying a factor of 2 does not affect SFInDecimal --> SFInBits calculation here.
                        jetHwPtMin_    = PtMin_iQuant_ * 2; 
                        jetHwPtMax_    = PtMax_iQuant_ * 2;

                        SFInBits_    = calculateSFInBits_addendZero(jetHwPt_, SFInDecimal_);

                        jetHwPtCorr_ = calibrateJet(jetHwPt_, SFInBits_)
                        SFInDecimal_fromJECLUT_ = jetHwPtCorr_ / jetHwPt_ if jetHwPt_ > 0 else 1


                        sComment_JECLUTHEX_ = ""
                        if ptBin_ == 0:
                            sComment_JECLUTHEX_ = "  # compresses ieta bin " + str(etaBin_) + ", pt " + str(ptBin_)
                        sSFClosure = "    %d %f %f    " % (SFInBits_, SFInDecimal_fromJECLUT_, SFInDecimal_fromJECLUT_ / SFInDecimal_)

                        fOut_LUT_calib_compress_HEX.write('%d %d%s\n' % (compBin_, SFInBits_, sComment_JECLUTHEX_))

                        compBin_last_JECLUT = compBin_

                    # JEc LUT Decimal
                    fOut_LUT_calib_compress.write('%d %d %d %d %d %d %f%s%s\n' % (int(CaloToolMPEtaBin_forLUT), int(PtMin_iPtQuant + 0.5), int(PtMax_iPtQuant + 0.5), int(avgPt_iPtQuant + 0.5), int(IEtaBin_forLUT_col1), int(PtBin_forLUT), avgSF_iPtQuant, sSFClosure,sComments ))

                


                ax.plot(data_calibSFs_compressedInEtaAndPt[CaloToolMPEtaBin_forLUT].keys(), data_calibSFs_compressedInEtaAndPt[CaloToolMPEtaBin_forLUT].values(),
                        label='pT-eta compressed SFs', color='black', marker='+', markersize=14, linestyle='', linewidth=3)
                    
                for iPt_quantile in range(len(Pt_quantiles)):
                    label1 = "Compressed pT bin edge" if iPt_quantile == 0 else None
                    ax.axvline(Pt_quantiles[iPt_quantile], color='dimgray', linestyle='dotted', label=label1)
                
                ax.set_xlabel('L1 jet pT [GeV]')
                ax.set_ylabel('SF')
                #ax.set_title("Simple Plot")  # Add a title to the axes.
                ax.legend()  # Add a legend.

                ax.set_xscale('log', base=10)

                # Set axes limit
                plt.xlim(10, 280)
                #plt.xlim(0, 40)
                
                #plt.show()
                plt.savefig('%s/SFs_ieta%dto%d.png' % (
                    sDirPlots_JECLUT, 
                    convert_CaloToolMPEta_to_IEta(CaloToolMPEtaRange[0]),
                    convert_CaloToolMPEta_to_IEta(CaloToolMPEtaRange[-1])
                    ))




            # dummy bins
            # nBinsMaxForEtaCompressionLUT
            #for IEtaBin_forLUT in range(len(EtaBins_sorted), nBinsMaxForEtaCompressionLUT):
            #    fOut_LUT_eta_compress.write('%d %d\n' % (int(IEtaBin_forLUT), int(IEtaBin_forLUT)))
            IEtaBin_forLUT_col0 += 1 # the next bin = 42
            while IEtaBin_forLUT_col0 < nBinsMaxForEtaCompressionLUT:
                IEtaBin_forLUT_col1 = IEtaBin_forLUT_col0 - 1 
                #fOut_LUT_eta_compress.write('%d %d\n' % (int(IEtaBin_forLUT_col0), int(IEtaBin_forLUT_col1)))
                fOut_LUT_eta_compress.write('%d %d\n' % (int(IEtaBin_forLUT_col0), int(0)))
                IEtaBin_forLUT_col0 += 1

            # dummy bins for JEC LUTs
            if MakeJECLUT:
                compBin_last_JECLUT = compBin_last_JECLUT + 1
                while compBin_last_JECLUT < nCompBinMax:
                    fOut_LUT_calib_compress_HEX.write('%d 0 # dummy\n' % (compBin_last_JECLUT))
                    compBin_last_JECLUT = compBin_last_JECLUT + 1

            
            fOut_LUT_pt_compress.close()
            fOut_LUT_eta_compress.close()
            if MakeJECLUT: fOut_LUT_calib_compress_HEX.close()
            fOut_LUT_calib_compress.close()

            print("Wrote %s/%s, %s/%s, %s/%s, %s/%s" % (
                sDir_LUT, sFOut_LUT_pt_compress,
                sDir_LUT, sFOut_LUT_eta_compress,
                sDir_LUT, sFOut_LUT_calib_compress_HEX,
                sDir_LUT, sFOut_LUT_calib_compress
            ))



            ## Compare data_calibSFs uncompressed SFs from bunch of iEta bins, and compare them with SFs_AvgOverAllIEtaBins -------------
            for iEtaBinsMerge in IEtaBinsMerge_forPlots:
                if PrintLevel >= 10:
                    print("iEtaBinsMerge: {}".format(iEtaBinsMerge))

                fig, ax = plt.subplots()  # Create a figure containing a single axes.

                x_array_tmp = None
                for Eta_bin0 in iEtaBinsMerge:
                    Eta_bin = int(Eta_bin0)
                    if Eta_bin not in list(data_calibSFs.keys()):
                        print("Eta {} not in data_calibSFs ".format(Eta_bin))
                        continue
                    
                    if PrintLevel >= 10:                    
                        print("Eta {}: ".format(Eta_bin))
                        print(" {} {}".format(len(list(data_calibSFs[Eta_bin].keys())), list(data_calibSFs[Eta_bin].keys())))
                        print(" {} {}".format(len(list(data_calibSFs[Eta_bin].values())), list(data_calibSFs[Eta_bin].values())))

                    # plot
                    x_array_tmp = np.array( list(data_calibSFs[Eta_bin].keys()) )
                    y_array_tmp = np.array( list(data_calibSFs[Eta_bin].values()) )

                    ax.plot(x_array_tmp,y_array_tmp, label='iEta=%s' % (Eta_bin)) # 'o'


                y_array_tmp = SFs_AvgOverEtaBins
                ax.plot(x_array_tmp,y_array_tmp, label='Avg. over all iEta', color='k') # 'o'

                #ax.plot(50,1, 50,1.2, color='k')
                #plt.plot(50,1, 50,1.2, color='k')
                #ax.axvline(50, color='k', linestyle='dashdot')
                #ax.axvline(50, color='k', linestyle=(0, (1, 10)), color='#FFFF00')
                #ax.axvline(50, color='#228B22', linestyle='dotted')
                #ax.axvline(50, color='mediumblue', linestyle='dotted')
                #x1, y1 = [50, 50], [1, 1.3]
                #ax.plot(x1, y1, color='k', ls='--')

                for iPt_quantile in range(len(Pt_quantiles)):
                    label1 = "Compressed pT bin edge" if iPt_quantile == 0 else None
                    ax.axvline(Pt_quantiles[iPt_quantile], color='dimgray', linestyle='dotted', label=label1)
               
                
                ax.set_xlabel('L1 jet pT [GeV]')
                ax.set_ylabel('SF')
                #ax.set_title("Simple Plot")  # Add a title to the axes.
                ax.legend()  # Add a legend.

                    
                #plt.show()
                plt.savefig('%s/ieta%dto%d_SFsOriginal.png' % (sDir1,iEtaBinsMerge[0], iEtaBinsMerge[-1]))
            #----------------------------------------------------------------------------------------------------------
                

            ### -----------------------------
            for Eta_bin in EtaBins_sorted:
                if Eta_bin == 29: continue
                
                fig, ax = plt.subplots()  # Create a figure containing a single axes.
                #plt.figure(figsize=(5,5))
                #fig = plt.figure(figsize=(18, 16), dpi=300)
                #fig.figure(figsize=(2, 2))
                
                # plot
                Pt_array_original  = np.array( list(data_calibSFs[Eta_bin].keys()) )
                SFs_array_original = np.array( list(data_calibSFs[Eta_bin].values()) )
                #
                SFs_array_AbgOverPtBinQuant = np.array( list(data_calibSFs_compressed_check[Eta_bin].values()) )
                #
                Pt_array_PtQuant  = np.array( list(data_calibSFs_compressed[Eta_bin].keys()) )
                SFs_array_AbgOverPtBinQuant_coarseBins = np.array( list(data_calibSFs_compressed[Eta_bin].values()) )
                

                ax.plot(Pt_array_original, SFs_array_original, label='iEta=%s' % (Eta_bin)) # 'o'

                ax.plot(Pt_array_original, SFs_array_AbgOverPtBinQuant, label='iEta=%s pT compressed - 1' % (Eta_bin), color='r') # 'o'

                #ax.plot(Pt_array_PtQuant, SFs_array_AbgOverPtBinQuant_coarseBins, label='iEta=%s pT compressed - 2' % (Eta_bin), color='None', marker='o', markerfacecolor='b', markeredgecolor='b') # 'o'
                

                for iPt_quantile in range(len(Pt_quantiles)):
                    label1 = "Compressed pT bin edge" if iPt_quantile == 0 else None
                    ax.axvline(Pt_quantiles[iPt_quantile], color='dimgray', linestyle='dotted', label=label1)
                
                ax.set_xlabel('L1 jet pT [GeV]')
                ax.set_ylabel('SF')
                #ax.set_title("Simple Plot")  # Add a title to the axes.
                ax.legend()  # Add a legend.

                ax.set_xscale('log', base=10)
                #plt.xscale('log', base=2)

                # Set axes limit
                plt.xlim(10, 280)
                
                #plt.show()
                plt.savefig('%s/SFs_pTCompressed_ieta%d.png' % (sDir1, Eta_bin))



            ## Compare Eta-Pt compresses SFs ---------------------------------------------------------------------------
            for CaloToolMPEtaBin_forLUT, CaloToolMPEtaRange in CaloToolMPEtaBinsMerge.items(): # run on each iEtaBinRange
                fig, ax = plt.subplots()  # Create a figure containing a single axes.

                for CaloTool_mpEta in CaloToolMPEtaRange: # runs each iEta in iEtaBinRange
                    idx_CaloTool_mpEta_inRange = CaloToolMPEtaRange.index(CaloTool_mpEta)
                    IEtaBin = convert_CaloToolMPEta_to_IEta(CaloTool_mpEta)

                    # plot
                    x_array_tmp = np.array( list(data_calibSFs_compressed[IEtaBin].keys()) )
                    y_array_tmp = np.array( list(data_calibSFs_compressed[IEtaBin].values()) )

                    ax.plot(x_array_tmp,y_array_tmp, label='iEta=%s' % (IEtaBin), marker='o', markerfacecolor=colors_matplotlib[idx_CaloTool_mpEta_inRange], markeredgecolor=colors_matplotlib[idx_CaloTool_mpEta_inRange], color='None') # 'o'

                y_array_tmp = np.array( list(data_calibSFs_compressedInEtaAndPt[CaloToolMPEtaBin_forLUT].values()) )
                ax.plot(x_array_tmp,y_array_tmp, label='Avg. over iEta range', color='None', marker='^', markerfacecolor='k') # 'o'

                for iPt_quantile in range(len(Pt_quantiles)):
                    label1 = "Compressed pT bin edge" if iPt_quantile == 0 else None
                    ax.axvline(Pt_quantiles[iPt_quantile], color='dimgray', linestyle='dotted', label=label1)
               
                
                ax.set_xlabel('L1 jet pT [GeV]')
                ax.set_ylabel('SF')
                #ax.set_title("Simple Plot")  # Add a title to the axes.
                ax.legend()  # Add a legend.

                ax.set_xscale('log', base=10)

                # Set axes limit
                plt.xlim(10, 280)
                
                #plt.show()
                plt.savefig('%s/EtaPtComptessesSF_iEta%02d_to_%02d.png' % (sDir1, convert_CaloToolMPEta_to_IEta(CaloToolMPEtaRange[0]), convert_CaloToolMPEta_to_IEta(CaloToolMPEtaRange[-1]) ))
            #----------------------------------------------------------------------------------------------------------
 

            for iEtaBinsMerge in IEtaBinsMerge_forPlots_1:

                fig, ax = plt.subplots()  # Create a figure containing a single axes.

                x_array_tmp = None
                for IEtaBin in iEtaBinsMerge:
                    if IEtaBin == 29: continue

                    # plot
                    x_array_tmp = np.array( list(data_calibSFs_compressed[IEtaBin].keys()) )
                    y_array_tmp = np.array( list(data_calibSFs_compressed[IEtaBin].values()) )

                    ax.plot(x_array_tmp,y_array_tmp, label='iEta=%s' % (IEtaBin)) # 'o'

                
                ax.set_xlabel('L1 jet pT [GeV]')
                ax.set_ylabel('SF')
                #ax.set_title("Simple Plot")  # Add a title to the axes.
                ax.legend()  # Add a legend.

                ax.set_xscale('log', base=10)

                # Set axes limit
                plt.xlim(10, 280)
                
                #plt.show()
                plt.savefig('%s/PtComptessesSF_iEta%02d_to_%02d.png' % (sDir1, iEtaBinsMerge[0], iEtaBinsMerge[-1] ))
            #----------------------------------------------------------------------------------------------------------
 

            for iEtaBinsMerge in IEtaBinsMerge_forPlots_2:

                fig, ax = plt.subplots()  # Create a figure containing a single axes.

                x_array_tmp = None
                for IEtaBin in iEtaBinsMerge:
                    if IEtaBin == 29: continue

                    # plot
                    x_array_tmp = np.array( list(data_calibSFs_compressed[IEtaBin].keys()) )
                    y_array_tmp = np.array( list(data_calibSFs_compressed[IEtaBin].values()) )

                    ax.plot(x_array_tmp,y_array_tmp, label='iEta=%s' % (IEtaBin)) # 'o'

                
                ax.set_xlabel('L1 jet pT [GeV]')
                ax.set_ylabel('SF')
                #ax.set_title("Simple Plot")  # Add a title to the axes.
                ax.legend()  # Add a legend.

                ax.set_xscale('log', base=10)

                # Set axes limit
                plt.xlim(10, 280)
                
                #plt.show()
                plt.savefig('%s/PtComptessesSF_v2_iEta%02d_to_%02d.png' % (sDir1, iEtaBinsMerge[0], iEtaBinsMerge[-1] ))
            #----------------------------------------------------------------------------------------------------------
 
                        

                

            print("Done !!!")

            
