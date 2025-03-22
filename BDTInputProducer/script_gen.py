
import time
import os
import sys
import math 
from array import array
import ROOT as R
R.gROOT.SetBatch(False)  ## Don't print histograms to screen while processing

from collections import OrderedDict, defaultdict
import csv
import glob
import json
import urllib.request
import argparse
import numpy as np 
from operator import xor
import uproot
import awkward as ak


PRT_EVT  = 1000  ## Print every Nth event
MAX_EVT  = -1     ## Number of events to process per chain
VERBOSE  = False  ## Verbose print-out
PrintLevel = 0
JetClustByHand =  True # True ## Run jet clustering by hand
#JetShapes = ['Default', '9x9', '7x9', '5x9', '3x9', '3x9_plus_0.5_times_9x9']
JetShapes = ['Default', '9x9', ]
JetShapesForML = ['9x9', ]
JetShapesType2 = [] # ['L1TauDefault']
#PUSAlgosAll      = ['L1JDefault', 'Raw', 'RawPUS', 'RawPUS_phiRingMin4', 'RawPUS_phiRingSide4', 'RawPUS_phiRingAdjacent', 'RawPUS_phiDefault']
PUSAlgosAll      = ['Raw', 'RawPUS', 'RawPUS_phiDefault']
PUSAlgosSelected = [] #['Raw', 'RawPUS', 'RawPUS_phiDefault']
PUSAlgosAllType2 = [] # ['Et', 'RawEt']
MatchEmulatedJetsWithUnpacked = False
HLT_Triggers_Required = [
    'IsoMu24_OneProng32' # HLT_IsoMu24_v15
]
# SingleJet180, IsoMu24_OneProng32
TrigThshs_OffMuPt = [ 24 ] # For e.g. for IsoMu24: [ 24 ], for DiMu24: [24, 24], for Mu24_Mu20: [24, 20]
TrigThshs_OffJetPt = [ 180 ] # For e.g. for SingleJet180: [ 180 ], for DiJet180: [180, 180], for Jet180_Jet120: [180, 120]

#GoldenJSONForData_list=["Cert_Collisions2022_eraG_362433_362760_Golden.json"]
GoldenJSONForData_list= ["https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions24/Cert_Collisions2024_378981_386951_Golden.json"]
useCutGenNVtxEq0 = False # Set False. Only for troubleshoot perfose. When set True: analyze (GEN.nVtx == 0) events from SinglePhoton_EpsilonPU sample to trouble-shoot high SFs in iEta 28
#offlineJetType = 'PUPPI' # 'CHS', 'PUPPI'  offlineCHSJet, offlinePUPPIJet. Set it as a command line argument 
nJetFilters = 14

runMode = 'makeInputForML' # '', 'CalCalibSF', 'CalibJetByHand', 'makeInputForML', 'trbshtPhiRingPUS'
# 'test'           # run on L1nano_*_1.root nano for tests
# ''               # 1st round to make jet resolution plots
# 'CalCalibSF'     # set true to fill PFjetPt vs L1jetPt histograms to calculate calibration SFs
# 'CalibJetByHand' # apply CaliSF to JetsByHand
# 'makeInputForML' # write jet information into .csv file to be input to MachineLearning
# 'makePUHisto'    # run quick to make PU histograms
# 'trbshtPhiRingPUS' # to troubleshoot PhiRing PU Et
sOutExt = ""


dataErasRunRange = {
    '2022C': [356426, 357482],
    '2022D': [357538, 359017], # [357538, 357733,  357734, 357930,  358381, 359017]
    '2022E': [359045, 360327], 
    '2022F': [360335, 362167], 
    '2022G': [362362, 362760], 
    #
    '2023B': [366386, 367079],
    '2023C': [367094, 369802],
    '2023D': [369827, 371227],
    #
    '2024A': [378919, 378968],
    '2024B': [378981, 379391],
    '2024C': [379415, 380238],
    '2024D': [380255, 380947],
    '2024E': [380949, 386753],
    # '2024F': [380949, 389999],
    # '2024G': [380949, 389999],
    # '2024H': [380949, 389999],
    '2024I': [386753, 389999],
    
}

L1TEffiTurnOn_TrigThrshs = [12.0, 35.0, 60.0, 90.0, 120.0, 180.0] # L1T eT threshold for L1T efficiency turn-on curve

MASS_MUON     = 0.1056584 # GeV https://pdglive.lbl.gov/Particle.action?node=S004&init=0
MASS_ELECTRON = 0.0005110 # GeV https://pdglive.lbl.gov/Particle.action?node=S003&init=0
SFToConvertInGeV = 0.5

# PU reweighting ----------------------------------------- 
usePUReweighting_0   = False;
#ipFilePUData = "/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v4_20210510/myStudies/run_1/2018_Data/l1analysis_def.root";
#ipFilePUData = "/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v4_20210510/myStudies/run_1/mimicMLJetRec/2018_Data/L1T_JetMET_Res_def_hadded.root"; # larger statistics
#ipFilePUData = "/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v4_20210510/myStudies/run_1/mimicMLJetRec/2018_SingleMu/L1T_JetMET_Res_def_hadded.root"; # larger statistics
ipFilePUData = "/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v4_20210510/myStudies/run_2/2018_SingleMu/IMPORTANT/L1T_JetMET_Res_def_hadded.root"
sHistoPUData = "nVtx";
ipFilePUMC   = "/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v4_20210510/myStudies/run_2/Run3_MC/IMPORTANT/L1T_JetMET_Res_def_hadded_v0.root"
sHistoPUMC   = "nVtx";

nPV_median = 55 # make avgLayer2SF plot for 50 < nPV <= nPV_median and nPV > nPV_median 

# Calib SFs -----------------------------------------------
PURangeUsedForCalibSF = "PU50to100" # "PU50to100", "PU1to25"
PtRangeForCalibSF = "Pt25To35" #  medPt, lowPt
#sipFileCalibSF = "L1T_JetMET_Layer2CalibSF_$OOTPUS.root";
#sHistoCalibSF = "IEta$ETA/h_jet_byHand_L1JetPt_vs_PFJetPt$JETSHAPE_$PUSALGORITHM_$ETA_PtAllBins_0_$L1MODE_fixBinWidthPFJetPt_ProfileAlongL1JetPt";
#sipFileCalibSF = "/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v4_20210510/myStudies/run_2/2018_SingleMu/IMPORTANT/L1TJetEtSF_2018_SingleMu_20210925.root";
#sipFileCalibSF = "/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v4_20210510/myStudies/run_2/2018_SingleMu/IMPORTANT/L1TJetEtSF_2018_SingleMu_20211002.root";
#sipFileCalibSF = "/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v4_20210510/myStudies/run_3/2018_SingleMu/IMPORTANT/L1TJetEtSF_2018_SingleMu_20211124.root";
#  sipFileCalibSF = "/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v4_20210510/myStudies/run_3/run_3_1_lowPt25to35/2018_SingleMu/IMPORTANT/L1TJetEtSF_2018_SingleMu_L1JetPt25To35_20211202.root"
#sHistoCalibSF  = "h_jet_byHand_res_vs_iEta_vs_nVtx_$JETSHAPE$PUSALGORITHM_HBEF_medPt_0_emu_$OOTPUS_PU50to100_SF_PFJByL1TJ"; # medPt: SF derived for PFJetpT in [60, 90] GeV
#  sHistoCalibSF  = "h_jet_byHand_res_vs_iEta_vs_nVtx$JETSHAPE_$PUSALGORITHM_HBEF_%s_0_emu_$OOTPUS_%s_SF_PFJByL1TJ" % (PtRangeForCalibSF, PURangeUsedForCalibSF); # medPt: SF derived for PFJetpT in [60, 90] GeV

# Read L1Jet CalibLayer2 SFs from csv files provided by Syed
icalibSF = 0 # 0, 1
calibSFLable = ['SF_RegressedTo_log(GenJetPt)_v6', 'SF_RegressedTo_log(L1JetPt)/log(GenJetPt)_v6'][icalibSF]  
sipFileCalibSF = {
    'Default': {
       'RawPUS': { # Chunky donut
           'fileName': '/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v5_20220311/run_1/IMPORTANT/L1JetLayer2SFs_Run3_QCD_Pt15to7000_TuneCP5_14TeV-pythia8_CMSSW1230pre1wPR37196_PFA1p_nVtxAll_L1JetRawPUS_v6_20220429.csv', 
           'SFLabel': ['ScaleFactor_Log(GenJet)', 'ScaleFactor_Log(L1JetRaw)/Log(GenJet)'][icalibSF],
           'L1JetPtVarName':'L1JetRawEnergy',
       },
       
       'RawPUS_phiDefault': {
           'fileName': '/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v5_20220311/run_1/IMPORTANT/L1JetLayer2SFs_Run3_QCD_Pt15to7000_TuneCP5_14TeV-pythia8_CMSSW1230pre1wPR37196_PFA1p_nVtxAll_PhiRingPUS_v6_20220429.csv', 
           'SFLabel': ['ScaleFactor_Log(GenJet)', 'ScaleFactor_Log(PhiRing)/Log(GenJet)'][icalibSF],
           'L1JetPtVarName':'PhiRingEnergy',
       },
    },


    '9x9': {
       'RawPUS': { # Chunky donut
           'fileName': '/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v5_20220311/run_1/IMPORTANT/L1JetLayer2SFs_Run3_QCD_Pt15to7000_TuneCP5_14TeV-pythia8_CMSSW1230pre1wPR37196_PFA1p_nVtxAll_L1JetRawPUS_v6_20220429.csv', 
           'SFLabel': ['ScaleFactor_Log(GenJet)', 'ScaleFactor_Log(L1JetRaw)/Log(GenJet)'][icalibSF],
           'L1JetPtVarName':'L1JetRawEnergy',
       },
       
       'RawPUS_phiDefault': {
           'fileName': '/afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v5_20220311/run_1/IMPORTANT/L1JetLayer2SFs_Run3_QCD_Pt15to7000_TuneCP5_14TeV-pythia8_CMSSW1230pre1wPR37196_PFA1p_nVtxAll_PhiRingPUS_v6_20220429.csv', 
           'SFLabel': ['ScaleFactor_Log(GenJet)', 'ScaleFactor_Log(PhiRing)/Log(GenJet)'][icalibSF],
           'L1JetPtVarName':'PhiRingEnergy',
       },
    },
    
}
#calibSF_additionalCorr = 8./7
#updateCalibSF_DefaultRaw = True # update calibSFs_additionalCorr['Default']['Raw'] depending on 
calibSFs_additionalCorr = {
    'Default': { # trigger tower Et was scalled by 0.5 while calculating l1Jet_pT by-hand and to calculate layer2 calibSFs
        'Raw': {'l1nanoChunkyDonut': 1. * 0.5,  'l1nanoPhiRing': 8./7 * 0.5 },
        'RawPUS': 1. * 0.5,
        'RawPUS_phiDefault': 8./7 * 0.5,
    },

    '9x9': {
        'Raw': 1.,
        'RawPUS': 1.,
        'RawPUS_phiDefault': 1.,
    },
    
}


# store selected events
sFOutSelectedEvents = "L1T_HCALL2Calib_stage1_SelectedEvents_tmp.txt" # "L1T_HCALL2Calib_stage1_SelectedEvents_tmp.txt"
# run on particular events
sFInEventsToRun ="" # "L1T_HCALL2Calib_stage1_SelectedEvents_1_tmp.txt"


calibSF_L1JetPtRange = [15., 250., 1.] # [<lowest pT>,  <hightest pT>,  <pT bin width>]



PT_MIN = 10 # 30   ## Minimum offline jet pT to consider
DR_MIN = 0.8 # 1.2  ## Minimum dR between offline jets considered for measurement
DR_MAX = 0.3  ## Maximum dR between L1T and offline jets for matching
DR_Jet_Ele_Min = 0.4 ## Minimum dR between PF jets and PF muons/electron
RATIO_PtEle_PtJet_Max = 0.5 ## Maximum pT fracttion overlapping electron or muon is allowed to have
PT_MAX_L1Jet = 1023 ## pt=1023.5 is assigned to l1jets with puEt=0 and jet's TP is saturated

PU_CUT = 'nVtxMin_45'  ## Pileup selection
ERA    = '2018AB'      ## Era, e.g. 2018AB or 2018D

PT_CAT = {}
PT_CAT['lowPt'] = [30,  60,   60]  ## Low pT, turn-on threshold, high pT
PT_CAT['medPt'] = [60,  75,   90]  ## Low pT, turn-on threshold, high pT   #[60,  90,   90]
PT_CAT['hiPt']  = [90, 120, 9999]  ## Low pT, turn-on threshold, high pT
#PT_CAT['modPt'] = [60,  75,   90]  ## Low pT, turn-on threshold, high pT
## for lowpt jet between 25 to 35 GeV
#PT_CAT['Below25Pt'] = [0,  15,   25]  ## Low pT, turn-on threshold, high pT
#PT_CAT['lowPt'] = [25,  30,   35]  ## Low pT, turn-on threshold, high pT
#PT_CAT['medPt'] = [35,  75,   90]  ## Low pT, turn-on threshold, high pT   #[60,  90,   90]
#PT_CAT['hiPt']  = [90, 120, 9999]  ## Low pT, turn-on threshold, high pT
## for lowpt jet between 25 to 35 GeV - version 2
PT_CAT['Ptlt25']   = [ 0,  15,   25]  ## Low pT, turn-on threshold, high pT
PT_CAT['Pt25To35'] = [25,  30,   35]  ## Low pT, turn-on threshold, high pT
PT_CAT['Pt35To60'] = [35,  55,   60]  ## Low pT, turn-on threshold, high pT
PT_CAT['Pt60To90'] = [60,  75,   90]  ## Low pT, turn-on threshold, high pT   #[60,  90,   90]
PT_CAT['Ptgt90']   = [90, 120, 9999]  ## Low pT, turn-on threshold, high pT
PtTrshForTau = {'Pt35To60': 40}

if runMode in ['CalCalibSF', 'makeInputForML']:
    PT_CAT['Below30Pt'] = [0,  16,   30]  ## Low pT, turn-on threshold, high pT

ETA_CAT = {}
ETA_CAT['HBEF'] = [0.000, 5.210]  ## Whole detector, 1 - 41
ETA_CAT['HB']   = [0.000, 1.392]  ## Trigger towers  1 - 16
ETA_CAT['HE1']  = [1.392, 1.740]  ## Trigger towers 17 - 20
ETA_CAT['HE2a'] = [1.740, 2.322]  ## Trigger towers 21 - 25
ETA_CAT['HE2b'] = [2.322, 3.000]  ## Trigger towers 26 - 28
ETA_CAT['HF']   = [3.000, 5.210]  ## Trigger towers 30 - 41

IETA_CAT = {}
IETA_CAT['HBEF'] = [ 1, 41]  ## Whole detector, 1 - 41
IETA_CAT['HB']   = [ 1, 16]  ## Trigger towers  1 - 16
IETA_CAT['HE1']  = [17, 20]  ## Trigger towers 17 - 20
IETA_CAT['HE2a'] = [21, 25]  ## Trigger towers 21 - 25
IETA_CAT['HE2b'] = [26, 28]  ## Trigger towers 26 - 28
IETA_CAT['HF']   = [30, 41]  ## Trigger towers 30 - 41

HCAL_ETA_CAT = ['HB', 'HE1', 'HE2', 'HF']

useAbsEtaBins = True
ETA_Bins = []
for iEta in range(-41,42):
    if iEta in [-29, 0, 29]:        continue;
    if useAbsEtaBins and iEta < 0:  continue;
    ETA_Bins.append(str(iEta))
#ETA_Bins.append('all')
ETA_Bins.append('HBEF')

map_iEta_Eta = {
    1:  [0.000, 0.087],
    2:  [0.087, 0.174],
    3:  [0.174, 0.261],
    4:  [0.261, 0.348],
    5:  [0.348, 0.435],
    6:  [0.435, 0.522],
    7:  [0.522, 0.609],
    8:  [0.609, 0.695],
    9:  [0.695, 0.783],
    10: [0.783, 0.870],

    11: [0.870, 0.957],
    12: [0.957, 1.044],
    13: [1.044, 1.131],
    14: [1.131, 1.218],
    15: [1.218, 1.305],
    16: [1.305, 1.392],
    17: [1.392, 1.479],
    18: [1.479, 1.566], # endcap starts
    19: [1.566, 1.653],    
    20: [1.653, 1.740],
    
    21: [1.740, 1.830],
    22: [1.830, 1.930],
    23: [1.930, 2.043],
    24: [2.043, 2.172],
    25: [2.172, 2.322],
    26: [2.322, 2.500],
    27: [2.500, 2.650],
    28: [2.650, 3.000], # endcap ends   
}

map_Jet_TowerIEta_IEta = {
    #0:  X,
    1:  1,
    2:  3,
    3:  5,
    4:  7,
    5:  9,
    6: 11,
    7: 13,
    8: 15,
    9: 17,
    10: 19,
    11: 21,
    12: 23,
    13: 25,
    14: 27,
    15: 29,
    16: 31,
    17: 33,
    18: 35,
    19: 37,
    20: 39,
    21: 41,
    22: 43,
    23: 46,
    24: 48,
    25: 52,
    26: 55,
    27: 59,
    28: 63,
    #29: XX,
    30: 69,
    31: 74,
    32: 78,
    33: 82,
    34: 86,
    35: 90,
    36: 94,
    37: 98,
    38: 102,
    39: 106,
    40: 110,
    41: 116
}

def convert_jetIPhi_to_jetTowerIPhi(jetIPhi):
    jetTowerIPhi = (jetIPhi + 1) / 2
    return jetTowerIPhi


jetPtBins_forCalibration = [3, 6, 9, 12, 15, 20, 25, 30, 35, 40, 45, 55, 70, 90, 120, 160, 200, 999.0]


def dIPhi(iPhiA, iPhiB):
    if iPhiA > iPhiB:
        if abs((iPhiA - 72) - iPhiB) < abs(iPhiA - iPhiB):
            return (iPhiA - 72) - iPhiB
        else:
            return iPhiA - iPhiB
    else:
        if abs(iPhiA - (iPhiB - 72)) < abs(iPhiA - iPhiB):
            return iPhiA - (iPhiB - 72)
        else:
            return iPhiA - iPhiB

## Delta-eta calculation taking into account there is no iEta = 0 or +/-29
def dIEta(iEtaA, iEtaB):
    dIEtaAbs = abs(iEtaA - iEtaB)
    if (iEtaA < -29) != (iEtaB < -29):
        dIEtaAbs -= 1
    if (iEtaA <   0) != (iEtaB <   0):
        dIEtaAbs -= 1
    if (iEtaA <  29) != (iEtaB <  29):
        dIEtaAbs -= 1
    if (iEtaA < iEtaB):
        dIEtaAbs *= -1
    return dIEtaAbs


def calculateJetIEta(eta):
    jetIEta_offlineJet = -50.0 # None # abs(vOff.Eta())
    for iEta, etaBinRange in map_iEta_Eta.items():
        if abs(eta) >= etaBinRange[0] and abs(eta) < etaBinRange[1]:
            jetIEta_offlineJet = float( iEta * math.copysign(1, eta) )
    
    return jetIEta_offlineJet

def convert_jetIEta_to_jetTowerIEta(jetIEta):
    jetTowerIEta = None
    jetIEtaAbs = abs(jetIEta)
    for jetTowerIEta_1, jetIEta_1 in map_Jet_TowerIEta_IEta.items():
        if jetIEtaAbs == jetIEta_1:
            jetTowerIEta = int( jetTowerIEta_1 * math.copysign(1, jetIEta) )
            break
    
    if not jetTowerIEta:
        print("convert_jetIEta_to_jetTowerIEta(jetIEta):: jetTowerIEta not set... why ?? \t\t **** ERROR ****")
        exit(0)
    
    return jetTowerIEta

def getJetPtCategory(pt):
    iPFJetPtCat = 'None'
    for iCat in PT_CAT.keys():
        if pt >= PT_CAT[iCat][0] and pt < PT_CAT[iCat][2]:
            iPFJetPtCat = iCat
            
    #if iPFJetPtCat == 'None':
    #    print "getJetPtCategory():: pt {}, iPFJetPtCat {}, PT_CAT {}".format(pt, iPFJetPtCat, PT_CAT)
    return iPFJetPtCat


def passGoldenJSON(goldenJSON, run, lumisection):
    run=str(run)
    passSelection = False
    if run not in goldenJSON.keys():
        if PrintLevel >= 1:
            print(f"passGoldenJSON():: run {run} not in GoldenJSON")
        return passSelection

    for LS_range_GJ in goldenJSON[run]:
        if lumisection >= LS_range_GJ[0] and lumisection <= LS_range_GJ[-1]:
            passSelection = True
            break

    return passSelection 


class BranchCollection(dict):
    def __init__(self, branches, prefix, count_key):
        super().__init__()
        for full_name, array in branches.items():
            if full_name == count_key:
                self["nObjects"] = array  # Store count separately
            else:
                attr_name = full_name.removeprefix(prefix)
                self[attr_name] = array

# Read ROOT file
def read_root_file(file_path):
    file = uproot.open(file_path)  # Open without 'with' to keep it open
    tree = file["Events"]  # Load Events tree
    return tree

# Extract branches efficiently
def extract_branches(tree):
    branch_dict = defaultdict(dict)

    for name in tree.keys():
        if name.startswith("Muon_") or name == "nMuon":
            branch_dict["muon"][name] = tree[name].array()
        elif name.startswith("Electron_") or name == "nElectron":
            branch_dict["electron"][name] = tree[name].array()
        elif name.startswith("L1UnpackedCaloTower_") or name == "nL1UnpackedCaloTower":
            branch_dict["uTT"][name] = tree[name].array()
        elif name.startswith("HcalUnpackedTPs_") or name == "nHcalUnpackedTPs":
            branch_dict["uTP_H"][name] = tree[name].array()
        elif name.startswith("EcalUnpackedTPs_") or name == "nEcalUnpackedTPs":
            branch_dict["uTP_E"][name] = tree[name].array()
        elif name.startswith("L1EmulCaloTower_") or name == "nL1EmulCaloTower":
            branch_dict["eTT"][name] = tree[name].array()
        elif name.startswith("L1EmulCaloCluster_") or name == "nL1EmulCaloCluster":
            branch_dict["eTC"][name] = tree[name].array()
        elif name.startswith("HcalEmulTPs_") or name == "nHcalEmulTPs":
            branch_dict["eTP_H"][name] = tree[name].array()
        elif name.startswith("EcalEmulTPs_") or name == "nEcalEmulTPs":
            branch_dict["eTP_E"][name] = tree[name].array()
        elif name.startswith("Jet_") or name == "nJet":
            branch_dict["jet"][name] = tree[name].array()
        elif name.startswith("GenJet_") or name == "nGenJet":
            branch_dict["Genjet"][name] = tree[name].array()
        elif name.startswith("L1EmulJet_") or name == "nL1EmulJet":
            branch_dict["emu"][name] = tree[name].array()
        elif name.startswith("L1Jet_") or name == "nL1Jet":
            branch_dict["unp"][name] = tree[name].array()
        elif name.startswith("PV_") or name == "PV_npvsGood":
            branch_dict["vtx"][name] = tree[name].array()
        elif name in ["event", "run", "luminosityBlock"]:
            branch_dict["evt"][name] = tree[name].array()
        # elif name.startswith("HLT_"):
        #     branch_dict["hlt"][name] = tree[name].array()

    return branch_dict



def run():
    start_time = time.time()

    # argument 1: OOT PU subtraction scheme name. For e.g. def, PFA2, PFA1p
    # argument 2: i/p L1 nano directory
    # argument 3: N.  Split total number of events into Nth part. For e.g. 10
    # argument 4: m. Run on mth quantile events out of spillted N events. For e.g. 0, 1, 2, 3 ..., 9


    parser = argparse.ArgumentParser()
    parser.add_argument('--l1nano',              type=str, dest='l1nanoPath', required=False, help="L1T nanos", default="../sample_root_files/nano_MC.root")
    # parser.add_argument('--l1nano',              type=str, dest='l1nanoPath', required=True, help="L1T nanos")
    parser.add_argument('--sampleName',            type=str, dest='sampleName'  , help="sampleName for o/p file name", default='JETMET')
    parser.add_argument('--HcalPUS',               type=str, dest='OOT_PU_scheme', help="HCAL OOT PUS scheme", default='PFA1p')
    parser.add_argument('--PUrangeTag',            type=str, dest='PUrangeTag', help="PU range tag", default='None')
    parser.add_argument('--N_parts',               type=int, dest='N_parts', help="Split i/p l1nanos into N_parts", default='1')
    parser.add_argument('--M_quantilesIpFilesSet', type=int, dest='M_quantilesIpFilesSet', help="Quantile of i/p l1nanos split", default='0')
    parser.add_argument('--outputname',              type=str, dest='output', required=False, help="output name", default="Nano_out")

    parseGroup1 = parser.add_mutually_exclusive_group(required=False)
    parseGroup1.add_argument('--l1MatchOffline', action='store_true',default= True)
    parseGroup1.add_argument('--l1MatchGen', action='store_true')
    parseGroup2 = parser.add_mutually_exclusive_group(required=False)
    parseGroup2.add_argument('--l1nanoChunkyDonut', action='store_true')
    parseGroup2.add_argument('--l1nanoPhiRing', action='store_true', default=True)
    parseGroup3 = parser.add_mutually_exclusive_group(required=False)
    parseGroup3.add_argument('--offlineCHSJet', action='store_true')
    parseGroup3.add_argument('--offlinePUPPIJet', action='store_true')

    args = parser.parse_args()
    print("args: {}".format(args))
    file_name          = args.l1nanoPath
    OOT_PU_scheme         = args.OOT_PU_scheme
    PUrangeTag            = args.PUrangeTag
    N_parts               = args.N_parts
    M_quantilesIpFilesSet = args.M_quantilesIpFilesSet
    l1MatchOffline        = args.l1MatchOffline
    l1MatchGen            = args.l1MatchGen
    l1nanoChunkyDonut   = args.l1nanoChunkyDonut
    l1nanoPhiRing       = args.l1nanoPhiRing
    sampleName            = args.sampleName
    offlineCHSJet         = args.offlineCHSJet
    offlinePUPPIJet       = args.offlinePUPPIJet
    output                = args.output


    print("Inputs: \n\t file_name/l1nanoPath: {}, \n\t sampleName: {}, \n\t  OOT_PU_scheme: {}, \n\t PUrangeTag: {}, \n\t N_parts: {}, \n\t M_quantilesIpFilesSet: {}, \n\t l1MatchOffline: {}, \n\t l1MatchGen: {}, \n\t l1nanoChunkyDonut: {}, \n\t l1nanoPhiRing: {})".format(
        file_name, sampleName, OOT_PU_scheme, PUrangeTag, N_parts, M_quantilesIpFilesSet, l1MatchOffline, l1MatchGen, l1nanoChunkyDonut, l1nanoPhiRing
    ))
    print(f"offlineCHSJet: {offlineCHSJet},  offlinePUPPIJet: {offlinePUPPIJet} ")

    tree = read_root_file(file_name)
    branch_data = extract_branches(tree)

    # Create collections dynamically
    Muon_br = BranchCollection(branch_data["muon"], "Muon_", "nMuon")
    Ele_br = BranchCollection(branch_data["electron"], "Electron_", "nElectron")
    uTT_br = BranchCollection(branch_data["uTT"], "L1UnpackedCaloTower_", "nL1UnpackedCaloTower")
    uTP_H_br = BranchCollection(branch_data["uTP_H"], "HcalUnpackedTPs_", "nHcalUnpackedTPs")
    uTP_E_br = BranchCollection(branch_data["uTP_E"], "EcalUnpackedTPs_", "nEcalUnpackedTPs")
    eTT_br = BranchCollection(branch_data["eTT"], "L1EmulCaloTower_", "nL1EmulCaloTower")
    eTC_br = BranchCollection(branch_data["eTC"], "L1EmulCaloCluster_", "nL1EmulCaloCluster")
    eTP_H_br = BranchCollection(branch_data["eTP_H"], "HcalEmulTPs_", "nHcalEmulTPs")
    eTP_E_br = BranchCollection(branch_data["eTP_E"], "EcalEmulTPs_", "nEcalEmulTPs")
    Jet_br = BranchCollection(branch_data["jet"], "Jet_", "nJet")
    GenJet_br = BranchCollection(branch_data["Genjet"], "GenJet_", "nGenJet")
    Emu_br = BranchCollection(branch_data["emu"], "L1EmulJet_", "nL1EmulJet")
    Unp_br = BranchCollection(branch_data["unp"], "L1Jet_", "nL1Jet")
    Vtx_br = BranchCollection(branch_data["vtx"], "PV_", "PV_npvsGood")
    Evt_br = BranchCollection(branch_data["evt"], "", "")
    # HLT_br = BranchCollection(branch_data["hlt"], "HLT_","")


    
    
    # in_file_names = [ l1nanoPath ]
    in_file_names = "N"
    sL1nano = "l1nanoChunkyDonut" if l1nanoChunkyDonut else "l1nanoPhiRing"
    out_file_str  = "L1Nano_%s_%s_%s_%s_%s.root" % (output, sampleName, sL1nano, OOT_PU_scheme, PUrangeTag)
    if runMode in ['CalCalibSF']:
        out_file_str = out_file_str.replace(".root", "_CalCalibSF.root")
    if runMode in ['CalibJetByHand']:
        #out_file_str = out_file_str.replace(".root", "_wCalibJetByHand.root")
        if isinstance(sipFileCalibSF, str) and ".root" in sipFileCalibSF: # read Layer2CalibSF from root files Siddhesh generated
            out_file_str = out_file_str.replace(".root", "_wCalibJetByHand.root")
        if isinstance(sipFileCalibSF, dict): #Read L1Jet CalibLayer2 SFs from csv files provided by Syed
            out_file_str = out_file_str.replace(".root", "_wCalibJetByHand_Calib%s.root" % (calibSFLable))
            out_file_str = out_file_str.replace(" ", "")
            out_file_str = out_file_str.replace("(", "_")
            out_file_str = out_file_str.replace(")", "_")
            out_file_str = out_file_str.replace("/", "DividedBy")            
    if sOutExt:
        out_file_str = out_file_str.replace(".root", "%s.root" % (sOutExt))
    if N_parts:
        #out_file_str  = "L1T_JetMET_Res_%s_part%d_of_%d.root" % (OOT_PU_scheme, M_quantilesIpFilesSet,N_parts)
        out_file_str = out_file_str.replace(".root", "_part%d_of_%d.root" % (M_quantilesIpFilesSet,N_parts))
    out_file_str = out_file_str.replace("__", "_")
    print("out_file_str: {}".format(out_file_str))

            

    usePUReweighting = False
    isMC = True
    for in_file in in_file_names:
        if "mc" in in_file.lower():
            isMC = True


    if usePUReweighting_0:
        '''
        if not isMC:
            print "input nanos are not for MC.... and still running with usePUReweighting == true ??? \t\t **** ERROR ****"
            exit(0)
        '''
        if isMC: usePUReweighting = True

    print("OOT_PU_scheme: {}".format(OOT_PU_scheme))
    '''
    if OOT_PU_scheme.lower() in ['def', 'pfa2']:            
        sipFileCalibSF_toUse = sipFileCalibSF.replace('$OOTPUS', 'PFA2')
    elif OOT_PU_scheme.lower() in ['pfa1p']:
        sipFileCalibSF_toUse = sipFileCalibSF.replace('$OOTPUS', 'PFA1p')
    '''
    sipFileCalibSF_toUse = sipFileCalibSF


    fOut_MLInputs = None
    fOut_MLInputs_writer = None
    isFirstEntry_WriteInputForML = True
    if runMode in ['makeInputForML']:        
        out_file_str_0 = out_file_str.replace('.root','.csv')
        print("runMode=makeInputForML:: output file: {}".format(out_file_str_0))
        fOut_MLInputs = open(out_file_str_0, mode='w')

        '''
        fieldnames = ['emp_name', 'dept', 'birth_month']
        writer = csv.DictWriter(fOut_MLInputs, fieldnames=fieldnames)

        writer.writeheader()
        writer.writerow({'emp_name': 'John Smith', 'dept': 'Accounting', 'birth_month': 'November'})
        writer.writerow({'emp_name': 'Erica Meyers', 'dept': 'IT', 'birth_month': 'March'})

        fOut_MLInputs.close()
        '''

    if runMode in ['makePUHisto']:
        print("Running with runMode = makePUHisto \n")
        
    # GoldenJSON -----------------------------------------------------------------------------
    goldenJSON = None
    #if type(GoldenJSONForData_list)
    # if isinstance(GoldenJSONForData_list, list) and len(GoldenJSONForData_list) > 0 and (not isMC):
    #     for GoldenJSONForData_ in GoldenJSONForData_list:
    #         #with open(GoldenJSONForData_) as fGoldenJSON_:
    #         #    goldenJSON_tmp = json.load(fGoldenJSON_)

    #         fGoldenJSON_ = None
    #         if GoldenJSONForData_.startswith('https://'):
    #             fGoldenJSON_ = urllib.request.urlopen(GoldenJSONForData_)
    #         else:
    #             fGoldenJSON_ = open(GoldenJSONForData_)

    #         goldenJSON_tmp = json.load(fGoldenJSON_)
    #         if not goldenJSON: goldenJSON = goldenJSON_tmp
    #         else: goldenJSON |= goldenJSON_tmp # append dict
    #         fGoldenJSON_.close()

    #     print(f"GoldenJSONForData_list: {GoldenJSONForData_list}")
    #     if PrintLevel >= 0:
    #         print(f"goldenJSON: {goldenJSON}")

    ###################
    ### Book histograms
    ###################
    if PrintLevel >= 1:
        print("Booking hitograms")
    hDummy = R.TH1D("hDummy","",1,0,1)
    # hDummy.SetDefaultSumw2()

    pt_bins  = [20,    0, 200]
    ptHigh_bins  = [200,    0, 1000]
    ptMid_bins   = [200,    0, 200]
    res_bins = [80, -1.5, 2.5]
    dR_bins  = [35,  0.0, 0.35]
    puEt_bins = [60,   0, 120]
    PUByRawPt_bins = [80, 0, 2]
    iEta_bins = [83, -41.5, 41.5]
    if useAbsEtaBins:
        iEta_bins = [41, 0.5, 41.5]
    CAL_TP_et_bins     = [ 513, 0,  513] # [1025, 0, 1025]
    CAL_TP_compEt_bins = [1025, 0, 1025]



    hnTotalEvents = R.TH1D("h_nTotalEvents", ";Channels;Events", 11,-0.5,10.5)

    hStat = R.TH1D("h_Stat", ";Conditions;Events", 51,-0.5,50.5)


    hnVts_vs_nTT_unp = R.TH2D("hnVts_vs_nTT_unp", "hnVts_vs_nTT_unp", 101,-0.5,100.5, 101,-0.5,2000.5)
    hnVts_vs_nTT_emu = R.TH2D("hnVts_vs_nTT_emu", "hnVts_vs_nTT_emu", 101,-0.5,100.5, 101,-0.5,2000.5)

    hdR_OffJet_OffMu_min  = R.TH1D("hdR_OffJet_OffMu_min",  "hdR_OffJet_OffMu_min",  60,0, 3.0)
    hdR_OffJet_OffEle_min = R.TH1D("hdR_OffJet_OffEle_min", "hdR_OffJet_OffEle_min", 60,0, 3.0)

    hdR_OffJet_OffMu_min_vs_vOffMuPtByvOffJetPt    = R.TH2D("hdR_OffJet_OffMu_min_vs_vOffMuPtByvOffJetPt",    "hdR_OffJet_OffMu_min_vs_vOffMuPtByvOffJetPt",    60,0, 3.0, 24, 0, 1.2)
    hdR_OffJet_OffEle_min_vs_vOffElePtByvOffJetPt  = R.TH2D("hdR_OffJet_OffEle_min_vs_vOffElePtByvOffJetPt",  "hdR_OffJet_OffEle_min_vs_vOffElePtByvOffJetPt",  60,0, 3.0, 24, 0, 1.2)

    hdR_OffJet_OffMu_min_forPtFracLt0p5  = R.TH1D("hdR_OffJet_OffMu_min_forPtFracLt0p5",  "hdR_OffJet_OffMu_min_forPtFracLt0p5",  60,0, 3.0)
    hdR_OffJet_OffEle_min_forPtFracLt0p5 = R.TH1D("hdR_OffJet_OffEle_min_forPtFracLt0p5", "hdR_OffJet_OffEle_min_forPtFracLt0p5", 60,0, 3.0)

    hdR_OffJet_L1Jet = {}
    for src in ['unp','emu']:
        hdR_OffJet_L1Jet[src] = R.TH1D("hdR_OffJet_L1Jet_%s" % src,  "hdR_OffJet_L1Jet_%s" % src,  60,0, 3.0)



    ## Variable distributions
    dists = []
    if runMode == '':
        dists = ['jet_eff', 'jet_num', 'jet_den', 'jet_res', 'jet_dR']

    hist = {}
    ## Loop over all distributions
    for dist in dists:
        hist[dist] = {}
        ## Loop over L1T jet algorithms
        for algo in ['PUS','noPUS','Raw','RawPUS']:
            hist[dist][algo] = {}
            ## Loop over eta regions
            for iEta in ETA_CAT.keys():
                hist[dist][algo][iEta] = {}
                ## Loop over pT ranges
                for iPt in list(PT_CAT.keys())+['PtAllBins']:
                    hist[dist][algo][iEta][iPt] = {}
                    ## Loop over unpacked and emulated
                    for src in ['unp','emu']:
                        hist[dist][algo][iEta][iPt][src] = []
                        ## Separate histogram for each input file (different TP reconstructions)
                        # for iTP in range(len(in_file_names)):

                        ## Pick binning for each histogram
                        h_bins = []
                        if   dist == 'jet_eff': h_bins = []
                        elif dist == 'jet_num': h_bins = pt_bins
                        elif dist == 'jet_den': h_bins = pt_bins
                        elif dist == 'jet_res': h_bins = res_bins
                        elif dist == 'jet_dR':  h_bins = dR_bins
                        else: print('\nInvalid distribution %s - no binning found!!! \nQuitting.'); sys.exit()

                        if dist == 'jet_eff':
                            hist[dist][algo][iEta][iPt][src].append(1 )
                        else:
                            hist[dist][algo][iEta][iPt][src].append( R.TH1D( 'h_%s_%s_%s_%s_%d_%s' % (dist, algo, iEta, iPt, src),
                                                                                'L1T %s %s %s in %s for %s' % (src, algo, dist, iEta, iPt),
                                                                                h_bins[0], h_bins[1], h_bins[2] ) )

                        ## End loop: for iTP in range(len(in_file_names))
                    ## End loop: for src in ['unp','emu']
                ## End loop: for iPt in PT_CAT.keys()
            ## End loop: for iEta in ETA_CAT.keys()
        ## End loop: for algo in ['PUS','noPUS','Raw','RawPUS']
    ## End loop: for dist in dists


    dists1 = [] # ['l1jetEt_vs_RawEtMinusPU']    
    hist1 = {}
    for dist in dists1:
        hist1[dist] = {}
        for iEta in ETA_Bins:            
            hist1[dist][iEta] = {}
            for src in ['unp','emu']:
                hist1[dist][iEta][src] = []

                ## Separate histogram for each input file (different TP reconstructions)
                # for iTP in range(len(in_file_names)):
                hist1[dist][iEta][src] .append( R.TH2D( 'h_%s_%s_%d_%s' % (dist, iEta, src),
                                                                                'L1T %s %s in ieta %s' % (src, dist, iEta),
                                                                                100,0,200, 100,0,3) )



    dists2 = [
        #'jet_byHand_eff', 'jet_byHand_num', 'jet_byHand_den', #'jet_byHand_res_vs_iEta', 'jet_byHand_PU_vs_iEta', 'jet_byHand_PUByRawPt_vs_iEta',
        # for calibration with PtAllBins 
        #'jet_byHand_L1JetPt_vs_PFJetPt',
        #'jet_byHand_L1JetPt_vs_DefaultL1JetPt',
        'jet_byHand_res_vs_iEta_vs_nVtx',
        'jet_byHand_res_woLayer2Calib_vs_iEta_vs_nVtx', 
    ]
    #dists2 = []
    hist2 = {}
    #print "\nSetting hist2::\ndists2: {}".format(dists2)
    #for jetShape in ['Default'] + JetShapes:
    for jetShape in JetShapes:
        # JetShape = "" plots are with the first version of code for 9x9 jets
        jetShape1 = jetShape
        if jetShape == 'Default':  jetShape1 = ""
        else:                      jetShape1 = "_%s" % (jetShape)

        ## Loop over all distributions
        for dist_1 in dists2:
            dist = '%s%s' % (dist_1, jetShape1)
            #print "    dist: {}".format(dist)

            hist2[dist] = {}
            ## Loop over L1T jet algorithms
            for algo in PUSAlgosAll: # ['Raw', 'RawPUS', 'RawPUS_phiRingMin4', 'RawPUS_phiRingSide4', 'RawPUS_phiRingAdjacent']: # ['PUS','noPUS','Raw','RawPUS']:
                if (algo == 'L1JDefault') and (jetShape != 'Default'): continue
                if (not l1nanoChunkyDonut) and (jetShape == 'Default') and (algo == 'RawPUS'):            continue
                if (not l1nanoPhiRing)     and (jetShape == 'Default') and (algo == 'RawPUS_phiDefault'): continue

                if dist_1 in ['jet_byHand_L1JetPt_vs_DefaultL1JetPt'] and algo != 'RawPUS': continue
                hist2[dist][algo] = {}                

                #print "   algo : {}".format(algo)
                ## Loop over eta regions
                for iEta in ETA_Bins: 
                    hist2[dist][algo][iEta] = {}
                    #print "    iEta: {}".format(iEta)
                    ## Loop over pT ranges
                    for iPt in list(PT_CAT.keys()) + ['PtAllBins']:
                        # L1JetPt_vs_PFJetPt plot to be plot with PtAllBins
                        if dist_1 in ['jet_byHand_L1JetPt_vs_PFJetPt', 'jet_byHand_L1JetPt_vs_DefaultL1JetPt'] and iPt != 'PtAllBins': continue

                        hist2[dist][algo][iEta][iPt] = {}
                        #print "    iPt: {}".format(iPt)
                        ## Loop over unpacked and emulated
                        for src in ['unp','emu']: #['unp','emu']:
                            hist2[dist][algo][iEta][iPt][src] = []
                            #print "    src: {}".format(src)
                            ## Separate histogram for each input file (different TP reconstructions)
                            # for iTP in range(len(in_file_names)):
                            ## Pick binning for each histogram
                            h_bins = []
                            if   'eff' in dist: h_bins = []
                            elif 'num' in dist: h_bins = pt_bins
                            elif 'den' in dist: h_bins = pt_bins
                            elif 'res' in dist: h_bins = res_bins
                            elif 'dR'  in dist: h_bins = dR_bins
                            elif '_PU_'  in dist: h_bins = puEt_bins
                            elif '_PUByRawPt_'  in dist: h_bins = PUByRawPt_bins
                            elif 'jet_byHand_L1JetPt_vs_PFJetPt' in dist: h_bins = pt_bins
                            elif 'jet_byHand_L1JetPt_vs_DefaultL1JetPt' in dist: h_bins = res_bins
                            else: print('\nInvalid distribution %s - no binning found!!! \nQuitting.'); sys.exit()

                            #print "      dist {}, algo {}, iEta {}, iPt {}, src {}, iTP {}".format(dist, algo, iEta, iPt, src, iTP)

                            if   'eff' in dist:
                                hist2[dist][algo][iEta][iPt][src].append(1 )
                            #elif dist in ['jet_byHand_num', 'jet_byHand_den']:
                            elif 'jet_byHand_num' in dist or 'jet_byHand_den' in dist:
                                # TH1D
                                hist2[dist][algo][iEta][iPt][src].append( R.TH1D( 'h_%s_%s_%s_%s_%s' % (dist, algo, iEta, iPt, src),
                                                                                    'L1T %s %s %s in %s for %s' % (src, algo, dist, iEta, iPt),
                                                                                    h_bins[0], h_bins[1], h_bins[2] ) )
                            elif dist_1 in ['jet_byHand_L1JetPt_vs_PFJetPt']:
                                # TH1D,  for calibration with PtAllBins 
                                hist2[dist][algo][iEta][iPt][src].append( R.TH2D( 'h_%s_%s_%s_%s_%s_fixBinWidthPFJetPt' % (dist, algo, iEta, iPt, src),
                                                                                    'L1T %s %s %s in %s for %s' % (src, algo, dist, iEta, iPt),
                                                                                    len(jetPtBins_forCalibration)-1, array('d',jetPtBins_forCalibration),
                                                                                    h_bins[0], h_bins[1], h_bins[2] ) )
                            elif dist_1 in ['jet_byHand_L1JetPt_vs_DefaultL1JetPt'] and jetShape == 'Default':
                                hist2[dist][algo][iEta][iPt][src].append( R.TH2D( 'h_%s_%s_%s_%s_%s' % (dist, algo, iEta, iPt, src),
                                                                                    'L1T %s %s %s in %s for %s' % (src, algo, dist, iEta, iPt),
                                                                                    len(jetPtBins_forCalibration)-1, array('d',jetPtBins_forCalibration),
                                                                                    h_bins[0], h_bins[1], h_bins[2]) )
                            elif dist_1 in ['jet_byHand_res_vs_iEta', 'jet_byHand_PU_vs_iEta', 'jet_byHand_PUByRawPt_vs_iEta']:
                                # TH2D
                                if iEta == 'HBEF':
                                    hist2[dist][algo][iEta][iPt][src].append( R.TH2D( 'h_%s_%s_%s_%s_%s' % (dist, algo, iEta, iPt, src),
                                                                                        'L1T %s %s %s in %s for %s' % (src, algo, dist, iEta, iPt),
                                                                                        iEta_bins[0],iEta_bins[1],iEta_bins[2], h_bins[0], h_bins[1], h_bins[2] ) )
                                else:
                                    continue
                            elif dist_1 in ['jet_byHand_res_vs_iEta_vs_nVtx', 'jet_byHand_res_woLayer2Calib_vs_iEta_vs_nVtx']:
                                if iEta == 'HBEF':
                                    hist2[dist][algo][iEta][iPt][src].append( R.TH3D( 'h_%s_%s_%s_%s_%s' % (dist, algo, iEta, iPt, src),
                                                                                        'L1T %s %s %s in %s for %s' % (src, algo, dist, iEta, iPt),
                                                                                        iEta_bins[0],iEta_bins[1],iEta_bins[2], 101,-0.5,100.5, h_bins[0], h_bins[1], h_bins[2] ) )

    jetRate_bins = [400, 0., 400.0];   jetRate_binWidth = (jetRate_bins[2] - jetRate_bins[1]) / jetRate_bins[0]; 
    eff_bins = [500, 0, 500]
    PU_bins  = [101, -0.5, 100.5]
    sAxexName_JetRate = ";Threshold E_{T} (GeV);nPV;Rate (Hz)"
    sAxexName_Eff     = ";E_{T} (GeV);nPV;Events / bin"

    dists3 = [
        'jet_byHand_eff_num_vs_PU', 'jet_byHand_eff_den_vs_PU'
    ]
    dists4 = [
        'jet_byHand_rates_singleJet', 'jet_byHand_rates_doubleJet', 'jet_byHand_rates_trippleJet', 'jet_byHand_rates_quadJet'
    ]

    hist3 = {}
    #print "\nSetting hist3::\ndists3: {}".format(dists3)
    #for jetShape in ['Default'] + JetShapes:
    for jetShape in JetShapes + JetShapesType2:
        # JetShape = "" plots are with the first version of code for 9x9 jets
        jetShape1 = jetShape
        if jetShape == 'Default':  jetShape1 = ""
        else:                      jetShape1 = "_%s" % (jetShape)

        ## Loop over all distributions
        for dist_1 in dists3:
            dist = '%s%s' % (dist_1, jetShape1)
            #print "    dist: {}".format(dist)

            hist3[dist] = {}
            ## Loop over L1T jet algorithms
            for algo in PUSAlgosAll + PUSAlgosAllType2: # ['Raw', 'RawPUS', 'RawPUS_phiRingMin4', 'RawPUS_phiRingSide4', 'RawPUS_phiRingAdjacent']: # ['PUS','noPUS','Raw','RawPUS']:
                # read proper jetShape and PUSAlgo combination
                if (jetShape in JetShapes      and algo not in PUSAlgosAll) or \
                (jetShape in JetShapesType2 and algo not in PUSAlgosAllType2 ):
                    continue

                if (algo == 'L1JDefault') and (jetShape != 'Default'): continue

                hist3[dist][algo] = {}                

                #print "   algo : {}".format(algo)
                ## Loop over eta categories
                for iEta in ETA_CAT.keys():
                    hist3[dist][algo][iEta] = {}
                    #print "    iEta: {}".format(iEta)

                    ## Loop over pT ranges
                    for iPt in PT_CAT.keys():                        
                        hist3[dist][algo][iEta][iPt] = {}
                        #print "    iPt: {}".format(iPt)

                        ## Loop over unpacked and emulated
                        for src in ['unp','emu']: #['unp','emu']:
                            hist3[dist][algo][iEta][iPt][src] = []
                            #print "    src: {}".format(src)

                            ## Separate histogram for each input file (different TP reconstructions)
                            # for iTP in range(len(in_file_names)):
                            ## Pick binning for each histogram
                            h_bins = []
                            sAxexName = None
                            if   'eff_num' in dist or 'eff_den' in dist:
                                h_bins = eff_bins
                                sAxexName = sAxexName_Eff
                            else: print('\nInvalid distribution %s - no binning found!!! \nQuitting.'); sys.exit()

                            #print "      dist {}, algo {}, iEta {}, iPt {}, src {}, iTP {}".format(dist, algo, iEta, iPt, src, iTP)

                            if 'jet_byHand_eff_num_vs_PU' in dist or 'jet_byHand_eff_den_vs_PU' in dist:
                                # TH2D
                                hist3[dist][algo][iEta][iPt][src].append( R.TH2D( 'h_%s_%s_%s_TrgTrsh%d_%s' % (dist, algo, iEta, PT_CAT[iPt][1], src),
                                                                                    sAxexName,
                                                                                    int(h_bins[0]), h_bins[1], h_bins[2],
                                                                                    int(PU_bins[0]),PU_bins[1],PU_bins[2] ) )
                                '''
                                print('h_%s_%s_%s_TrgTrsh%d_%d_%s : jetShape %s, dist_1 %s, algo %s, iEta %s, iPt %s, src %s, iTP %d' % \
                                        (dist, algo, iEta, PT_CAT[iPt][1], iTP, src,\
                                        jetShape,dist_1,algo,iEta,iPt,src,iTP
                                        )); sys.stdout.flush();
                                '''


    hist4 = {}
    #print "\nSetting hist4::\ndists4: {}".format(dists4)
    #for jetShape in ['Default'] + JetShapes:
    for jetShape in JetShapes + JetShapesType2:
        # JetShape = "" plots are with the first version of code for 9x9 jets
        jetShape1 = jetShape
        if jetShape == 'Default':  jetShape1 = ""
        else:                      jetShape1 = "_%s" % (jetShape)

        ## Loop over all distributions
        for dist_1 in dists4:
            dist = '%s%s' % (dist_1, jetShape1)
            #print "    dist: {}".format(dist)

            hist4[dist] = {}
            ## Loop over L1T jet algorithms
            for algo in PUSAlgosAll + PUSAlgosAllType2: # ['Raw', 'RawPUS', 'RawPUS_phiRingMin4', 'RawPUS_phiRingSide4', 'RawPUS_phiRingAdjacent']: # ['PUS','noPUS','Raw','RawPUS']:
                # read proper jetShape and PUSAlgo combination
                if (jetShape in JetShapes      and algo not in PUSAlgosAll) or \
                (jetShape in JetShapesType2 and algo not in PUSAlgosAllType2 ):
                    continue

                if (algo == 'L1JDefault') and (jetShape != 'Default'): continue

                hist4[dist][algo] = {}                

                #print "   algo : {}".format(algo)
                ## Loop over eta categories
                for iEta in ETA_CAT.keys():
                    hist4[dist][algo][iEta] = {}
                    #print "    iEta: {}".format(iEta)

                    ## Loop over unpacked and emulated
                    for src in ['unp','emu']: #['unp','emu']:
                        hist4[dist][algo][iEta][src] = []
                        #print "    src: {}".format(src)

                        ## Separate histogram for each input file (different TP reconstructions)
                        # for iTP in range(len(in_file_names)):
                            ## Pick binning for each histogram
                        h_bins = []
                        sAxexName = None
                        if 'jet_byHand_rates_' in dist:
                            h_bins = jetRate_bins
                            sAxexName = sAxexName_JetRate
                        else: print('\nInvalid distribution %s - no binning found!!! \nQuitting.'); sys.exit()

                        #print "      dist {}, algo {}, iEta {}, src {}, iTP {}".format(dist, algo, iEta, src, iTP)

                        if 'jet_byHand_rates_' in dist:
                            # TH2D
                            hist4[dist][algo][iEta][src].append( R.TH2D( 'h_%s_%s_%s_%s' % (dist, algo, iEta, src),
                                                                            sAxexName,
                                                                            int(h_bins[0]), h_bins[1], h_bins[2],
                                                                            int(PU_bins[0]),PU_bins[1],PU_bins[2] ) )


    # JetShapesType2, PUSAlgosAllType2
    dists5 = [
        #'l1TauMatchingPFJet_res_vs_iEta_vs_nVtx',
        #'jet_byHand_res_vs_iEta_vs_nVtx_l1TauDefault',
        'jet_byHand_res_vs_iEta_vs_nVtx', 
    ]
    dists5 = []
    hist5 = {}
    for jetShape in JetShapesType2:
        jetShape1 = "_%s" % (jetShape)        

        for dist_1 in dists5:
            dist = '%s%s' % (dist_1, jetShape1)
            hist5[dist] = {}

            for algo in PUSAlgosAllType2: # ['Et', 'RawEt']
                hist5[dist][algo] = {}

                ## Loop over eta categories
                for iEta in ETA_CAT.keys():
                    hist5[dist][algo][iEta] = {}

                    ## Loop over pT ranges
                    for iPt in PT_CAT.keys()+['PtAllBins']:                        
                        hist5[dist][algo][iEta][iPt] = {}

                        ## Loop over unpacked and emulated
                        for src in ['unp','emu']: #['unp','emu']:
                            hist5[dist][algo][iEta][iPt][src] = []

                            ## Separate histogram for each input file (different TP reconstructions)
                            # for iTP in range(len(in_file_names)):
                            ## Pick binning for each histogram
                            h_bins = []
                            sAxexName = None
                            if 'res' in dist:
                                h_bins = res_bins
                            else: print('\nInvalid distribution %s - no binning found!!! \nQuitting.'); sys.exit()

                            if dist_1 in ['jet_byHand_res_vs_iEta_vs_nVtx']:
                                if iEta == 'HBEF':
                                    hist5[dist][algo][iEta][iPt][src].append( R.TH3D( 'h_%s_%s_%s_%s_%s' % (dist, algo, iEta, iPt, src),
                                                                                            'L1T %s %s %s in %s for %s' % (src, algo, dist, iEta, iPt),
                                                                                            iEta_bins[0],iEta_bins[1],iEta_bins[2], 101,-0.5,100.5, h_bins[0], h_bins[1], h_bins[2] ) )



    #Calibrate CAL TPs:
    dists6 = [
        'ECAL_TP_et_vs_iEta_vs_nVts',  'ECAL_TP_compEt_vs_iEta_vs_nVts',
        'HCAL_TP_et_vs_iEta_vs_nVts',  'HCAL_TP_compEt_vs_iEta_vs_nVts',
    ]
    dists6 = []
    hist6 = {}
    for dist in dists6:
        hist6[dist] = {}

        for src in ['unp','emu']:
            hist6[dist][src] = []

            # for iTP in range(len(in_file_names)):
            h_bins = []
            if '_TP_et_' in dist:
                h_bins = CAL_TP_et_bins
            elif '_TP_compEt_' in dist:
                h_bins = CAL_TP_compEt_bins

            hist6[dist][src].append( R.TH2D( 'h_%s_%s' % (dist, src),
                                                'L1T %s in %s' % (dist, src),
                                                iEta_bins[0],iEta_bins[1],iEta_bins[2], 101,-0.5,100.5, h_bins[0], h_bins[1], h_bins[2] )  )

    #print "hist5: {}".format(hist5)


    # dists7
    dists7 = [
        'CheckPUEt_iEta_vs_PUEtRatioWrtDefault_9x9_RawPUS_phiDefault'
    ]
    hists7 = {}
    for dist in dists7:
        hists7[dist] = {}

        ## Loop over unpacked and emulated
        for src in ['emu']: #['unp','emu']:
            hists7[dist][src] = []

            ## Separate histogram for each input file (different TP reconstructions)
            # for iTP in range(len(in_file_names)):
            hists7[dist][src].append( R.TH2D( 'h_%s_%s' % (dist, src),
                                                'h_%s_%s' % (dist, src),
                                                iEta_bins[0],iEta_bins[1],iEta_bins[2],
                                                100, 0, 2) )


    # dists8
    dists8 = [
        'L1JetPtRawPUS', 'TT_iet', 'ECALTP_et', 'ECALTP_compEt', 'HCALTP_et', 'HCALTP_compEt',
        'L1JetPtRawPUS_vs_RefJetPt',

        'TT_iet_nRefJetsEq1', 'ECALTP_et_nRefJetsEq1', 'ECALTP_compEt_nRefJetsEq1', 'HCALTP_et_nRefJetsEq1', 'HCALTP_compEt_nRefJetsEq1',
        'TT_iet_RefJetEtalt2p5', 'ECALTP_et_RefJetEtalt2p5', 'ECALTP_compEt_RefJetEtalt2p5', 'HCALTP_et_RefJetEtalt2p5', 'HCALTP_compEt_RefJetEtalt2p5',
        'TT_iet_RefJetEtaBtw2p5And3p1', 'ECALTP_et_RefJetEtaBtw2p5And3p1', 'ECALTP_compEt_RefJetEtaBtw2p5And3p1', 'HCALTP_et_RefJetEtaBtw2p5And3p1', 'HCALTP_compEt_RefJetEtaBtw2p5And3p1',        
    ]
    hist8 = {}
    for dist in dists8:
        hist8[dist] = {}

        ## Loop over unpacked and emulated
        for src in ['unp','emu']:
            hist8[dist][src] = {}

            for iEta in ETA_Bins:

                if   dist in [
                        'L1JetPtRawPUS', 'TT_iet',
                        'TT_iet_nRefJetsEq1',
                        'TT_iet_RefJetEtalt2p5',
                        'TT_iet_RefJetEtaBtw2p5And3p1',                        
                ]:
                    hist8[dist][src][iEta] = R.TH1D( 'h%s_iEta%s_%s' % (dist, iEta, src),
                                                    '%s iEta%s %s' % (dist, iEta, src),
                                                    ptHigh_bins[0], ptHigh_bins[1], ptHigh_bins[2] )

                if   dist in [
                        'ECALTP_et', 'ECALTP_compEt', 'HCALTP_et', 'HCALTP_compEt',
                        'ECALTP_et_nRefJetsEq1', 'ECALTP_compEt_nRefJetsEq1', 'HCALTP_et_nRefJetsEq1', 'HCALTP_compEt_nRefJetsEq1',
                        'ECALTP_et_RefJetEtalt2p5', 'ECALTP_compEt_RefJetEtalt2p5', 'HCALTP_et_RefJetEtalt2p5', 'HCALTP_compEt_RefJetEtalt2p5',
                        'ECALTP_et_RefJetEtaBtw2p5And3p1', 'ECALTP_compEt_RefJetEtaBtw2p5And3p1', 'HCALTP_et_RefJetEtaBtw2p5And3p1', 'HCALTP_compEt_RefJetEtaBtw2p5And3p1',
                ]:
                    hist8[dist][src][iEta] = R.TH1D( 'h%s_iEta%s_%s' % (dist, iEta, src),
                                                    '%s iEta%s %s' % (dist, iEta, src),
                                                    ptMid_bins[0], ptMid_bins[1], ptMid_bins[2] )

                elif dist in ['L1JetPtRawPUS_vs_RefJetPt']:
                    hist8[dist][src][iEta] = R.TH2D( 'h%s_iEta%s_%s' % (dist, iEta, src),
                                                    '%s iEta%s %s' % (dist, iEta, src),
                                                    ptHigh_bins[0], ptHigh_bins[1], ptHigh_bins[2],
                                                    ptHigh_bins[0], ptHigh_bins[1], ptHigh_bins[2])



    hnL1JetUnp_0 = R.TH1D('hnL1JetUnp_0', 'hnL1JetUnp_0', 31, -0.5, 30.5)
    hL1JetUnp_Pt_0 = R.TH1D('hL1JetUnp_Pt_0', 'hL1JetUnp_Pt_0', 100, 0, 300)
    hL1JetUnp_Eta_0 = R.TH1D('hL1JetUnp_Eta_0', 'hL1JetUnp_Eta_0', 100, -6, 6)
    hL1JetUnp_Phi_0 = R.TH1D('hL1JetUnp_Phi_0', 'hL1JetUnp_Phi_0', 100, -3.14, 3.14)

    hnL1JetEmu_0 = R.TH1D('hnL1JetEmu_0', 'hnL1JetEmu_0', 31, -0.5, 30.5)
    hL1JetEmu_Pt_0 = R.TH1D('hL1JetEmu_Pt_0', 'hL1JetEmu_Pt_0', 100, 0, 300)
    hL1JetEmu_Eta_0 = R.TH1D('hL1JetEmu_Eta_0', 'hL1JetEmu_Eta_0', 100, -6, 6)
    hL1JetEmu_Phi_0 = R.TH1D('hL1JetEmu_Phi_0', 'hL1JetEmu_Phi_0', 100, -3.14, 3.14)

    hnOfflineJet_0 = R.TH1D('hnhOfflineJet_0', 'hnhOfflineJet_0', 31, -0.5, 30.5)
    hOfflineJet_Pt_0 = R.TH1D('hhOfflineJet_Pt_0', 'hhOfflineJet_Pt_0', 100, 0, 300)
    hOfflineJet_Eta_0 = R.TH1D('hhOfflineJet_Eta_0', 'hhOfflineJet_Eta_0', 100, -6, 6)
    hOfflineJet_Phi_0 = R.TH1D('hhOfflineJet_Phi_0', 'hhOfflineJet_Phi_0', 100, -3.14, 3.14)

    #
    hnL1JetUnp_1 = R.TH1D('hnL1JetUnp_1', 'hnL1JetUnp_1', 31, -0.5, 30.5)
    hL1JetUnp_Pt_1 = R.TH1D('hL1JetUnp_Pt_1', 'hL1JetUnp_Pt_1', 100, 0, 300)
    hL1JetUnp_Eta_1 = R.TH1D('hL1JetUnp_Eta_1', 'hL1JetUnp_Eta_1', 100, -6, 6)
    hL1JetUnp_Phi_1 = R.TH1D('hL1JetUnp_Phi_1', 'hL1JetUnp_Phi_1', 100, -3.14, 3.14)

    hnL1JetEmu_1 = R.TH1D('hnL1JetEmu_1', 'hnL1JetEmu_1', 31, -0.5, 30.5)
    hL1JetEmu_Pt_1 = R.TH1D('hL1JetEmu_Pt_1', 'hL1JetEmu_Pt_1', 100, 0, 300)
    hL1JetEmu_Eta_1 = R.TH1D('hL1JetEmu_Eta_1', 'hL1JetEmu_Eta_1', 100, -6, 6)
    hL1JetEmu_Phi_1 = R.TH1D('hL1JetEmu_Phi_1', 'hL1JetEmu_Phi_1', 100, -3.14, 3.14)

    #
    hnL1JetUnp_2 = R.TH1D('hnL1JetUnp_2', 'hnL1JetUnp_2', 31, -0.5, 30.5)
    hL1JetUnp_Pt_2 = R.TH1D('hL1JetUnp_Pt_2', 'hL1JetUnp_Pt_2', 100, 0, 300)
    hL1JetUnp_Eta_2 = R.TH1D('hL1JetUnp_Eta_2', 'hL1JetUnp_Eta_2', 100, -6, 6)
    hL1JetUnp_Phi_2 = R.TH1D('hL1JetUnp_Phi_2', 'hL1JetUnp_Phi_2', 100, -3.14, 3.14)

    hnL1JetEmu_2 = R.TH1D('hnL1JetEmu_2', 'hnL1JetEmu_2', 31, -0.5, 30.5)
    hL1JetEmu_Pt_2 = R.TH1D('hL1JetEmu_Pt_2', 'hL1JetEmu_Pt_2', 100, 0, 300)
    hL1JetEmu_Eta_2 = R.TH1D('hL1JetEmu_Eta_2', 'hL1JetEmu_Eta_2', 100, -6, 6)
    hL1JetEmu_Phi_2 = R.TH1D('hL1JetEmu_Phi_2', 'hL1JetEmu_Phi_2', 100, -3.14, 3.14)


    hCaloTTi_iEta_vs_Phi_tmplate = R.TH2D( 'hCaloTT_iEta_vs_iPhi', '', 83,-41.5,41.5, 72,0.5,72.5, ) # iEta on xaxis and iPhi on y-axis

    hCaloTowers_iEta_vs_iPhi_list = []
    hCaloTTs_iEta_vs_iPhi_list    = []


    hTTEtMax_forL1JetPUEt0 = None
    hL1JetRawEt_vs_L1JetEt_forL1JetPUEt0 = None
    if runMode in ['trbshtPhiRingPUS']: 
        hTTEtMax_forL1JetPUEt0 = R.TH1D('hTTEtMax_forL1JetPUEt0', '', 400, 0, 800)
        hL1JetRawEt_vs_L1JetEt_forL1JetPUEt0 = R.TH2D('hL1JetRawEt_vs_L1JetEt_forL1JetPUEt0', '', 400,0,1200, 400,0,1200)


    # Check JEC SFs
    hJEC_iEta_vs_Pt = R.TH2D('hJEC_iEta_vs_Pt', 'hJEC_iEta_vs_Pt', 41,0.5,41.5, 257*2,-0.25,256.75 )

    hRefJet_pt_0  = R.TH1D('hRefJet_pt_0',  'hRefJet_pt_0',  100, 0, 400 )
    hRefJet_eta_0 = R.TH1D('hRefJet_eta_0', 'hRefJet_eta_0', 100, -6, 6 )
    hRefJet_phi_0 = R.TH1D('hRefJet_phi_0', 'hRefJet_phi_0', 100, -3.14, 3.14 )

    hRefJet_pt_0_1  = R.TH1D('hRefJet_pt_0_1',  'hRefJet_pt_0_1',  100, 0, 400 )
    hRefJet_eta_0_1 = R.TH1D('hRefJet_eta_0_1', 'hRefJet_eta_0_1', 100, -6, 6 )
    hRefJet_phi_0_1 = R.TH1D('hRefJet_phi_0_1', 'hRefJet_phi_0_1', 100, -3.14, 3.14 )

    hRefJet_pt_0_1_test = {}
    hRefJet_eta_0_1_test = {}
    hRefJet_phi_0_1_test = {}
    for idx_ in range( nJetFilters + 1 ):
        hRefJet_pt_0_1_test[idx_]  = R.TH1D('hRefJet_pt_0_1_test_%s' % idx_,  'hRefJet_pt_0_1_test_%s' % idx_,  100, 0, 400 )
        hRefJet_eta_0_1_test[idx_] = R.TH1D('hRefJet_eta_0_1_test_%s' % idx_, 'hRefJet_eta_0_1_test_%s' % idx_, 100, -6, 6 )
        hRefJet_phi_0_1_test[idx_] = R.TH1D('hRefJet_phi_0_1_test_%s' % idx_, 'hRefJet_phi_0_1_test_%s' % idx_, 100, -3.14, 3.14 )

    hRefJet_pt_0_2  = R.TH1D('hRefJet_pt_0_2',  'hRefJet_pt_0_2',  100, 0, 400 )
    hRefJet_eta_0_2 = R.TH1D('hRefJet_eta_0_2', 'hRefJet_eta_0_2', 100, -6, 6 )
    hRefJet_phi_0_2 = R.TH1D('hRefJet_phi_0_2', 'hRefJet_phi_0_2', 100, -3.14, 3.14 )

    hRefJet_pt_0_3  = R.TH1D('hRefJet_pt_0_3',  'hRefJet_pt_0_3',  100, 0, 400 )
    hRefJet_eta_0_3 = R.TH1D('hRefJet_eta_0_3', 'hRefJet_eta_0_3', 100, -6, 6 )
    hRefJet_phi_0_3 = R.TH1D('hRefJet_phi_0_3', 'hRefJet_phi_0_3', 100, -3.14, 3.14 )

    hRefJet_pt_0_4  = R.TH1D('hRefJet_pt_0_4',  'hRefJet_pt_0_4',  100, 0, 400 )
    hRefJet_eta_0_4 = R.TH1D('hRefJet_eta_0_4', 'hRefJet_eta_0_4', 100, -6, 6 )
    hRefJet_phi_0_4 = R.TH1D('hRefJet_phi_0_4', 'hRefJet_phi_0_4', 100, -3.14, 3.14 )

    hRefJet_pt_1  = R.TH1D('hRefJet_pt_1',  'hRefJet_pt_1',  100, 0, 400 )
    hRefJet_eta_1 = R.TH1D('hRefJet_eta_1', 'hRefJet_eta_1', 100, -6, 6 )
    hRefJet_phi_1 = R.TH1D('hRefJet_phi_1', 'hRefJet_phi_1', 100, -3.14, 3.14 )

    hRefJet_pt_2  = R.TH1D('hRefJet_pt_2',  'hRefJet_pt_2',  100, 0, 400 )
    hRefJet_eta_2 = R.TH1D('hRefJet_eta_2', 'hRefJet_eta_2', 100, -6, 6 )
    hRefJet_phi_2 = R.TH1D('hRefJet_phi_2', 'hRefJet_phi_2', 100, -3.14, 3.14 )

    hRefJet_pt_2p01  = R.TH1D('hRefJet_pt_2p01',  'hRefJet_pt_2p01',  100, 0, 400 )
    hRefJet_eta_2p01 = R.TH1D('hRefJet_eta_2p01', 'hRefJet_eta_2p01', 100, -6, 6 )
    hRefJet_phi_2p01 = R.TH1D('hRefJet_phi_2p01', 'hRefJet_phi_2p01', 100, -3.14, 3.14 )

    hRefJet_pt_2p02  = R.TH1D('hRefJet_pt_2p02',  'hRefJet_pt_2p02',  100, 0, 400 )
    hRefJet_eta_2p02 = R.TH1D('hRefJet_eta_2p02', 'hRefJet_eta_2p02', 100, -6, 6 )
    hRefJet_phi_2p02 = R.TH1D('hRefJet_phi_2p02', 'hRefJet_phi_2p02', 100, -3.14, 3.14 )

    hRefJet_pt_2p03  = R.TH1D('hRefJet_pt_2p03',  'hRefJet_pt_2p03',  100, 0, 400 )
    hRefJet_eta_2p03 = R.TH1D('hRefJet_eta_2p03', 'hRefJet_eta_2p03', 100, -6, 6 )
    hRefJet_phi_2p03 = R.TH1D('hRefJet_phi_2p03', 'hRefJet_phi_2p03', 100, -3.14, 3.14 )


    if PrintLevel >= 1: 
        print("Histograms booked"); sys.stdout.flush();

    # nVtx
    hnVtx       = R.TH1D("nVtx",       "nVtx",       201,-0.5,200.5);
    hnVtx_ReWtd = R.TH1D("nVtx_ReWtd", "nVtx_ReWtd", 201,-0.5,200.5);
    hPUWt       = None
    hnVtxData   = None
    hnVtxMC     = None
    hnVtxMCWtd  = None

    if usePUReweighting:
        print("ipFilePUData: %s \nipFilePUMC: %s" % (ipFilePUData, ipFilePUMC))
        fInPUData = R.TFile(ipFilePUData)
        fInPUMC   = R. TFile(ipFilePUMC)
        if  ( not fInPUData.IsOpen() ):
            print("PUReweighting: input file %s couldn't open" % (ipFilePUData))
            exit(0)
        if ( not fInPUMC.IsOpen() ):
            print("PUReweighting: input file %s couldn't open" % (ipFilePUMC))
            exit(0)

        hnVtxData = fInPUData.Get(sHistoPUData);
        hnVtxMC   = fInPUMC.Get(sHistoPUMC);
        if ( not hnVtxData):
            print("PUReweighting: Couldn't fetch histogram %s  from %s  \t\t ******** ERROR *********\n" % (sHistoPUData, ipFilePUData))
            exit(0)

        if ( not hnVtxMC):
            print("PUReweighting: Couldn't fetch histogram %s  from %s  \t\t ******** ERROR *********\n" % (sHistoPUMC, ipFilePUMC))
            exit(0)

        if (hnVtxData.GetNbinsX() != hnVtxMC.GetNbinsX()): 
            print("PUReweighting: hnVtxData->GetNbinsX() != hnVtxMC->GetNbinsX() \t\t ******** ERROR *********\n");
            exit(0)

        hnVtxData.SetName("%s_Data" % (hnVtxData.GetName()));
        hnVtxMC.SetName("%s_MC" % (hnVtxMC.GetName()));

        # scale hnVtxData, hnVtxMC to have area = 1
        print("Area before scaling: hnVtxData %g, \t hnVtxMC %g " % (hnVtxData.Integral(1,hnVtxData.GetNbinsX()), hnVtxMC.Integral(1,hnVtxMC.GetNbinsX())))
        hnVtxData.Scale(1. / hnVtxData.Integral(1,hnVtxData.GetNbinsX()))
        hnVtxMC.Scale(1. / hnVtxMC.Integral(1,hnVtxMC.GetNbinsX()))    
        print("Area after scaling: hnVtxData %g, \t hnVtxMC %g " % (hnVtxData.Integral(1,hnVtxData.GetNbinsX()), hnVtxMC.Integral(1,hnVtxMC.GetNbinsX())))

        # calculate PU weights
        hPUWt = R.TH1D("hPUWeight","PU weight", hnVtxData.GetNbinsX(), hnVtxData.GetXaxis().GetXmin(), hnVtxData.GetXaxis().GetXmax())
        print("PUWeights::\n")
        for iBin in range(1, hPUWt.GetNbinsX()+1):
            nVtx = hPUWt.GetBinCenter(iBin);
            nCountsData = hnVtxData.GetBinContent(iBin);
            nCountsMC   = hnVtxMC.GetBinContent(iBin);
            PUwt = 0; 

            if nCountsData > 0 and nCountsMC > 0:
                PUwt = nCountsData / nCountsMC

            hPUWt.SetBinContent(iBin, PUwt);
            print("iBin %3d, nVtx %3g, PUwt %g, \t nCountsData %g,  nCountsMC %g  " % (iBin,nVtx,PUwt, nCountsData,nCountsMC))

        hnVtxMCWtd = hnVtxMC.Clone("%s_wtd" % (hnVtxMC.GetName()))
        for iBin in range(1, hnVtxMCWtd.GetNbinsX()+1):
            wt = hPUWt.GetBinContent(iBin);
            hnVtxMCWtd.SetBinContent(iBin, hnVtxMCWtd.GetBinContent(iBin) * wt);




    ### Read Layer2 CalibSF -----------------------------------------------------------------
    if PrintLevel >= 1:
        print("Read PUwgts. Now read Layer2CalibSF"); sys.stdout.flush();

    calibSFs = OrderedDict()
    if runMode in ['CalibJetByHand'] and isinstance(sipFileCalibSF, str) and ".root" in sipFileCalibSF:
        print("Now read Layer2CalibSF from root files Siddhesh generated ")

        fInFileCalibSF = R.TFile(sipFileCalibSF_toUse)
        if not fInFileCalibSF.IsOpen():
            print("CalibJetByHand: input file %s couldn't open" % (sipFileCalibSF_toUse))
            exit(0)
        print("CalibJetByHand: input calibration SF reading from %s " % (sipFileCalibSF_toUse))

        if PrintLevel >= 1:
            print("CalibJetByHand:: Calibration SFs::")
        #for jetShape in ['Default'] + JetShapes:
        for jetShape in JetShapes + JetShapesType2:
            # JetShape = "" plots are with the first version of code for 9x9 jets
            jetShape1 = jetShape
            if jetShape == 'Default':  jetShape1 = ""
            else:                      jetShape1 = "_%s" % (jetShape)
            if PrintLevel >= 1:
                print("jetShape: %s" % (jetShape))


            #calibSFs[jetShape] = {}
            calibSFs_jetShape_tmp1 = OrderedDict()
            print("0: jetShape: {}, calibSFs_jetShape_tmp1: {}".format(jetShape, calibSFs_jetShape_tmp1))
            for PUSAlgo in PUSAlgosSelected + PUSAlgosAllType2:

                # read proper jetShape and PUSAlgo combination
                if (jetShape in JetShapes      and PUSAlgo not in PUSAlgosSelected ) or \
                (jetShape in JetShapesType2 and PUSAlgo not in PUSAlgosAllType2 ):
                    continue

                #calibSFs[jetShape][PUSAlgo] = {}                
                if PrintLevel >= 1:
                    print("%4s%s" % (" ", PUSAlgo))

                sHistoName = sHistoCalibSF
                sHistoName = sHistoName.replace('$JETSHAPE',     jetShape1)
                sHistoName = sHistoName.replace('$PUSALGORITHM', PUSAlgo)
                if OOT_PU_scheme.lower() in ['def', 'pfa2']:
                    sHistoName = sHistoName.replace('$OOTPUS', 'PFA2')
                elif OOT_PU_scheme.lower() in ['pfa1p']:
                    sHistoName = sHistoName.replace('$OOTPUS', 'PFA1p')
            # print "CalibJetByHand: input calibration SF histogram %s" % (sHistoName)

                hHisto_SF_vs_iEta = fInFileCalibSF.Get(sHistoName)
                if not hHisto_SF_vs_iEta: # L1T calib SF are measured only for a few jetShape and PUSAlgorithms
                    continue

                print("CalibJetByHand: input calibration SF histogram %s" % (sHistoName))
                calibSFs_jetShape_tmp1[PUSAlgo] = OrderedDict()

                for iEta in ETA_Bins:
                    if iEta == 'HBEF': continue
                    if PrintLevel >= 1:
                        print("%8s%3s" % (" ", iEta))

                    binNumber_iEta = hHisto_SF_vs_iEta.FindBin(float(iEta))
                    SF_tmp         = hHisto_SF_vs_iEta.GetBinContent(binNumber_iEta)
                    calibSFs_jetShape_tmp1[PUSAlgo][iEta] = []
                    calibSFs_jetShape_tmp1[PUSAlgo][iEta].append([0., 999., SF_tmp]) # ['L1T pT bin low-edge',  'L1T pT bin high-edge',  'L1T layer2 calib SF']

                    '''
                    sHistoName = sHistoCalibSF
                    sHistoName = sHistoName.replace('$JETSHAPE',     jetShape1)
                    sHistoName = sHistoName.replace('$PUSALGORITHM', PUSAlgo)
                    sHistoName = sHistoName.replace('$ETA',          iEta)
                    sHistoName = sHistoName.replace('$L1MODE',       'emu')
                    hProfile = fInFileCalibSF.Get(sHistoName)
                    if not hProfile:
                        print "CalibJetByHand: hProfile %s couldn't read" % (sHistoName)
                        exit(0)


                    calibSFs[jetShape][PUSAlgo][iEta] = []
                    for iBin in range(1, hProfile.GetNbinsX()+1):
                        binLowEdge = hProfile.GetXaxis().GetBinLowEdge(iBin)
                        binUpEdge  = hProfile.GetXaxis().GetBinUpEdge(iBin)
                        binCenter  = hProfile.GetXaxis().GetBinCenter(iBin)
                        AvgPFJetPt = hProfile.GetBinContent(iBin)
                        calibSF    = AvgPFJetPt / binCenter
                        if calibSF < 1e-3: calibSF = 1 # 
                        calibSFs[jetShape][PUSAlgo][iEta].append([binLowEdge, binUpEdge, calibSF])
                        if PrintLevel >= 1:
                            print "%12s[%g, %g]: SF = %g / %g = %g" % (" ", binLowEdge,binUpEdge, AvgPFJetPt,binCenter,calibSF)
                    '''

            print("1: jetShape: {}, calibSFs_jetShape_tmp1: {}".format(jetShape, calibSFs_jetShape_tmp1))
            if calibSFs_jetShape_tmp1: # SF for the corresponding jetShape exist
                calibSFs[jetShape] = calibSFs_jetShape_tmp1

        print("\n\nCalibJetByHand: calibSFs: {}".format(json.dumps(calibSFs, sort_keys=True, indent=4)))


    #if updateCalibSF_DefaultRaw and l1nanoChunkyDonut:
    #    if 'Default' in calibSFs_additionalCorr.keys() and 'Raw' in calibSFs_additionalCorr['Default'].keys():
    #        calibSFs_additionalCorr['Default']['Raw'] = 1


    print("calibSFs_additionalCorr {} ".format(calibSFs_additionalCorr)); sys.stdout.flush();



    ### Read L1Jet CalibLayer2 SFs from csv files provided by Syed -----------------------------------------------------------------
    if runMode in ['CalibJetByHand'] and isinstance(sipFileCalibSF, dict) :
        print("Now read Layer2CalibSF from csv files provided by Syed")

        calibSF_L1JetPtRangeMin = calibSF_L1JetPtRange[0]
        calibSF_L1JetPtRangeMax = calibSF_L1JetPtRange[1]
        calibSF_L1JetPtBinWidth = calibSF_L1JetPtRange[2]

        for jetShape in sipFileCalibSF.keys():
            #print("sipFileCalibSF: jetShape {}".format(jetShape))
            if jetShape not in JetShapes + JetShapesType2:
                print("sipFileCalibSF: jetShape {} not in JetShapes or JetShapesType2".format(jetShape))
                continue

            calibSFs_jetShape_tmp1 = OrderedDict()
            for PUSAlgo in sipFileCalibSF[jetShape].keys():
                #print("    PUSAlgo {}".format(PUSAlgo))
                if PUSAlgo not in PUSAlgosAll:
                    print("sipFileCalibSF: PUSAlgo {} not in PUSAlgosAll".format(PUSAlgo))
                    continue


                sipFileCalibSF_toRead = sipFileCalibSF[jetShape][PUSAlgo]['fileName']
                print("sipFileCalibSF[{}][{}]: {} ".format(jetShape, PUSAlgo, sipFileCalibSF_toRead))
                calibSFs_jetShape_tmp1[PUSAlgo] = OrderedDict()
                with open(sipFileCalibSF_toRead, mode='r') as fipFileCalibSF_toRead:
                    calibSF_csv_reader = csv.DictReader(fipFileCalibSF_toRead)

                    #print("\n\ncalibSF_csv_reader: {} {}\n\n\n".format(type(calibSF_csv_reader), calibSF_csv_reader))
                    for calibSF_csv_row in calibSF_csv_reader:
                        iEta_tmp = calibSF_csv_row['L1JetTowerIEtaAbs'] # str
                        #l1JetPt_tmp = float(calibSF_csv_row['PhiRingEnergy'])
                        l1JetPt_tmp = float(calibSF_csv_row[ sipFileCalibSF[jetShape][PUSAlgo]['L1JetPtVarName'] ])
                        SF_tmp = float(calibSF_csv_row[ sipFileCalibSF[jetShape][PUSAlgo]['SFLabel'] ])
                        if PrintLevel >= 11:
                            print("calibSF_csv_row: {}, iEta {} {}, l1JetPt {} {}, SF {} {}={}".format(
                                calibSF_csv_row, type(iEta_tmp),iEta_tmp, type(l1JetPt_tmp),l1JetPt_tmp, type(SF_tmp),SF_tmp,calibSF_csv_row[calibSFLable]))

                        if l1JetPt_tmp < calibSF_L1JetPtRangeMin or l1JetPt_tmp > calibSF_L1JetPtRangeMax: continue

                        pT_bin_min = l1JetPt_tmp - (calibSF_L1JetPtBinWidth / 2.)
                        pT_bin_max = l1JetPt_tmp + (calibSF_L1JetPtBinWidth / 2.)
                        if l1JetPt_tmp < calibSF_L1JetPtRangeMin + (calibSF_L1JetPtBinWidth / 2.):
                            pT_bin_min = 0
                        if l1JetPt_tmp > calibSF_L1JetPtRangeMax - (calibSF_L1JetPtBinWidth / 2.):
                            pT_bin_max = 9999.
                        if PrintLevel >= 11:
                            print("pT_bin_min {}, pT_bin_max {}".format(pT_bin_min, pT_bin_max))

                        if iEta_tmp not in calibSFs_jetShape_tmp1[PUSAlgo].keys():
                            calibSFs_jetShape_tmp1[PUSAlgo][iEta_tmp] = []
                        calibSFs_jetShape_tmp1[PUSAlgo][iEta_tmp].append([pT_bin_min, pT_bin_max, SF_tmp])

            if calibSFs_jetShape_tmp1:
                calibSFs[jetShape] = calibSFs_jetShape_tmp1

        print("\n\nCalibJetByHand: calibSFs: {}".format(json.dumps(calibSFs, sort_keys=True, indent=4)))



    # save selected event numbers into 'selectedEvents_list'    
    selectedEvents_list = []


    # run on selected events
    eventsToRun_list = []

    if sFInEventsToRun:
        with open(sFInEventsToRun) as fInEventsToRun:
            for line in fInEventsToRun.readlines():
                line1 = line.replace('\n', '')
                #print(" line: >>%s<<" % line1)
                eventsToRun_list.append( line1 )

        print("eventsToRun_list: {} \n\t {} ".format(len(eventsToRun_list), eventsToRun_list))
    #return


    if PrintLevel >= 0:
        print(f"Start event loop: isMC: {isMC}"); sys.stdout.flush();

    nTotalEvents_byChains = []


    Ntot=len(Jet_br["pt"])
    print(f"Total events: {Ntot}")  
    print(list(Jet_br.keys()))

    l1MatchGen= True
    l1MatchOffline=False
    nTotalEvents_byChains=[]
    nTotalEvents_byChains.append(0)
    for iEvent in range(Ntot):

        l1JetRef_br = None
        nRefJets    = 0
        et_RefJets  = None
        eta_RefJets = None
        phi_RefJets = None
        if l1MatchOffline:
            nRefJets= Jet_br['nObjects'][iEvent]
            nOffJets= Jet_br['nObjects'][iEvent]
            et_RefJets=Jet_br["pt"][iEvent]
            # print(et_RefJets)
            eta_RefJets=Jet_br["eta"][iEvent]
            phi_RefJets=Jet_br["phi"][iEvent]

        elif l1MatchGen:
            nRefJets= GenJet_br['nObjects'][iEvent]
            nOffJets= GenJet_br['nObjects'][iEvent]
            et_RefJets=GenJet_br["pt"][iEvent]
            # print(et_RefJets)
            eta_RefJets=GenJet_br["eta"][iEvent]
            phi_RefJets=GenJet_br["phi"][iEvent]

        # print("Event: %d, nRefJets: %d" % (iEvent, nRefJets))

        nTotalEvents_byChains[0] += 1
        hStat.Fill(0)

        #Analyze (GEN.nVtx == 0) events from SinglePhoton_EpsilonPU sample to trouble-shoot high SFs in iEta 28 ----
        if isMC and useCutGenNVtxEq0:
            if Vtx_br['nObjects'][iEvent] > 0: continue
        # ----------------------------------------------------------------------------------------------------------
            

    
        # if not isMC and len(GoldenJSONForData_list) > 0:
        #     if not passGoldenJSON(goldenJSON, int(Evt_br['run'][iEvent]), int(Evt_br['luminosityBlock'][iEvent])):
        #         # print(f"Run:LS:Event:  %d:%d:%d   fails GoldenJSON " %(int(Evt_br['run']), int(Evt_br['luminosityBlock']), int(Evt_br['event']))); sys.stdout.flush();
        #         continue


        dataEra = ''
        if not isMC:
            for Era_, eraRunRange_ in dataErasRunRange.items():
                if int(Evt_br['run'][iEvent]) >= eraRunRange_[0] and int(Evt_br['run'][iEvent]) <= eraRunRange_[1]:
                    dataEra = Era_
                    break

        hStat.Fill(1)
        # print(list(HLT_br.keys()))
        # Apply HLT triggers requirements -------------------------------------------------
        # if not isMC and len(HLT_Triggers_Required) > 0:
        #     passHLTTrgs = False
        #     # if PrintLevel >= 20:
        #     #     print(f"Evt_br.hlt.size(): {Evt_br.hlt.size()}"); sys.stdout.flush();
        #     for HLT_TRG_name_required in HLT_Triggers_Required:
        #         if any(HLT_TRG_name_required in name for name in list(HLT_br.keys())):
        #             # print("Passed")
        #             passHLTTrgs = True
        #             break

        #     if not passHLTTrgs: 
        #         continue

        hStat.Fill(2)

        
        hCaloTowers_iEta_vs_iPhi = None
        hCaloTTs_iEta_vs_iPhi    = None

        if sFInEventsToRun: # and runMode in ['trbshtPhiRingPUS']:
            sRunLSEvent = "%s:%s:%s" % (int(Evt_br['run'][iEvent]), int(Evt_br['luminosityBlock'][iEvent]), int(Evt_br['event'][iEvent]))
            sRunLSEvent_toUse = sRunLSEvent.replace(":", "_") 
            if sRunLSEvent not in eventsToRun_list: continue
            print("Ruuning on selected event %s" % (sRunLSEvent))

        # print("Hi")

        puWeight = 1.0
        nVtx = Vtx_br['nObjects'][iEvent]

        if usePUReweighting:
            bin_puWeight = hPUWt.FindBin(nVtx);
            puWeight = hPUWt.GetBinContent(bin_puWeight);

        hnVtx.Fill(nVtx);    
        hnVtx_ReWtd.Fill(nVtx, puWeight);

        nOffJets  = int(Jet_br['nObjects'][iEvent])
        nOffMuons = int(Muon_br['nObjects'][iEvent])
        nOffEles  = int(Ele_br['nObjects'][iEvent])
        nUnpJets  = int(Unp_br['nObjects'][iEvent])
        nEmuJets  = int(Emu_br['nObjects'][iEvent])
        # nEmuHTPs  = int(eTP_br.nHCALTP[iEvent])
        # nEmuETPs  = int(eTP_br.nECALTP[iEvent])
        nEmuTTs   = int(eTT_br['nObjects'][iEvent])
        nUnpTTs   = int(uTT_br['nObjects'][iEvent])
        # nEmuTCs   = int(eTC_br.nCluster[iEvent])
        # nGenJets  = int(Gen_br.nJet[iEvent])

        offlinePUPPIJet = False

        if   l1MatchOffline:
            #if offlineJetType == 'PUPPI':
            if offlinePUPPIJet:
                nRefJets     = Jet_br.puppi_nJets
                et_RefJets   = Jet_br.puppi_etCorr
                eta_RefJets  = Jet_br.puppi_eta
                phi_RefJets  = Jet_br.puppi_phi
                #print(f"PUPPI jets")
            else:
                nRefJets     = Jet_br['nObjects'][iEvent]
                et_RefJets   = Jet_br['pt'][iEvent]
                eta_RefJets  = Jet_br['eta'][iEvent]
                phi_RefJets  = Jet_br['phi'][iEvent]                    
        elif l1MatchGen:
            print("USing Gen")
            nRefJets     = GenJet_br['nObjects'][iEvent]
            et_RefJets   = GenJet_br['pt'][iEvent]
            eta_RefJets  = GenJet_br['eta'][iEvent]
            phi_RefJets  = GenJet_br['phi'][iEvent]      

        hnVts_vs_nTT_unp.Fill(nVtx, nUnpTTs)
        hnVts_vs_nTT_emu.Fill(nVtx, nEmuTTs)

        nUnpJets_Bx0 = 0
        for iJ in range(nUnpJets): 
            if Unp_br['bx'][iEvent][iJ] != 0: continue   
            nUnpJets_Bx0 += 1        
            hL1JetUnp_Pt_0.Fill(Unp_br['pt'][iEvent][iJ])
            hL1JetUnp_Eta_0.Fill(Unp_br['eta'][iEvent][iJ])
            hL1JetUnp_Phi_0.Fill(Unp_br['phi'][iEvent][iJ])
        hnL1JetUnp_0.Fill(nUnpJets_Bx0)

        nEmuJets_Bx0 = 0
        for iJ in range(nEmuJets):
            if Emu_br['bx'][iEvent][iJ] != 0: continue
            nEmuJets_Bx0 += 1
            hL1JetEmu_Pt_0.Fill(Emu_br['pt'][iEvent][iJ])
            hL1JetEmu_Eta_0.Fill(Emu_br['eta'][iEvent][iJ])
            hL1JetEmu_Phi_0.Fill(Emu_br['phi'][iEvent][iJ])
        hnL1JetEmu_0.Fill(nEmuJets_Bx0)

        for iJ in range(nOffJets):
            if offlinePUPPIJet:
                hOfflineJet_Pt_0.Fill(Jet_br.puppi_etCorr[iJ])
                hOfflineJet_Eta_0.Fill(Jet_br.puppi_eta[iJ])
                hOfflineJet_Phi_0.Fill(Jet_br.puppi_phi[iJ])
            else:
                hOfflineJet_Pt_0.Fill(Jet_br['pt'][iEvent][iJ])
                hOfflineJet_Eta_0.Fill(Jet_br['eta'][iEvent][iJ])
                hOfflineJet_Phi_0.Fill(Jet_br['phi'][iEvent][iJ])                    
        hnOfflineJet_0.Fill(nOffJets)

        #Mimic trigger condition
        # if not isMC and len(HLT_Triggers_Required) > 0:
        #     passingTrigThshs = True  
        #     # Iterate through all triggers in HLT_Triggers_Required
        #     for trigger_type in HLT_Triggers_Required:
        #         if trigger_type == 'IsoMu24_OneProng32':
        #             nOffMuons_passingTrigThsh = [0] * len(TrigThshs_OffMuPt)
        #             for iMu in range(nOffMuons):
        #                 if not Muon_br['tightId'][iEvent][iMu]: continue
                        
        #                 for iTrigThsh in range(len(TrigThshs_OffMuPt)):
        #                     if Muon_br['pt'][iEvent][iMu] > TrigThshs_OffMuPt[iTrigThsh]:
        #                         nOffMuons_passingTrigThsh[iTrigThsh] += 1

        #             if not all(x > 0 for x in nOffMuons_passingTrigThsh):  
        #                 passingTrigThshs = False

        #         elif trigger_type == 'SingleJet180':
        #             nOffJets_passingTrigThsh = [0] * len(TrigThshs_OffJetPt)
        #             for iJet in range(nOffJets):
        #                 for iTrigThsh in range(len(TrigThshs_OffJetPt)):
        #                     if Jet_br['pt'][iEvent][iJet] > TrigThshs_OffJetPt[iTrigThsh]:
        #                         nOffJets_passingTrigThsh[iTrigThsh] += 1

        #             if not all(x > 0 for x in nOffJets_passingTrigThsh):  
        #                 passingTrigThshs = False

        #         else:
        #             print(f"HLT_Triggers_Required: {trigger_type} not implemented")
        #             passingTrigThshs = False

        #     if not passingTrigThshs:  # If one trigger fails, exit early
        #         break

        #     hStat.Fill(3)
        
        # -----------------------------------------------------------------------------------------------------------------------------

        ### save l1jets reconstructed with different JetShape+PUS for single/double/tripple/qud-jet trigger rates ---------------------

        l1JetCollection = OrderedDict()
        for src in ['unp','emu']:
            l1JetCollection[src] = OrderedDict()
            
            #for jetShape in ['Default'] + JetShapes:
            for jetShape in JetShapes + JetShapesType2:
                l1JetCollection[src][jetShape] = OrderedDict()
                
                for algo1 in PUSAlgosAll + PUSAlgosAllType2: # ['Raw', 'RawPUS', 'RawPUS_phiRingMin4', 'RawPUS_phiRingSide4', 'RawPUS_phiRingAdjacent']:
                    # read proper jetShape and PUSAlgo conbination
                    if (jetShape in JetShapes      and algo1 not in PUSAlgosAll) or \
                        (jetShape in JetShapesType2 and algo1 not in PUSAlgosAllType2 ):
                        continue
            
                    if (algo1 == 'L1JDefault') and (jetShape != 'Default'): continue
                    
                    l1JetCollection[src][jetShape][algo1] = OrderedDict()
                    
                    for ieta_cat in IETA_CAT.keys():
                        if ieta_cat == 'HBEF': continue
                        
                        l1JetCollection[src][jetShape][algo1][ieta_cat] = []

        ## Create list of offfline RECO jets which are too close to other RECO jets
        bad_off_jets = []
        ## Loop over all offline RECO jets
        #for iOff in range(nOffJets):
        for iOff in range(nRefJets):
            iOff_vec = R.TLorentzVector()
            iOff_vec.SetPtEtaPhiM(et_RefJets[iOff], eta_RefJets[iOff], phi_RefJets[iOff], 0)
            
            ## Loop over all offline RECO jets with higher pT
            for jOff in range(iOff):
                jOff_vec = R.TLorentzVector()
                jOff_vec.SetPtEtaPhiM(et_RefJets[jOff], eta_RefJets[jOff], phi_RefJets[jOff], 0)

                if iOff_vec.DeltaR(jOff_vec) < DR_MIN:
                    # print '\n  * Removing offline jet pT = %.1f, eta = %.2f, phi = %.2f' % (iOff_vec.Pt(), iOff_vec.Eta(), iOff_vec.Phi())
                    # print '  * Has dR = %.2f to jet pT = %.1f, eta = %.2f, phi = %.2f' % (iOff_vec.DeltaR(jOff_vec), jOff_vec.Pt(), jOff_vec.Eta(), jOff_vec.Phi())
                    bad_off_jets.append(iOff)
                    break

        isFirstRefJet = True
        for iOff in range(nRefJets):

            hStat.Fill(4)
            hRefJet_pt_0.Fill(et_RefJets[iOff])
            hRefJet_eta_0.Fill(eta_RefJets[iOff])
            hRefJet_phi_0.Fill(phi_RefJets[iOff])
            
            ## Remove offline jets which overlap other jets
            if iOff in bad_off_jets: continue

            hRefJet_pt_0_1.Fill(et_RefJets[iOff])
            hRefJet_eta_0_1.Fill(eta_RefJets[iOff])
            hRefJet_phi_0_1.Fill(phi_RefJets[iOff])                
            
            hStat.Fill(5)

            vOff = R.TLorentzVector()
            vOff.SetPtEtaPhiM(et_RefJets[iOff], eta_RefJets[iOff], phi_RefJets[iOff], 0)               
            ### save l1jets reconstructed with different JetShape+PUS for single/double/tripple/qud-jet trigger rates ---------------------

            
            if   l1MatchOffline:
                selectPFJet = True
                ## PF jet filters as recommended by Aaron: https://github.com/cms-l1t-offline/cms-l1t-analysis/blob/master/cmsl1t/filters/jets.py
                abs_eta = abs(eta_RefJets[iOff])
                reject_if = None

                if dataEra in [
                    '2022F', '2022G', 
                    '2023B', '2023C', '2023D', 
                    '2024A', '2024B', '2024C', '2024D', '2024E', '2024F', '2024G', '2024H', '2024I' ]:
                    # https://twiki.cern.ch/twiki/bin/view/CMS/JetID13p6TeV#Recommendations_for_the_13_6_AN1
                    # https://github.com/bundocka/cmssw/blob/7d536e034f7dd0773eec3f306508c80c67fb1960/L1Trigger/L1Tnanos/plugins/L1JetRecoTreeProducer.cc#L689-L715

                    isCentralJet          =  abs_eta <= 2.6
                    isForwardCentralJet_1 = (abs_eta > 2.6 and abs_eta <= 2.7)
                    isForwardCentralJet_2 = (abs_eta > 2.7 and abs_eta <= 3.0)
                    isForwardJet          =  abs_eta > 3.0
                    reject_if = [
                        isCentralJet          and Jet_br['neHEF'][iEvent][iOff]  >= 0.99, # neutralHadronEnergyFraction()
                        isCentralJet          and Jet_br['neEmEF'][iEvent][iOff] >= 0.90, # jet_data->nemef.push_back(it->neutralEmEnergyFraction());
                        isCentralJet          and (Jet_br['chMultiplicity'][iEvent][iOff] + Jet_br['neMultiplicity'][iEvent][iOff]) <= 1, # jet_data->cMult.push_back(it->chargedMultiplicity()); jet_data->nMult.push_back(it->neutralMultiplicity());
                        isCentralJet          and Jet_br['muEF'][iEvent][iOff]   >= 0.80, # jet_data->mef.push_back(it->muonEnergyFraction());
                        isCentralJet          and Jet_br['chHEF'][iEvent][iOff]  <= 0.01, # jet_data->chef.push_back(it->chargedHadronEnergyFraction());
                        isCentralJet          and Jet_br['chMultiplicity'][iEvent][iOff] == 0, # jet_data->cMult.push_back(it->chargedMultiplicity());
                        isCentralJet          and Jet_br['chEmEF'][iEvent][iOff] >= 0.80, # jet_data->cemef.push_back(it->chargedEmEnergyFraction());

                        isForwardCentralJet_1 and Jet_br['neHEF'][iEvent][iOff]  >= 0.90, # neutralHadronEnergyFraction()
                        isForwardCentralJet_1 and Jet_br['neEmEF'][iEvent][iOff] >= 0.99, # jet_data->nemef.push_back(it->neutralEmEnergyFraction());
                        isForwardCentralJet_1 and Jet_br['muEF'][iEvent][iOff]   >= 0.80, # jet_data->mef.push_back(it->muonEnergyFraction());
                        isForwardCentralJet_1 and Jet_br['chEmEF'][iEvent][iOff] >= 0.80, # jet_data->cemef.push_back(it->chargedEmEnergyFraction());

                        isForwardCentralJet_2 and Jet_br['neEmEF'][iEvent][iOff] >= 0.99, # jet_data->nemef.push_back(it->neutralEmEnergyFraction());                             

                        ##isForwardJet          and Jet_br.nhef[iOff]  <= 0.20, # neutralHadronEnergyFraction()
                        isForwardJet          and Jet_br['neEmEF'][iEvent][iOff] >= 0.40, # jet_data->nemef.push_back(it->neutralEmEnergyFraction());
                        isForwardJet          and Jet_br['neMultiplicity'][iEvent][iOff] <= 2, # jet_data->nMult.push_back(it->neutralMultiplicity());
                    ]
                    # print(reject_if)

                if reject_if and any(reject_if):  # Ensure reject_if is not None
                    selectPFJet = False

                if not selectPFJet: continue

                hRefJet_pt_0_2.Fill(et_RefJets[iOff])
                hRefJet_eta_0_2.Fill(eta_RefJets[iOff])
                hRefJet_phi_0_2.Fill(phi_RefJets[iOff])     

                # OfflineJet check with OfflineMuon for overlap ---------------------------------------
                passingJetMuOverlap = False
                dR_OffJet_OffMu_min = 99999.0
                idx_OffMu_nearestToOffJet = -1
                nOffMuons= Muon_br['nObjects'][iEvent]
                nOffEles= Ele_br['nObjects'][iEvent]
                for iMu in range(nOffMuons):
                    # muon selection
                    # if not Muon_br.isMediumMuon[iEvent][iMu]: continue
                    if not Muon_br['mediumId'][iEvent][iMu]: continue
                    
                    vOffMu = R.TLorentzVector()
                    vOffMu.SetPtEtaPhiM(Muon_br['pt'][iEvent][iMu], Muon_br['eta'][iEvent][iMu], Muon_br['phi'][iEvent][iMu], MASS_MUON)

                    # MASS_ELECTRON
                    dr_tmp = vOff.DeltaR(vOffMu)
                    if dr_tmp < dR_OffJet_OffMu_min:
                        dR_OffJet_OffMu_min = dr_tmp
                        idx_OffMu_nearestToOffJet = iMu

                    if dr_tmp < DR_Jet_Ele_Min and Muon_br['pt'][iEvent][iMu] / vOff.Pt() > RATIO_PtEle_PtJet_Max:
                        passingJetMuOverlap = True
                        

                if idx_OffMu_nearestToOffJet != -1:
                    pTFraction_tmp = Muon_br['pt'][iEvent][idx_OffMu_nearestToOffJet] / vOff.Pt()
                    hdR_OffJet_OffMu_min.Fill(dR_OffJet_OffMu_min)
                    hdR_OffJet_OffMu_min_vs_vOffMuPtByvOffJetPt.Fill(dR_OffJet_OffMu_min, pTFraction_tmp)
                    if pTFraction_tmp < 0.5:
                        hdR_OffJet_OffMu_min_forPtFracLt0p5.Fill(dR_OffJet_OffMu_min)
                    
                    
                # OfflineJet check with OfflineElectron for overlap ---------------------------------------
                passingJetEleOverlap = False
                dR_OffJet_OffEle_min = 99999.0
                idx_OffEle_nearestToOffJet = -1
                for iEle in range(nOffEles):
                    # electron selection
                    # if not Ele_br.isMediumElectron[iEvent][iEle]: continue
                    
                    vOffEle = R.TLorentzVector()
                    vOffEle.SetPtEtaPhiM(Ele_br['pt'][iEvent][iEle], Ele_br['eta'][iEvent][iEle], Ele_br['phi'][iEvent][iEle], MASS_ELECTRON)

                    # MASS_ELECTRON
                    dr_tmp = vOff.DeltaR(vOffEle)
                    if dr_tmp < dR_OffJet_OffEle_min:
                        dR_OffJet_OffEle_min = dr_tmp
                        idx_OffEle_nearestToOffJet = iEle

                    if dr_tmp < DR_Jet_Ele_Min and Ele_br['pt'][iEvent][iEle] / vOff.Pt() > RATIO_PtEle_PtJet_Max:
                        passingJetEleOverlap = True


                if idx_OffEle_nearestToOffJet != -1:
                    pTFraction_tmp = Ele_br['pt'][iEvent][idx_OffEle_nearestToOffJet] / vOff.Pt()
                    hdR_OffJet_OffEle_min.Fill(dR_OffJet_OffEle_min)
                    hdR_OffJet_OffEle_min_vs_vOffElePtByvOffJetPt.Fill(dR_OffJet_OffEle_min, pTFraction_tmp)
                    if pTFraction_tmp < 0.5:
                        hdR_OffJet_OffEle_min_forPtFracLt0p5.Fill(dR_OffJet_OffEle_min)


                if passingJetMuOverlap: continue

                hStat.Fill(7)

                hRefJet_pt_0_3.Fill(et_RefJets[iOff])
                hRefJet_eta_0_3.Fill(eta_RefJets[iOff])
                hRefJet_phi_0_3.Fill(phi_RefJets[iOff])                

                if passingJetEleOverlap: continue

                hStat.Fill(8)

                hRefJet_pt_0_4.Fill(et_RefJets[iOff])
                hRefJet_eta_0_4.Fill(eta_RefJets[iOff])
                hRefJet_phi_0_4.Fill(phi_RefJets[iOff])

            # 2018 data: HE- dead zone (HEM15/16)
            if not isMC:
                # no PUPPI jets were stored in L1nanos for run2 data
                if Evt_br['run'][iEvent] > 319077 and Evt_br['run'][iEvent] < 340000 and \
                    Jet_br['eta'][iEvent][iOff] > -3.4  and Jet_br['eta'][iEvent][iOff] < -1.17 and \
                    Jet_br['phi'][iEvent][iOff] > -1.97 and Jet_br['phi'][iEvent][iOff] < -0.47:
                    continue

            hStat.Fill(9)

            # save Reco jet pT for Gen jet ------------------------------------------------------------

            if l1MatchGen:
                off_vec_matchedTo_gen_vec = R.TLorentzVector()
                off_vec_matchedTo_gen_vec.SetPtEtaPhiM(0, 0, 0, 0)
                for iOff in range(nOffJets):
                    iOff_vec = R.TLorentzVector()
                    iOff_vec.SetPtEtaPhiM(Jet_br['pt'][iEvent][iOff], Jet_br['eta'][iEvent][iOff], Jet_br['phi'][iEvent][iOff], 0)
                    if vOff.DeltaR(iOff_vec) < DR_MAX and iOff_vec.Pt() > off_vec_matchedTo_gen_vec.Pt():
                        off_vec_matchedTo_gen_vec = iOff_vec

            # TT, TP plots --------------------------------------------------------------------------------------------------------------------
            if isMC and useCutGenNVtxEq0 and (nRefJets == 1):
                for src in ['emu']: #['unp','emu']:
                    TT_br_toUse = None
                    TP_br_toUse = None
                    TC_br_toUse = None
                    if src == 'unp':
                        TT_br_toUse = uTT_br 
                        TP_br_e_toUse = uTP_e_br
                        TP_br_h_toUse = uTP_h_br
                        TC_br_toUse = uTC_br
                    else:
                        TT_br_toUse = eTT_br
                        TP_br_e_toUse = eTP_e_br
                        TP_br_h_toUse = eTP_h_br
                        TC_br_toUse = eTC_br

                    TT28Abv125 = False
                    RefJetEtaAbs = abs(vOff.Eta())
                    if RefJetEtaAbs > 2.5 and RefJetEtaAbs < 3.1:
                        for iTT in range(TT_br_toUse['nObjects'][iEvent]):
                            sIEta = str(abs(TT_br_toUse['ieta'][iEvent][iTT]))
                            hist8['TT_iet_RefJetEtaBtw2p5And3p1'][src][sIEta].Fill(TT_br_toUse['iet'][iEvent][iTT])
                            
                            if abs(TT_br_toUse['ieta'][iEvent][iTT]) in [ 28] and TT_br_toUse['iet'][iEvent][iTT] > 125:
                                    TT28Abv125 = True
                                    

                        # for iTP in range(TP_br_e_toUse.nECALTP[iEvent]):
                        #     sIEta = str(abs(TP_br_e_toUse['ieta'][iEvent][iTP]))
                        #     hist8['ECALTP_et_RefJetEtaBtw2p5And3p1'    ][src][sIEta].Fill(TP_br_e_toUse.et[iEvent][iTP])
                        #     hist8['ECALTP_compEt_RefJetEtaBtw2p5And3p1'][src][sIEta].Fill(TP_br_e_toUse.compEt[iEvent][iTP])

                        # for iTP in range(TP_br_h_toUse.nHCALTP[iEvent]):
                        #     sIEta = str(abs(TP_br_h_toUse['ieta'][iEvent][iTP]))
                        #     hist8['HCALTP_et_RefJetEtaBtw2p5And3p1'    ][src][sIEta].Fill(TP_br_h_toUse.et[iEvent][iTP])
                        #     hist8['HCALTP_compEt_RefJetEtaBtw2p5And3p1'][src][sIEta].Fill(TP_br_h_toUse.compEt[iEvent][iTP])

                            
                    if RefJetEtaAbs < 2.5:
                        for iTT in range(TT_br_toUse['nObjects'][iEvent]):
                            sIEta = str(abs(TT_br_toUse['ieta'][iEvent][iTT]))
                            hist8['TT_iet_RefJetEtalt2p5'][src][sIEta].Fill(TT_br_toUse['iet'][iEvent][iTT])

                        # for iTP in range(TP_br_e_toUse.nECALTP[iEvent]):
                        #     sIEta = str(abs(TP_br_e_toUse['ieta'][iEvent][iTP]))
                        #     hist8['ECALTP_et_RefJetEtalt2p5'    ][src][sIEta].Fill(TP_br_e_toUse.et[iEvent][iTP])
                        #     hist8['ECALTP_compEt_RefJetEtalt2p5'][src][sIEta].Fill(TP_br_e_toUse.compEt[iEvent][iTP])

                        # for iTP in range(TP_br_h_toUse.nHCALTP[iEvent]):
                        #     sIEta = str(abs(TP_br_h_toUse['ieta'][iEvent][iTP]))
                        #     hist8['HCALTP_et_RefJetEtalt2p5'    ][src][sIEta].Fill(TP_br_h_toUse.et[iEvent][iTP])
                        #     hist8['HCALTP_compEt_RefJetEtalt2p5'][src][sIEta].Fill(TP_br_h_toUse.compEt[iEvent][iTP])

                    # if PrintLevel >= 0:
                    #     if RefJetEtaAbs > 2.5 and RefJetEtaAbs < 3.1  and TT28Abv125:
                    #         print(f"RefJet: pt {vOff.Pt()}, eta {vOff.Eta()}, phi {vOff.Phi()}.   {src}, {nEmuTTs = }, {nEmuTCs = }, {nEmuETPs = }, {nEmuHTPs = } ")
                    #         for iTT in range(TT_br_toUse.nTower):
                    #             #if abs(TT_br_toUse.ieta[iTT]) not in [27, 28]: continue
                    #             print(f"    TT {iTT}: iet {TT_br_toUse.iet[iTT]},  ieta {TT_br_toUse.ieta[iTT]}, iphi {TT_br_toUse.iphi[iTT]},  iem {TT_br_toUse.iem[iTT]}, ihad {TT_br_toUse.ihad[iTT]}, iqual {TT_br_toUse.iqual[iTT]}")
                    #         print(" ")

                    #         for iTC in range(TC_br_toUse.nCluster):
                    #             #if abs(TC_br_toUse.ieta[iTC]) not in [27, 28]: continue
                    #             print(f"    TC {iTC}: iet {TC_br_toUse.iet[iTC]}, ieta {TC_br_toUse.ieta[iTC]}, iphi {TC_br_toUse.iphi[iTC]},")
                    #         print(" ")
                                    
                    #         for iTP in range(TP_br_toUse.nECALTP):
                    #             #if abs(TP_br_toUse.ecalTPieta[iTP]) not in [27, 28]: continue
                    #             print(f"   TPECAL {iTP}: ecalTPet {TP_br_toUse.ecalTPet[iTP]}, ecalTPcompEt {TP_br_toUse.ecalTPcompEt[iTP]}, ecalTPieta {TP_br_toUse.ecalTPieta[iTP]} , ecalTPiphi {TP_br_toUse.ecalTPiphi[iTP]}, ecalTPCaliphi {TP_br_toUse.ecalTPCaliphi[iTP]} ")
                    #         print(" ")
                                
                    #         for iTP in range(TP_br_toUse.nHCALTP):
                    #             #if abs(TP_br_toUse.hcalTPieta[iTP]) not in [27, 28]: continue
                    #             print(f"   TPHCAL {iTP}: hcalTPet {TP_br_toUse.hcalTPet[iTP]}, hcalTPcompEt {TP_br_toUse.hcalTPcompEt[iTP]}, hcalTPieta {TP_br_toUse.hcalTPieta[iTP]}, hcalTPiphi {TP_br_toUse.hcalTPiphi[iTP]}, hcalTPCaliphi {TP_br_toUse.hcalTPCaliphi[iTP]}")

                                
            # # -----------------------------------------------------------------------------------------------------------------------------
            

            hRefJet_pt_1.Fill(vOff.Pt())
            hRefJet_eta_1.Fill(vOff.Eta())
            hRefJet_phi_1.Fill(vOff.Phi())

            if isFirstRefJet:
                isFirstRefJet = False

                nUnpJets_Bx0 = 0
                for iJ in range(nUnpJets): 
                    if Unp_br['bx'][iEvent][iJ] != 0: continue   
                    nUnpJets_Bx0 += 1        
                    hL1JetUnp_Pt_1.Fill(Unp_br['pt'][iEvent][iJ])
                    hL1JetUnp_Eta_1.Fill(Unp_br['eta'][iEvent][iJ])
                    hL1JetUnp_Phi_1.Fill(Unp_br['phi'][iEvent][iJ])
                hnL1JetUnp_1.Fill(nUnpJets_Bx0)

                nEmuJets_Bx0 = 0
                for iJ in range(nEmuJets):
                    if Emu_br['bx'][iEvent][iJ] != 0: continue
                    nEmuJets_Bx0 += 1
                    hL1JetEmu_Pt_1.Fill(Emu_br['pt'][iEvent][iJ])
                    hL1JetEmu_Eta_1.Fill(Emu_br['eta'][iEvent][iJ])
                    hL1JetEmu_Phi_1.Fill(Emu_br['phi'][iEvent][iJ])
                hnL1JetEmu_1.Fill(nEmuJets_Bx0)

            # plot pT, eta, phi of L1JetUnp jet matched to RefJet 
            idxL1Jet_matchedRefJet = -1
            L1JetPt_matchedRefJet = 0
            for iJ in range(nUnpJets): 
                if Unp_br['bx'][iEvent][iJ] != 0: continue 
                vL1Jet_ = R.TLorentzVector()
                vL1Jet_.SetPtEtaPhiM(Unp_br['pt'][iEvent][iJ], Unp_br['eta'][iEvent][iJ], Unp_br['phi'][iEvent][iJ], 0) 
                if vL1Jet_.DeltaR(vOff) < DR_MAX and vL1Jet_.Pt() > L1JetPt_matchedRefJet:
                    idxL1Jet_matchedRefJet = iJ
                    L1JetPt_matchedRefJet = vL1Jet_.Pt()
            if idxL1Jet_matchedRefJet >= 0:
                hL1JetUnp_Pt_2.Fill(Unp_br['pt'][iEvent][idxL1Jet_matchedRefJet])
                hL1JetUnp_Eta_2.Fill(Unp_br['eta'][iEvent][idxL1Jet_matchedRefJet])
                hL1JetUnp_Phi_2.Fill(Unp_br['phi'][iEvent][idxL1Jet_matchedRefJet])                                            

            # plot pT, eta, phi of L1JetEmu jet matched to RefJet 
            idxL1Jet_matchedRefJet = -1
            L1JetPt_matchedRefJet = 0
            for iJ in range(nEmuJets): 
                if Emu_br['bx'][iEvent][iJ] != 0: continue 
                vL1Jet_ = R.TLorentzVector()
                vL1Jet_.SetPtEtaPhiM(Emu_br['pt'][iEvent][iJ], Emu_br['eta'][iEvent][iJ], Emu_br['phi'][iEvent][iJ], 0) 
                if vL1Jet_.DeltaR(vOff) < DR_MAX and vL1Jet_.Pt() > L1JetPt_matchedRefJet:
                    idxL1Jet_matchedRefJet = iJ
                    L1JetPt_matchedRefJet = vL1Jet_.Pt()
            if idxL1Jet_matchedRefJet >= 0:
                hL1JetEmu_Pt_2.Fill(Emu_br['pt'][iEvent][idxL1Jet_matchedRefJet])
                hL1JetEmu_Eta_2.Fill(Emu_br['eta'][iEvent][idxL1Jet_matchedRefJet])
                hL1JetEmu_Phi_2.Fill(Emu_br['phi'][iEvent][idxL1Jet_matchedRefJet])       

            jetIEta_offlineJet     = calculateJetIEta(vOff.Eta())
            jetIEtaAbs_offlineJet  = abs(jetIEta_offlineJet)
            sjetIEta_offlineJet    = "%d" % (int(jetIEta_offlineJet))
            sjetIEtaAbs_offlineJet = "%d" % (int(jetIEtaAbs_offlineJet))
            #if sjetIEta_offlineJet in hist_PFJetPt_iEtawise:
            #    hist_PFJetPt_iEtawise[sjetIEta_offlineJet].Fill(vOff.Pt(), puWeight )


            data_dict = OrderedDict()
            #if VERBOSE and iEvt % PRT_EVT is 0: print '  * Run %d, LS %d, event %d, nVtx %d' % (int(Evt_br['run']), int(Evt_br['luminosityBlock']), int(Evt_br['event']), int(nVtx))
            data_dict['runNumber']                = int(Evt_br['run'][iEvent])
            data_dict['lumiSectionNumber']        = int(Evt_br['luminosityBlock'][iEvent])
            data_dict['eventNumber']              = int(Evt_br['event'][iEvent])
            data_dict['nVertexReco']              = int(nVtx)
            # data_dict['nTT_Unpacked']             = nUnpTTs
            data_dict['nTT_Emulated']             = nEmuTTs
            if   l1MatchOffline:
                data_dict['PFJetEtCorr']          = vOff.Pt()
            elif l1MatchGen:
                data_dict['GenJetEt']             = vOff.Pt()
                data_dict['nVertexGen']           = int(nVtx)
                # data_dict['nMeanPUGen']           = int(Gen_br.nMeanPU)
                data_dict['matchedPFJetEtCorr']   = off_vec_matchedTo_gen_vec.Pt()
        
            if runMode not in ['CalCalibSF', 'makeInputForML'] and vOff.Pt() < PT_MIN: continue
        
            hStat.Fill(10)

            PFJetEtaCat = 'None'
            for iCat in ETA_CAT.keys():
                if iCat == 'HBEF': continue
                if abs(vOff.Eta()) > ETA_CAT[iCat][0] and abs(vOff.Eta()) < ETA_CAT[iCat][1]:
                    PFJetEtaCat = iCat
            if PFJetEtaCat == 'None' or PFJetEtaCat == 'HBEF':
                if l1MatchOffline:
                    print('\n\nSUPER-BIZZARE JET THAT FALLS INTO NO ETA CATEGORIES!!!  eta = %.3f\n\n' % vOff.Eta())
                continue
            sEtaCat_PFJet = PFJetEtaCat
            
            iPFJetPtCat = 'None'
            for iCat in PT_CAT.keys():
                if vOff.Pt() > PT_CAT[iCat][0] and vOff.Pt() < PT_CAT[iCat][2]:
                    iPFJetPtCat = iCat
            if iPFJetPtCat == 'None':
                if vOff.Pt() > PT_CAT['lowPt'][0] and vOff.Pt() < PT_CAT['hiPt'][2]:
                    print('\n\nSUPER-BIZZARE JET THAT FALLS INTO NO PT CATEGORIES!!!  pT = %.3f\n\n' % vOff.Pt())
                continue
            
            if VERBOSE or PrintLevel >= 1:
                print("    offlineJet: eta: {}, phi: {}, etCorr: {}".format(eta_RefJets[iOff], phi_RefJets[iOff], et_RefJets[iOff]))
        
            hStat.Fill(11)

            ## Find highest-pT Level-1 jet with good dR matching to unpacked jet
            max_pt = {}
            vMax   = {}
            matchedEmuIdx = {}
            for src in ['unp','emu']:
                max_pt[src] = {}
                vMax  [src] = {}
                matchedEmuIdx[src] = {}
                for algo in ['PUS','noPUS','Raw','RawPUS']:
                    max_pt[src][algo] = -99
                    vMax  [src][algo] = R.TLorentzVector()
                    matchedEmuIdx[src][algo] = -1
            ## Loop over all L1T unpacked jets
            if PrintLevel >= 12:
                print("  * UnpJets ({}):: ".format(nUnpJets))
            for iUnp in range(nUnpJets):

                if Unp_br['bx'][iEvent][iUnp] != 0: continue  ## Use only jets in BX 0
                
                hStat.Fill(12)

                vUnp = {}
                for algo in ['PUS','noPUS','Raw','RawPUS']:
                    vUnp[algo] = R.TLorentzVector()  ## Create a 4-vector of the L1T jet

                vUnp['PUS']   .SetPtEtaPhiM(Unp_br['pt'][iEvent][iUnp],                           Unp_br['eta'][iEvent][iUnp], Unp_br['phi'][iEvent][iUnp], 0)
                vUnp['noPUS'] .SetPtEtaPhiM(Unp_br['pt'][iEvent][iUnp]    + Unp_br['puEt'][iEvent][iUnp], Unp_br['pt'][iEvent][iUnp], Unp_br['phi'][iEvent][iUnp], 0)
                vUnp['Raw']   .SetPtEtaPhiM(Unp_br['rawEt'][iEvent][iUnp],                        Unp_br['eta'][iEvent][iUnp], Unp_br['phi'][iEvent][iUnp], 0)
                vUnp['RawPUS'].SetPtEtaPhiM(Unp_br['rawEt'][iEvent][iUnp] - Unp_br['puEt'][iEvent][iUnp], Unp_br['eta'][iEvent][iUnp], Unp_br['phi'][iEvent][iUnp], 0)

                for algo in ['PUS','noPUS','Raw','RawPUS']:
                    if vUnp[algo].DeltaR(vOff) < DR_MAX and vUnp[algo].Pt() > max_pt['unp'][algo]:
                        max_pt['unp'][algo] = vUnp[algo].Pt()
                        vMax  ['unp'][algo] = vUnp[algo]
                        matchedEmuIdx['unp'][algo] = iUnp
                
                #hist_L1Jet_unp_TowerIEta_vs_IEta.Fill(Unp_br.jetTowerIEta[iUnp], Unp_br.jetIEta[iUnp])
                #hist_L1Jet_unp_TowerIPhi_vs_IPhi.Fill(Unp_br.jetTowerIPhi[iUnp], Unp_br.jetIPhi[iUnp]) 
                         
                
            ## End loop: for iUnp in range(nUnpJets)

            ## Loop over all L1T emulated jets
            for iEmu in range(nEmuJets):
                if Emu_br['bx'][iEvent][iEmu] != 0: continue  ## Use only jets in BX 0

                hStat.Fill(13)
                
                vEmu = {}
                for algo in ['PUS','noPUS','Raw','RawPUS']:
                    vEmu[algo] = R.TLorentzVector()  ## Create a 4-vector of the L1T jet

                vEmu['PUS']   .SetPtEtaPhiM(Emu_br['pt'][iEvent][iEmu],                           Emu_br['eta'][iEvent][iEmu], Emu_br['phi'][iEvent][iEmu], 0)
                vEmu['noPUS'] .SetPtEtaPhiM(Emu_br['pt'][iEvent][iEmu]    + Emu_br['puEt'][iEvent][iEmu], Emu_br['eta'][iEvent][iEmu], Emu_br['phi'][iEvent][iEmu], 0)
                vEmu['Raw']   .SetPtEtaPhiM(Emu_br['rawEt'][iEvent][iEmu],                        Emu_br['eta'][iEvent][iEmu], Emu_br['phi'][iEvent][iEmu], 0)
                vEmu['RawPUS'].SetPtEtaPhiM(Emu_br['rawEt'][iEvent][iEmu] - Emu_br['puEt'][iEvent][iEmu], Emu_br['eta'][iEvent][iEmu], Emu_br['phi'][iEvent][iEmu], 0)

                for algo in ['PUS','noPUS','Raw','RawPUS']:
                    if vEmu[algo].DeltaR(vOff) < DR_MAX and vEmu[algo].Pt() > max_pt['emu'][algo]:
                        max_pt['emu'][algo] = vEmu[algo].Pt()
                        vMax  ['emu'][algo] = vEmu[algo]
                        matchedEmuIdx['emu'][algo] = iEmu

                #hist_L1Jet_emu_TowerIEta_vs_IEta.Fill(Emu_br.jetTowerIEta[iEvent][iEmu], Emu_br.jetIEta[iEvent][iEmu])
                #hist_L1Jet_emu_TowerIPhi_vs_IPhi.Fill(Emu_br.jetTowerIPhi[iEvent][iEmu], Emu_br.jetIPhi[iEvent][iEmu])
        
            ## Re-set the |eta| categories based on emulated and unpacked L1T jet eta, if there is a match
            etaCat = {}
            for src in ['unp','emu']:
                etaCat[src] = {}
                for algo in ['PUS','noPUS','Raw','RawPUS']:
                    etaCat[src][algo] = 'None'

                    if max_pt[src][algo] > 0:
                        for iCat in ETA_CAT.keys():
                            if abs(vMax[src][algo].Eta()) > ETA_CAT[iCat][0] and abs(vMax[src][algo].Eta()) < ETA_CAT[iCat][1]:
                                etaCat[src][algo] = iCat
                    else:       etaCat[src][algo] = PFJetEtaCat

                    if etaCat[src][algo] == None:
                        print('\n\nSUPER-BIZZARE JET THAT FALLS INTO NO ETA CATEGORIES!!!  eta[%s][%s] = %.3f\n\n' % (vMax[src][algo].Eta(), src, algo))
                        continue
                    # if etaCat[src][algo] != PFJetEtaCat:
                    #     print '  * L1T jet (eta = %.3f) not in same category as RECO jet (eta = %.3f)' % (vMax[src][algo].Eta(), vOff.Eta())


            for src in ['unp','emu']:
                for algo in ['PUS','noPUS','Raw','RawPUS']:
                    for jPt in PT_CAT.keys():
                        if 'jet_den' in hist.keys():
                            hist['jet_den'][algo][etaCat[src][algo]][jPt][src][iEvent].Fill( vOff.Pt() )
                            hist['jet_den'][algo]['HBEF'           ][jPt][src][iEvent].Fill( vOff.Pt() )
                        if vMax[src][algo].Pt() > PT_CAT[jPt][1]:
                            if 'jet_num' in hist.keys():
                                hist['jet_num'][algo][etaCat[src][algo]][jPt][src][iEvent].Fill( vOff.Pt() )
                                hist['jet_num'][algo]['HBEF'           ][jPt][src][iEvent].Fill( vOff.Pt() )

                    if max_pt[src][algo] > 0:
                        # iL1JetPtCat = getJetPtCategory( vMax[src][algo].Pt() )
                        # etaCat[src][algo]
                        if 'jet_res' in hist.keys():
                            hist['jet_res'][algo][PFJetEtaCat][iPFJetPtCat][src][iEvent].Fill( (vMax[src][algo].Pt() - vOff.Pt()) / vOff.Pt(), puWeight )
                            hist['jet_res'][algo]['HBEF'     ][iPFJetPtCat][src][iEvent].Fill( (vMax[src][algo].Pt() - vOff.Pt()) / vOff.Pt(), puWeight )
                        if 'jet_dR' in hist.keys():
                            hist['jet_dR'] [algo][PFJetEtaCat][iPFJetPtCat][src][iEvent].Fill( (vMax[src][algo].DeltaR(vOff)), puWeight )
                            hist['jet_dR'] [algo]['HBEF'     ][iPFJetPtCat][src][iEvent].Fill( (vMax[src][algo].DeltaR(vOff)), puWeight )

                ## End loop: for algo in ['PUS','noPUS','Raw','RawPUS']
            hStat.Fill(14)

            # for src in ['unp','emu']: 
            #     l1TP_br = None
            #     if src == 'unp':
            #         l1TP_br = uTP_br
            #     elif src == 'emu':
            #         l1TP_br = eTP_br                               

            #     # if PrintLevel >= 8:
            #     #     print("{} {} HCAL TP ({})".format(" "*4,src, l1TP_br.nECALTP))
            #     #     for iTP in range(l1TP_br.nECALTP):
            #     #         print(" {} ecalTPieta {}".format( iTP, l1TP_br.ecalTPieta[iTP]))
            #     #         print(" {} ecalTPiphi {}".format( iTP, l1TP_br.ecalTPiphi[iTP]))
            #     #         print(" {} ecalTPiCaliphi {}".format( iTP, l1TP_br.ecalTPCaliphi[iTP]))
            #     #         print(" {} ecalTPet {}".format( iTP, l1TP_br.ecalTPet[iTP]))
            #     #         print(" {} ecalTPcompEt {}".format( iTP, l1TP_br.ecalTPcompEt[iTP]))
            #     #         print(" {} ecalTPfineGrain {}".format( iTP, l1TP_br.ecalTPfineGrain[iTP]))
                        
            #     #         print("{} {}: ecalTPieta {}, ecalTPiphi {}, ecalTPCaliphi {}, ecalTPet {}, ecalTPcompEt {}, ecalTPfineGrain {}".format(" "*6, iTP, l1TP_br.ecalTPieta[iTP], l1TP_br.ecalTPiphi[iTP], l1TP_br.ecalTPCaliphi[iTP], l1TP_br.ecalTPet[iTP], l1TP_br.ecalTPcompEt[iTP], l1TP_br.ecalTPfineGrain[iTP]))
                    
            #     #     print("{} {} HCAL TP ({})".format(" "*4,src, l1TP_br.nECALTP))
            #     #     for iTP in range(l1TP_br.nHCALTP):
            #     #         print("{} {}: hcalTPieta {}, hcalTPiphi {}, hcalTPCaliphi {}, hcalTPet {}, hcalTPcompEt {}, hcalTPfineGrain {}".format(" "*6, iTP, l1TP_br.hcalTPieta[iTP], l1TP_br.hcalTPiphi[iTP], l1TP_br.hcalTPCaliphi[iTP], l1TP_br.hcalTPet[iTP], l1TP_br.hcalTPcompEt[iTP], l1TP_br.hcalTPfineGrain[iTP]))
                    
            #     for iTP in range(l1TP_br.nECALTP):
            #         if 'ECAP_TP_et_vs_iEta_vs_nVts' in dists6:
            #             hist6['ECAP_TP_et_vs_iEta_vs_nVts'][src][iCh].    Fill( l1TP_br.ecalTPieta[iTP], nVtx, l1TP_br.ecalTPet[iTP],     puWeight )
            #         if 'ECAP_TP_compEt_vs_iEta_vs_nVts' in dists6:
            #             hist6['ECAP_TP_compEt_vs_iEta_vs_nVts'][src][iCh].Fill( l1TP_br.ecalTPieta[iTP], nVtx, l1TP_br.ecalTPcompEt[iTP], puWeight )
                
            #     for iTP in range(l1TP_br.nHCALTP):
            #         if 'HCAP_TP_et_vs_iEta_vs_nVts' in dists6:
            #             hist6['HCAP_TP_et_vs_iEta_vs_nVts'][src][iCh].    Fill( l1TP_br.hcalTPieta[iTP], nVtx, l1TP_br.hcalTPet[iTP],     puWeight )
            #         if 'HCAP_TP_compEt_vs_iEta_vs_nVts' in dists6:
            #             hist6['HCAP_TP_compEt_vs_iEta_vs_nVts'][src][iCh].Fill( l1TP_br.hcalTPieta[iTP], nVtx, l1TP_br.hcalTPcompEt[iTP], puWeight )

                


            hStat.Fill(30)        

            if not JetClustByHand:
                continue

            # run jet clustring by hand ---------------------------------------------
            if VERBOSE or PrintLevel >= 1:
                sTmp = "matchedEmuIdx: "
                for src in ['unp','emu']:
                    for algo in ['PUS','noPUS','Raw','RawPUS']:
                        sTmp += "  %s_%s %d" % (src,algo,matchedEmuIdx[src][algo])
                print("        {}".format(sTmp))

            # check clusters/TTs of the emulated/unpacked jet that matched to the offline jet
            # emulated jet index: l1jet_idx = matchedEmuIdx[src]['PUS']
            #for src in ['unp', 'emu']: # ['unp','emu'] # 'unp' doesn't work as jetTOwerIEta=0 variables are stored in unpacked branch
            for src in ['emu']:
                l1jet_br = None
                l1TC_br  = None
                l1TT_br  = None
                if src == 'unp':
                    l1jet_br = Unp_br
                    l1TT_br  = uTT_br
                elif src == 'emu':
                    l1jet_br = Emu_br
                    l1TC_br  = eTC_br
                    l1TT_br  = eTT_br

                # use l1 jet, leading in pT with algo='PUS' that matches to offline jet, 
                # as a reference (for jetToweriEta, jetTowerIPhi) to form cluster around (jetToweriEta, jetTowerIPhi).   
                l1jet_idx = matchedEmuIdx[src]['PUS']
                if l1jet_idx < 0: # no dR matching between emulated/unpacked jet and offline jet is found
                    res_dummy          = -1.49  # (l1jet_pt - vOff.Pt()) / vOff.Pt()
                    #jetIEta_offlineJet = calculateJetIEta(vOff.Eta())  # -50. # None # abs(vOff.Eta())
                    '''
                    for iEta, etaBinRange in map_iEta_Eta.items():
                        if abs(vOff.Eta()) >= etaBinRange[0] and abs(vOff.Eta()) < etaBinRange[1]:
                            jetIEta_offlineJet = float( iEta * math.copysign(1, vOff.Eta()) )
                    #print "    * Matched emulated jet not find. Offline jet eta: {}, iEta: {}".format(vOff.Eta(), jetIEta_offlineJet)
                    '''
                    
                    # fill dummy value in jet resolution histograms
                    if jetIEtaAbs_offlineJet <= 41:  # skip jetIEta_offlineJet in HF, hence not set
                        #for jetShape in ['Default'] + JetShapes:
                        for jetShape in JetShapes:
                            # JetShape = "" plots are with the first version of code for 9x9 jets
                            jetShape1 = jetShape
                            if jetShape == 'Default':  jetShape1 = ""
                            else:                      jetShape1 = "_%s" % (jetShape)
                        
                            for algo1 in PUSAlgosAll: # ['Raw', 'RawPUS', 'RawPUS_phiRingMin4', 'RawPUS_phiRingSide4', 'RawPUS_phiRingAdjacent']: # ['PUS','noPUS','Raw','RawPUS']:
                                if (algo1 == 'L1JDefault') and (jetShape != 'Default'): continue
                                if (not l1nanoChunkyDonut) and (jetShape == 'Default') and (algo1 == 'RawPUS'):            continue
                                if (not l1nanoPhiRing)     and (jetShape == 'Default') and (algo1 == 'RawPUS_phiDefault'): continue
                                
                                if jetIEta_offlineJet < -41: break # jetIEta_offlineJet is in HF, hence not set

                                jetIEta_offlineJet_tmp = jetIEtaAbs_offlineJet if useAbsEtaBins else jetIEta_offlineJet
                                #hist2['jet_byHand_res_vs_iEta%s' % (jetShape1)][algo1]['HBEF'      ][iPFJetPtCat        ][src][iEvent].Fill(jetIEta_offlineJet_tmp, res_dummy, puWeight)
                                #hist2['jet_byHand_res_vs_iEta%s' % (jetShape1)][algo1]['HBEF'      ]['PtAllBins'][src][iEvent].Fill(jetIEta_offlineJet_tmp, res_dummy, puWeight)                                    
                                # if 'jet_byHand_res_vs_iEta_vs_nVtx' in dists2:
                                    # hist2['jet_byHand_res_vs_iEta_vs_nVtx%s' % (jetShape1)][algo1]['HBEF'        ][iPFJetPtCat][src][iEvent].Fill(jetIEta_offlineJet_tmp, nVtx, res_dummy, puWeight)
                                    # hist2['jet_byHand_res_vs_iEta_vs_nVtx%s' % (jetShape1)][algo1]['HBEF'        ]['PtAllBins'][src][iEvent].Fill(jetIEta_offlineJet_tmp, nVtx, res_dummy, puWeight)

                    hRefJet_pt_2p01.Fill(vOff.Pt())
                    hRefJet_eta_2p01.Fill(vOff.Eta())
                    hRefJet_phi_2p01.Fill(vOff.Phi())                                
                    
                    continue 

                hStat.Fill(32)
                hRefJet_pt_2p02.Fill(vOff.Pt())
                hRefJet_eta_2p02.Fill(vOff.Eta())
                hRefJet_phi_2p02.Fill(vOff.Phi()) 

                if not isMC and src == 'emu':
                    isMatchToUnp = False
                # Unpacked L1T jets in L1TNuples have jetEt information stored and not rawEt,  puEt etc, hence they can not be used for JEC derivation.
                # So for JEC derivation, use L1T emulated jets with the same jetEt, jetEta and jetPhi as of the unpacked L1T jet
                # Don't apply it when emulating with diffent layer 1 SF as that was used during data taking
                    for iUnp in range(nUnpJets):
                        if Unp_br['bx'][iEvent][iUnp] != 0: continue  ## Use only jets in BX 0

                        # Unp_br.jetEt[iUnp], Unp_br.jetEta[iUnp], Unp_br.jetPhi[iUnp]
                        # Emu_br.jetEt[iEmu], Emu_br.jetEta[iEmu], Emu_br.jetPhi[iEmu]
                        # l1jet_idx
                        # l1jet_br
                        if  ( (abs(Unp_br['pt'][iEvent][iUnp]  - l1jet_br['pt'][iEvent][l1jet_idx])  < 1e-8) and \
                                (abs(Unp_br['eta'][iEvent][iUnp] - l1jet_br['eta'][iEvent][l1jet_idx]) < 1e-8) and \
                                (abs(Unp_br['phi'][iEvent][iUnp] - l1jet_br['phi'][iEvent][l1jet_idx]) < 1e-8) ):
                            isMatchToUnp = True
                            break

                    if not isMatchToUnp and MatchEmulatedJetsWithUnpacked : continue

                hRefJet_pt_2p03.Fill(vOff.Pt())
                hRefJet_eta_2p03.Fill(vOff.Eta())
                hRefJet_phi_2p03.Fill(vOff.Phi())

                if src in ['emu']:
                    jetIEta_tmp    = convert_jetIEta_to_jetTowerIEta(l1jet_br['hwEta'][iEvent][l1jet_idx])
                    jetHwPt_tmp    = (l1jet_br['rawEt'][iEvent][l1jet_idx] - l1jet_br['puEt'][iEvent][l1jet_idx])
                    jetPt_tmp      = jetHwPt_tmp * 0.5 # 0.5 factor to conver hardware pT to GeV unit
                    jetIEtaBin_tmp = hJEC_iEta_vs_Pt.GetXaxis().FindBin( abs(jetIEta_tmp) )
                    jetPtBin_tmp   = hJEC_iEta_vs_Pt.GetYaxis().FindBin(jetPt_tmp)
                    JEC_tmp            =  float(2*l1jet_br['pt'][iEvent][l1jet_idx]) / float(jetHwPt_tmp)
                    hJEC_iEta_vs_Pt.SetBinContent(jetIEtaBin_tmp, jetPtBin_tmp, JEC_tmp)


                if l1jet_br['pt'][iEvent][l1jet_idx] > PT_MAX_L1Jet: continue
                hStat.Fill(34)

                hRefJet_pt_2.Fill(vOff.Pt())
                hRefJet_eta_2.Fill(vOff.Eta())
                hRefJet_phi_2.Fill(vOff.Phi())

                vL1Jet = R.TLorentzVector()
                vL1Jet.SetPtEtaPhiM(l1jet_br['pt'][iEvent][l1jet_idx], l1jet_br['eta'][iEvent][l1jet_idx], l1jet_br['phi'][iEvent][l1jet_idx], 0) # Aaron: use jet.etCorr instead of jet.et

                hdR_OffJet_L1Jet[src].Fill(vOff.DeltaR(vL1Jet))
                
                jetEtPUS_L1JetDefault = l1jet_br['pt'][iEvent][l1jet_idx]
                jetIEta = None
                jetIPhi = None
                if src == 'emu':
                    jetIEta = l1jet_br['towerIEta'][iEvent][l1jet_idx]
                    jetIPhi = l1jet_br['towerIPhi'][iEvent][l1jet_idx]
                elif src == 'unp':
                    jetIEta = l1jet_br['towerIEta'][iEvent][l1jet_idx]
                    jetIPhi = l1jet_br['towerIPhi'][iEvent][l1jet_idx]
                jetEt       = l1jet_br['pt'][iEvent][l1jet_idx]
                jetIEtaAbs  = abs(jetIEta)
                sjetIEta    = str(jetIEta)
                sjetIEtaAbs = str(jetIEtaAbs)
                sjetIEta_toUse = sjetIEtaAbs if useAbsEtaBins else sjetIEta;
                jetIEta_toUse  =  jetIEtaAbs if useAbsEtaBins else  jetIEta;
                sL1JetEtaCat = None
                for ieta_cat in IETA_CAT.keys():
                    if ieta_cat == 'HBEF': continue
                    if jetIEtaAbs >= IETA_CAT[ieta_cat][0] and jetIEtaAbs <= IETA_CAT[ieta_cat][1]:
                        sL1JetEtaCat = ieta_cat


                hist8['L1JetPtRawPUS'][src][sjetIEtaAbs].Fill( (l1jet_br['rawEt'][iEvent][l1jet_idx] - l1jet_br['puEt'][iEvent][l1jet_idx])*SFToConvertInGeV )
                hist8['L1JetPtRawPUS'][src]['HBEF'     ].Fill( (l1jet_br['rawEt'][iEvent][l1jet_idx] - l1jet_br['puEt'][iEvent][l1jet_idx])*SFToConvertInGeV )
                
                hist8['L1JetPtRawPUS_vs_RefJetPt'][src][sjetIEtaAbs].Fill(vOff.Pt(), (l1jet_br['rawEt'][iEvent][l1jet_idx] - l1jet_br['puEt'][iEvent][l1jet_idx])*SFToConvertInGeV )
                hist8['L1JetPtRawPUS_vs_RefJetPt'][src]['HBEF'     ].Fill(vOff.Pt(), (l1jet_br['rawEt'][iEvent][l1jet_idx] - l1jet_br['puEt'][iEvent][l1jet_idx])*SFToConvertInGeV )
                                    

                l1J_computationStarted = False # flag used to stoge L1J PU as histogram name
                
                data_dict['L1JetType']          = src
                data_dict['L1JetDefault_Et']    = jetEtPUS_L1JetDefault 
                data_dict['L1JetTowerIEtaAbs']  = jetIEtaAbs 
                #data_dict['L1JetTowerIEta'] = jetIEta
                #data_dict['L1JetTowerIPhi'] = jetIPhi

                TTiet_max_JetShapewise = {}
                for jetShape in JetShapes:
                    if jetShape == 'Default': continue
                    TTiet_max_JetShapewise[jetShape] = -1

                ## Compare pT to sum from ECAL and HCAL TPs
                Raw_HTP_ET  = 0
                Raw_HTP_iET = 0
                Raw_ETP_ET  = 0
                Raw_ETP_iET = 0
                Raw_TT_iET  = 0
                Raw_TT_nET  = 0
                PUS_TT_iET  = [0,0,0,0]
                PUS_TT_ring = [0,0,0,0,0] # phi ring excluding central and adjacant clusters
                PUS_TT_ring2 = [0,0,0,0,0,0,0] # phi ring excluding central cluster
                
                # Different shape jets
                Raw_TT_iET_ByJetShape   = {}
                PUS_TT_ring_ByJetShape  = {} # phi ring excluding central and adjacant clusters
                PUS_TT_ring2_ByJetShape = {} # phi ring excluding central cluster
                for jetShape in JetShapes:
                    if jetShape == 'Default': continue
                    Raw_TT_iET_ByJetShape[jetShape]   = 0
                    PUS_TT_ring_ByJetShape[jetShape]  = [0,0,0,0,0] # phi ring excluding central and adjacant clusters
                    PUS_TT_ring2_ByJetShape[jetShape] = [0,0,0,0,0,0,0] # phi ring excluding central cluster
                
                for iTT in range(l1TT_br["nObjects"][iEvent]):
                    TTieta = l1TT_br['ieta'][iEvent][iTT]
                    TTiphi = l1TT_br['iphi'][iEvent][iTT]
                    TTiet  = l1TT_br['iet'][iEvent][iTT] * 0.5 # multiply each "TT.iet" value by 0.5, to convert to units of GeV
                    #TTiet  = l1TT_br.iet [iEvent][iTT] # multiply each "TT.iet" value by 0.5, to convert to units of GeV <<<<<<<<<<<<< WRONG ?? <<<<<<<<<<<<<<<<<<<<<<<
                                                

                    if hCaloTowers_iEta_vs_iPhi:
                        TTieta_bin = hCaloTowers_iEta_vs_iPhi.GetXaxis().FindBin( TTieta )
                        TTiphi_bin = hCaloTowers_iEta_vs_iPhi.GetYaxis().FindBin( TTiphi )
                        hCaloTowers_iEta_vs_iPhi.SetBinContent(TTieta_bin, TTiphi_bin, TTiet)

                    
                    dIEta_TT_Seed    = dIEta(TTieta, jetIEta)
                    dIPhi_TT_Seed    = dIPhi(TTiphi, jetIPhi)
                    AbsdIEta_TT_Seed = abs( dIEta_TT_Seed )
                    AbsdIPhi_TT_Seed = abs( dIPhi_TT_Seed )
                    
                    if PrintLevel >= 7:
                        print("    %s %s TTieta %g, TTiphi %g, TTiet %g, dIEta_TT_Seed %g, dIPhi_TT_Seed %g" % \
                            (" " * 6, src,
                                TTieta,TTiphi,TTiet, dIEta_TT_Seed,dIPhi_TT_Seed ))
                    
                    # variable shape jets with phi ring PU subtraction                       
                    for jetShape in JetShapes:
                        if jetShape == 'Default': continue
                        # jetShape:
                        # '3x9':                    jet Et = sum Et in 3x9 region aroud seed
                        # '3x9_plus_0.5_times_9x9': jet Et = (sum Et in 3x9 aroud seed) + 0.5*(sum Et in 9x9 aroud seed, but outside 3x9 central part)
                        outerClusterEtWgt = None # this is also used as a flag if cluster has central and outer parts with lower weights to outer parts
                        
                        # cluster's central part. For e.g. for jetShape = '3x9_plus_0.5_times_9x9', it read 3x9 
                        centralClustSizeInEta = int(jetShape.split('_')[0].split('x')[0])
                        centralClustSizeInPhi = int(jetShape.split('_')[0].split('x')[1]) 
                        # cluster's outer part, if that is the case.  For e.g. for jetShape = '3x9_plus_0.5_times_9x9', it read 9x9 
                        outerClustSizeInEta   = int(jetShape.split('_')[-1].split('x')[0])
                        outerClustSizeInPhi   = int(jetShape.split('_')[-1].split('x')[1])
                        if '_plus_' in jetShape:
                            outerClusterEtWgt = float(jetShape.split('_plus_')[1].split('_times_')[0])
                        
                        # consider full cluster size. i.e. including outer cluster if that is the case
                        # for case e.g. jetShape = '3x9'   or   central part for '3x9_plus_0.5_times_9x9'
                        clustSizeInEta        = centralClustSizeInEta
                        clustSizeInPhi        = centralClustSizeInPhi
                        if outerClusterEtWgt: # cluster's outer part with lower Et weights to consider. for case e.g. jetShape = '3x9_plus_0.5_times_9x9'
                            clustSizeInEta    = outerClustSizeInEta
                            clustSizeInPhi    = outerClustSizeInPhi
                        
                        clustSizeAroundSeedInPhi                                   = clustSizeInPhi / 2
                        clustSizeAroundSeedInEtaLow = clustSizeAroundSeedInEtaHigh = clustSizeInEta / 2
                        if (clustSizeInEta % 2) == 0: # even TT size cluster
                            if jetIEta > 0:
                                clustSizeAroundSeedInEtaLow  = (clustSizeInEta / 2)
                                clustSizeAroundSeedInEtaHigh = (clustSizeInEta / 2) - 1
                            else:
                                clustSizeAroundSeedInEtaLow  = (clustSizeInEta / 2) - 1
                                clustSizeAroundSeedInEtaHigh = (clustSizeInEta / 2)
                        
                        centralClustSizeAroundSeedInEtaLow = centralClustSizeAroundSeedInEtaHigh = centralClustSizeInEta / 2
                        if outerClusterEtWgt: # cluster's outer part with lower Et weights to consider:
                            if (centralClustSizeInEta % 2) == 0: # even TT size cluster
                                if jetIEta > 0:
                                    centralClustSizeAroundSeedInEtaLow  = (centralClustSizeInEta / 2)
                                    centralClustSizeAroundSeedInEtaHigh = (centralClustSizeInEta / 2) - 1
                                else:
                                    centralClustSizeAroundSeedInEtaLow  = (centralClustSizeInEta / 2) - 1
                                    centralClustSizeAroundSeedInEtaHigh = (centralClustSizeInEta / 2)


                        
                        if ( dIEta_TT_Seed >= -1*clustSizeAroundSeedInEtaLow and dIEta_TT_Seed <= clustSizeAroundSeedInEtaHigh ): # within cluster phi ring
                            TTiet_toUse = TTiet
                            if outerClusterEtWgt: # cluster's outer part with lower Et weights to consider:
                                if not ( dIEta_TT_Seed >= -1*centralClustSizeAroundSeedInEtaLow and dIEta_TT_Seed <= centralClustSizeAroundSeedInEtaHigh ): # outside centralCluster phi ring
                                    TTiet_toUse = TTiet * outerClusterEtWgt
                            
                            ## "Phi ring" PU subtraction sides
                            ring_dPhi     = dIPhi_TT_Seed 
                            ring_dPhi_pos = ring_dPhi
                            if ring_dPhi < 0:
                                ring_dPhi_pos = 72 - abs(ring_dPhi)  ## Shift into all-positive scale

                            if abs(ring_dPhi) > (4 + 9):  ## Not adjacent in phi
                                PUS_TT_ring_ByJetShape[jetShape][int((ring_dPhi_pos - 14) / 9)] += TTiet_toUse  ## Fill 5 9x9 regions

                            if abs(ring_dPhi) > 4: ## ring starting from adjacent phi cluster
                                #print(f"Siddh here1: jetShape {jetShape},  ring_dPhi_pos {ring_dPhi_pos}, (ring_dPhi_pos - 5) / 9: {(ring_dPhi_pos - 5) / 9} TTiet_toUse: {TTiet_toUse}"); sys.stdout.flush();
                                PUS_TT_ring2_ByJetShape[jetShape][int((ring_dPhi_pos - 5) / 9)] += TTiet_toUse  ## Fill 7 9x9 regions
                                if VERBOSE or PrintLevel >= 5:
                                    print("    %s %s: ieta %d, iphi %d, iet %g, iet_toUse %g,  dIEta_TT_Seed %g, dIPhi_TT_Seed %g -- TT in PhiRing PU" % \
                                    (" " * 8, jetShape, TTieta,TTiphi,TTiet,TTiet_toUse,
                                        dIEta_TT_Seed,dIPhi_TT_Seed )); sys.stdout.flush();

                                if jetShape == '9x9':
                                    if hCaloTTs_iEta_vs_iPhi:
                                        TTieta_bin = hCaloTTs_iEta_vs_iPhi.GetXaxis().FindBin( TTieta )
                                        TTiphi_bin = hCaloTTs_iEta_vs_iPhi.GetYaxis().FindBin( TTiphi )
                                        hCaloTTs_iEta_vs_iPhi.SetBinContent(TTieta_bin, TTiphi_bin, TTiet)

                                if TTiet > TTiet_max_JetShapewise[jetShape]: TTiet_max_JetShapewise[jetShape] = TTiet


                            if ( AbsdIPhi_TT_Seed <= clustSizeAroundSeedInPhi ): # within cluster
                                Raw_TT_iET_ByJetShape[jetShape] += TTiet_toUse
                                if VERBOSE or PrintLevel >= 5:
                                    print("    %s %s: ieta %d, iphi %d, iet %g, iet_toUse %g,  dIEta_TT_Seed %g, dIPhi_TT_Seed %g -- TT within cluster" % \
                                    (" " * 8, jetShape, TTieta,TTiphi,TTiet,TTiet_toUse,
                                        dIEta_TT_Seed,dIPhi_TT_Seed ))

                                if jetShape == '9x9':
                                    if hCaloTTs_iEta_vs_iPhi:
                                        TTieta_bin = hCaloTTs_iEta_vs_iPhi.GetXaxis().FindBin( TTieta )
                                        TTiphi_bin = hCaloTTs_iEta_vs_iPhi.GetYaxis().FindBin( TTiphi )
                                        hCaloTTs_iEta_vs_iPhi.SetBinContent(TTieta_bin, TTiphi_bin, TTiet)
                                
                                if TTiet > TTiet_max_JetShapewise[jetShape]: TTiet_max_JetShapewise[jetShape] = TTiet
                                        
                    
                    ## "Phi ring" PU subtraction sides
                    if abs( dIEta(TTieta, jetIEta) ) <= 4:  ## In the same 9-ring region as the jet
                        ring_dPhi     = dIPhi(TTiphi, jetIPhi)
                        ring_dPhi_pos = ring_dPhi
                        if ring_dPhi < 0:
                            ring_dPhi_pos = 72 - abs(ring_dPhi)  ## Shift into all-positive scale

                        if abs(ring_dPhi) > (4 + 9):  ## Not adjacent in phi
                            PUS_TT_ring[int((ring_dPhi_pos - 14) / 9)] += TTiet  ## Fill 5 9x9 regions
                            
                        if abs(ring_dPhi) > 4: ## ring starting from adjacent phi cluster
                            PUS_TT_ring2[int((ring_dPhi_pos - 5) / 9)] += TTiet  ## Fill 7 9x9 regions

                    ## "Chunky doughnut" PU subtraction sides
                    if abs( dIEta(TTieta, jetIEta) ) > 7 or  abs( dIPhi(TTiphi, jetIPhi) ) > 7: continue
                    if abs( dIEta(TTieta, jetIEta) ) > 4 and abs( dIPhi(TTiphi, jetIPhi) ) > 4: continue

                    if dIEta(TTieta, jetIEta) >  4 and dIEta(TTieta, jetIEta) <  8 and abs( dIPhi(TTiphi, jetIPhi) ) < 5:
                        PUS_TT_iET[0] += TTiet
                        # print '      ** Added to PUS[0]'
                    if dIEta(TTieta, jetIEta) < -4 and dIEta(TTieta, jetIEta) > -8 and abs( dIPhi(TTiphi, jetIPhi) ) < 5:
                        PUS_TT_iET[1] += TTiet
                        # print '      ** Added to PUS[1]'
                    if dIPhi(TTiphi, jetIPhi) >  4 and dIPhi(TTiphi, jetIPhi) <  8 and abs( dIEta(TTieta, jetIEta) ) < 5:
                        PUS_TT_iET[2] += TTiet
                        # print '      ** Added to PUS[2]'
                    if dIPhi(TTiphi, jetIPhi) < -4 and dIPhi(TTiphi, jetIPhi) > -8 and abs( dIEta(TTieta, jetIEta) ) < 5:
                        PUS_TT_iET[3] += TTiet
                        # print '      ** Added to PUS[3]'

                    ## Central jet sum
                    if abs( dIEta(TTieta, jetIEta) ) > 4 or abs( dIPhi(TTiphi, jetIPhi) ) > 4: continue
                    Raw_TT_iET += TTiet
                    Raw_TT_nET += 1
                    # print '    - Trigger tower iEta = %d, iPhi = %d, pT = %.1f' % (TTieta, TTiphi, TTiet)
                    # print '      dIEta(%d, %d) = %d' % (TTieta, jetIEta, dIEta(TTieta, jetIEta))
                    # print '      ** Added to central sum'
                    if VERBOSE:
                        print("    %sieta = %d, iphi = %d, iet = %g, iem = %g, ihad = %g, iratio = %g, iqual = %g, et = %g, eta = %g, phi = %g" % (" " * 8,eTT_br['ieta'][iEvent][iTT], eTT_br['iphi'][iEvent][iTT], eTT_br['iet'][iEvent][iTT], eTT_br.iem[iEvent][iTT], eTT_br.ihad[iEvent][iTT], eTT_br.iratio[iEvent][iTT], eTT_br.iqual[iEvent][iTT], eTT_br.et[iEvent][iTT], eTT_br['eta'][iEvent][iTT], eTT_br['phi'][iEvent][iTT]))

                PUet_ChunkyDonut     =  sum(PUS_TT_iET)  - max(PUS_TT_iET)
                PUS_TT_ring_Min4     = (sum(PUS_TT_ring) - max(PUS_TT_ring)) / 4.0
                PUS_TT_ring_Side4    = (sum(PUS_TT_ring) - PUS_TT_ring[2])   / 4.0
                PUS_TT_ring_Adjacent = (PUS_TT_ring2[0] + PUS_TT_ring2[-1]) / 2.0
                PUS_TT_ring_Full     = sum(PUS_TT_ring2) + Raw_TT_iET
                
                if VERBOSE or PrintLevel >= 6:
                    #print "%8s       Raw_TT_iET: %g, PUet_ChunkyDonut: %g, PUS_TT_ring_Min4: %g, PUS_TT_ring_Side4: %g, PUS_TT_ring_Adjacent: %g" % (" ", Raw_TT_iET, PUet_ChunkyDonut, PUS_TT_ring_Min4, PUS_TT_ring_Side4, PUS_TT_ring_Adjacent)
                    print("%8s   Default    Raw_TT_iET: %g, PUet_ChunkyDonut: %g" % (" ", Raw_TT_iET, PUet_ChunkyDonut))
                    for jetShape in JetShapes:
                        if jetShape == 'Default': continue
                        PUS_TT_ring_Min4_tmp     = (sum(PUS_TT_ring_ByJetShape[jetShape]) - max(PUS_TT_ring_ByJetShape[jetShape])) / 4.0
                        PUS_TT_ring_Side4_tmp    = (sum(PUS_TT_ring_ByJetShape[jetShape]) - PUS_TT_ring_ByJetShape[jetShape][2])   / 4.0
                        PUS_TT_ring_Adjacent_tmp = (PUS_TT_ring2_ByJetShape[jetShape][0] + PUS_TT_ring2_ByJetShape[jetShape][-1]) / 2.0                            
                        #print "%8s %s: Raw_TT_iET: %g, PUet_ChunkyDonut: %g, PUS_TT_ring_Min4: %g, PUS_TT_ring_Side4: %g, PUS_TT_ring_Adjacent: %g" % (" ", jetShape, Raw_TT_iET_ByJetShape[jetShape], -1, PUS_TT_ring_Min4_tmp, PUS_TT_ring_Side4_tmp, PUS_TT_ring_Adjacent_tmp)

                        Raw_TT_iET_tmp           = Raw_TT_iET_ByJetShape[jetShape]
                        PUet_ChunkyDonut_tmp                 = PUet_ChunkyDonut if jetShape == '9x9' else 0
                        PUS_TT_ring_Full_tmp     = sum(PUS_TT_ring2_ByJetShape[jetShape]) + Raw_TT_iET_ByJetShape[jetShape]
                        
                        l1jet_PU_pt_PhiRing = Raw_TT_iET_tmp - (8.0/7*Raw_TT_iET_tmp - 1./7*PUS_TT_ring_Full_tmp) #  l1jet_pt = (8.0/7*JetRaw - 1./7*SumFullPhiRing)
                        print("%8s   %s   RAWPUS: Raw_TT_iET: %g, PUet_ChunkyDonut: %g, PUS: %g \t PhiRIng: PU: %g, PUS: %g " % \
                            (" ", JetShapes, \
                                Raw_TT_iET_tmp, PUet_ChunkyDonut_tmp, ( Raw_TT_iET_tmp - PUet_ChunkyDonut_tmp), \
                                l1jet_PU_pt_PhiRing, (Raw_TT_iET_tmp - l1jet_PU_pt_PhiRing)
                                ))


                if l1jet_br['puEt'][iEvent][l1jet_idx] == 0 and l1jet_br['rawEt'][iEvent][l1jet_idx] > 0 and runMode in ['trbshtPhiRingPUS'] and \
                    '9x9' in TTiet_max_JetShapewise.keys():
                    hTTEtMax_forL1JetPUEt0.Fill( TTiet_max_JetShapewise['9x9'] )
                    hL1JetRawEt_vs_L1JetEt_forL1JetPUEt0.Fill(l1jet_br['rawEt'][iEvent][l1jet_idx], l1jet_br['pt'][iEvent][l1jet_idx])
                    print("puEt {}, rawEt {}, jetSeedEt {}, jetEt {} \t TTiet_max {}".format(l1jet_br['puEt'][iEvent][l1jet_idx], l1jet_br['rawEt'][iEvent][l1jet_idx], l1jet_br.seedEt[iEvent][l1jet_idx], l1jet_br['pt'][iEvent][l1jet_idx],  TTiet_max_JetShapewise['9x9']))
                    
                #print "jetIEta {}, src {}, iEvent {}".format(jetIEta, src, iEvent)
                if 'l1jetEt_vs_RawEtMinusPU' in dists1:
                    hist1['l1jetEt_vs_RawEtMinusPU'][sjetIEta_toUse][src][iEvent].Fill(jetEt, (Raw_TT_iET - PUet_ChunkyDonut) / jetEt, puWeight)
                    hist1['l1jetEt_vs_RawEtMinusPU']['HBEF'        ][src][iEvent].Fill(jetEt, (Raw_TT_iET - PUet_ChunkyDonut) / jetEt, puWeight)

                
                
                
                l1jet_pt_JetShapeAndAlgoWise = {}
                # jet resolution plots for different jet shapes   
                #for jetShape in ['Default'] + JetShapes:
                for jetShape in JetShapes:
                    # JetShape = "" plots are with the first version of code for 9x9 jets
                    jetShape1 = jetShape
                    if jetShape == 'Default':  jetShape1 = ""
                    else:                      jetShape1 = "_%s" % (jetShape)
                    
                    # assign cluster and PU energy for jetShape under consideration
                    Raw_TT_iET_tmp           = 0
                    PUet_tmp                 = 0
                    PUS_TT_ring_Min4_tmp     = 0
                    PUS_TT_ring_Side4_tmp    = 0
                    PUS_TT_ring_Adjacent_tmp = 0
                    PUS_TT_ring_Full_tmp     = 0
                    if jetShape == 'Default':
                        #Raw_TT_iET_tmp           = Raw_TT_iET
                        #PUet_tmp                 = PUet_ChunkyDonut
                        Raw_TT_iET_tmp           = l1jet_br['rawEt'][iEvent][l1jet_idx]
                        PUet_tmp                 = l1jet_br['puEt'][iEvent][l1jet_idx]
                        PUS_TT_ring_Min4_tmp     = PUS_TT_ring_Min4
                        PUS_TT_ring_Side4_tmp    = PUS_TT_ring_Side4
                        PUS_TT_ring_Adjacent_tmp = PUS_TT_ring_Adjacent
                        PUS_TT_ring_Full_tmp     = PUS_TT_ring_Full
                    else:
                        Raw_TT_iET_tmp           = Raw_TT_iET_ByJetShape[jetShape]
                        PUet_tmp                 = PUet_ChunkyDonut if jetShape == '9x9' else 0
                        PUS_TT_ring_Min4_tmp     = (sum(PUS_TT_ring_ByJetShape[jetShape]) - max(PUS_TT_ring_ByJetShape[jetShape])) / 4.0 
                        PUS_TT_ring_Side4_tmp    = (sum(PUS_TT_ring_ByJetShape[jetShape]) - PUS_TT_ring_ByJetShape[jetShape][2])   / 4.0 
                        PUS_TT_ring_Adjacent_tmp = (PUS_TT_ring2_ByJetShape[jetShape][0] + PUS_TT_ring2_ByJetShape[jetShape][-1]) / 2.0
                        PUS_TT_ring_Full_tmp     = sum(PUS_TT_ring2_ByJetShape[jetShape]) + Raw_TT_iET_ByJetShape[jetShape]      
                    
                    if jetShape == 'Default':
                        #data_dict['L1JetDefault_RawEtPUS']                  = Raw_TT_iET_tmp - PUet_tmp
                        #data_dict['L1JetDefault_PU']                     = PUet_tmp
                        data_dict['L1JetDefault_RawEt']                  = l1jet_br['rawEt'][iEvent][l1jet_idx] * 0.5 # multiply  0.5, to convert to units of GeV
                        if   l1nanoChunkyDonut:
                            data_dict['L1JetDefault_PUEt_ChunkyDonut']   = l1jet_br['puEt'][iEvent][l1jet_idx] * 0.5 # multiply  0.5, to convert to units of GeV
                        elif l1nanoPhiRing:
                            data_dict['L1JetDefault_PUEt_PhiRing']       = l1jet_br['puEt'][iEvent][l1jet_idx] * 0.5 # multiply  0.5, to convert to units of GeV
                            
                    elif jetShape in JetShapesForML: # for machine learning training 
                        data_dict['L1Jet%s_RawEt' % (jetShape)]          = Raw_TT_iET_ByJetShape[jetShape]
                        data_dict['L1Jet%s_EtSum7PUTowers' % (jetShape)] = sum(PUS_TT_ring2_ByJetShape[jetShape])
                        if jetShape == '9x9':
                            data_dict['L1Jet%s_PUEt_ChunkyDonut' % (jetShape)]   = PUet_ChunkyDonut 
                        

                    l1jet_pt_JetShapeAndAlgoWise[jetShape] = {}  

                    for algo1 in PUSAlgosAll: # ['Raw', 'RawPUS', 'RawPUS_phiRingMin4', 'RawPUS_phiRingSide4', 'RawPUS_phiRingAdjacent']: # ['PUS','noPUS','Raw','RawPUS']:
                        if (algo1 == 'L1JDefault') and (jetShape != 'Default'): continue
                        if (not l1nanoChunkyDonut) and (jetShape == 'Default') and (algo1 == 'RawPUS'):            continue
                        if (not l1nanoPhiRing)     and (jetShape == 'Default') and (algo1 == 'RawPUS_phiDefault'): continue
                        
                        if Raw_TT_iET_tmp <= 0: continue # skip cluster with eT=0

                        l1jet_PU_pt = 0
                        if jetShape == 'Default':
                            if algo1  not in ['L1JDefault', 'Raw', 'RawPUS',  'RawPUS_phiDefault']: continue
                            
                            if algo1 == 'Raw':                     l1jet_PU_pt = 0.
                            
                            if   l1nanoChunkyDonut:
                                if algo1 == 'RawPUS':                  l1jet_PU_pt = PUet_tmp
                                if algo1 == 'RawPUS_phiDefault':       continue
                            elif l1nanoPhiRing:
                                if algo1 == 'RawPUS':                  continue
                                if algo1 == 'RawPUS_phiDefault':       l1jet_PU_pt = PUet_tmp
                                    
                        else:
                            if algo1 == 'Raw':                     l1jet_PU_pt = 0.
                            if algo1 == 'RawPUS':                  l1jet_PU_pt = PUet_tmp
                            if algo1 == 'RawPUS_phiRingMin4':      l1jet_PU_pt = PUS_TT_ring_Min4_tmp
                            if algo1 == 'RawPUS_phiRingSide4':     l1jet_PU_pt = PUS_TT_ring_Side4_tmp
                            if algo1 == 'RawPUS_phiRingAdjacent':  l1jet_PU_pt = PUS_TT_ring_Adjacent_tmp
                            if algo1 == 'RawPUS_phiDefault':       l1jet_PU_pt = Raw_TT_iET_tmp - (8.0/7*Raw_TT_iET_tmp - 1./7*PUS_TT_ring_Full_tmp) #  l1jet_pt = (8.0/7*JetRaw - 1./7*SumFullPhiRing)                            
                        l1jet_pt = Raw_TT_iET_tmp - l1jet_PU_pt
                        l1jet_pt_woLayer2Calib = l1jet_pt
                        


                        sPrintTmp1 = ""
                        if runMode in ['CalibJetByHand']:
                            sPrintTmp1 += "%4s%8s, %24s l1j pT = %g " % (' ', jetShape, algo1, l1jet_pt)


                        # calibSF_additionalCorr -------------
                        if runMode in ['CalibJetByHand'] and \
                            jetShape in calibSFs_additionalCorr.keys() and algo1 in calibSFs_additionalCorr[jetShape].keys():
                            calibSF_additionalCorr_tmp = 1.

                            # handle jetShape == 'Default' and algo1 == 'Raw' case first, whiEvent has seprate additional SFs for l1nanoChunkyDonut and l1nanoPhiRing
                            if jetShape == 'Default' and algo1 == 'Raw':
                                sL1nanoPUS_tmp = 'l1nanoChunkyDonut' if l1nanoChunkyDonut else 'l1nanoPhiRing'
                                if sL1nanoPUS_tmp not in calibSFs_additionalCorr[jetShape][algo1].keys():
                                    print("sL1nanoPUS_tmp {} not in calibSFs_additionalCorr[{}][{}]: {} \t\t\t **** ERROR ****".format(sL1nanoPUS_tmp, jetShape,algo1, calibSFs_additionalCorr[jetShape][algo1].keys()))
                                else:
                                    calibSF_additionalCorr_tmp = calibSFs_additionalCorr[jetShape][algo1][sL1nanoPUS_tmp]
                            else:
                                calibSF_additionalCorr_tmp = calibSFs_additionalCorr[jetShape][algo1]

                            l1jet_pt *= calibSF_additionalCorr_tmp
                            l1jet_pt_woLayer2Calib *= calibSF_additionalCorr_tmp

                            sPrintTmp1 += " * %g = %g " % (calibSF_additionalCorr_tmp, l1jet_pt)


                            
                        # l1jet_pt calibration: layer2CalibSF ----------                            
                        if runMode in ['CalibJetByHand'] and \
                            jetShape in calibSFs.keys() and algo1 in calibSFs[jetShape].keys() and sjetIEta_toUse in calibSFs[jetShape][algo1].keys():
                            for layer2CalibSF_list in calibSFs[jetShape][algo1][sjetIEta_toUse]:
                                # layer2CalibSF_list: [bin_pt_low, bin_pt_up, layer2CalibSF]                                    
                                if not (l1jet_pt_woLayer2Calib >= layer2CalibSF_list[0] and l1jet_pt_woLayer2Calib < layer2CalibSF_list[1] ): continue
                                layer2CalibSF_tmp = layer2CalibSF_list[2]
                                
                                l1jet_pt = l1jet_pt_woLayer2Calib * layer2CalibSF_tmp

                                sPrintTmp1 += " * %g = %g " % (layer2CalibSF_tmp, l1jet_pt)

                                

                            

                        if PrintLevel >= 1 and runMode in ['CalibJetByHand']:
                            print(sPrintTmp1)


                        

                        if PrintLevel >= 3 :
                            sTmp = ""                                
                            if algo1 == 'RawPUS_phiDefault' and l1jet_br['puEt'][iEvent][l1jet_idx] > 0:
                                sTmp = "\t PUEt myCal/default: (%g / %g) %.2f" % (PUS_TT_ring_Full_tmp/ 8, l1jet_br['puEt'][iEvent][l1jet_idx], PUS_TT_ring_Full_tmp/ 8 /l1jet_br['puEt'][iEvent][l1jet_idx])
                            print("%4s%8s, %24s, %3s, pT: Raw: %7.2f, PU: %7.2f, l1j: %7.2f,  %7.2f,    diff: %g %s" % \
                                (' ',jetShape,algo1,sjetIEta_toUse,  \
                                    Raw_TT_iET_tmp, l1jet_PU_pt, l1jet_pt_woLayer2Calib, l1jet_pt, (l1jet_pt - l1jet_pt_woLayer2Calib), sTmp) )
                        

                        # troubleshoot histograms
                        if (jetShape == '9x9') and (algo1 == 'RawPUS_phiDefault') and l1jet_br['puEt'][iEvent][l1jet_idx] > 0 and 1==0:
                            #hists7['CheckPUEt_iEta_vs_PUEtRatioWrtDefault_9x9_RawPUS_phiDefault'][src][iEvent].Fill(jetIEta_toUse, l1jet_PU_pt / l1jet_br['puEt'][iEvent][l1jet_idx])
                            hists7['CheckPUEt_iEta_vs_PUEtRatioWrtDefault_9x9_RawPUS_phiDefault'][src][iEvent].Fill(jetIEta_toUse, float(int(PUS_TT_ring_Full_tmp/8.0)) / l1jet_br['puEt'][iEvent][l1jet_idx])
                            #if  l1jet_PU_pt / l1jet_br['puEt'][iEvent][l1jet_idx] < 0.10: # and :
                            if  abs( (PUS_TT_ring_Full_tmp/ 8 / l1jet_br['puEt'][iEvent][l1jet_idx] ) - 1 ) > 0.01: # and :
                                selectedEvents_list.append( "%d:%d:%d" % (int(Evt_br['run']), int(Evt_br['luminosityBlock']), int(Evt_br['event'])) )

                            if hCaloTTs_iEta_vs_iPhi:
                                hCaloTTs_iEta_vs_iPhi.SetTitle( "%s, %.1f/%g" % (hCaloTTs_iEta_vs_iPhi.GetTitle(), l1jet_PU_pt,l1jet_br['puEt'][iEvent][l1jet_idx]) )


                        l1jet_pt_JetShapeAndAlgoWise[jetShape][algo1] = l1jet_pt
                        # troubleshoot l1jet_pt_phiRing
                        if 1==0 and jetShape == '9x9' and algo1 == 'RawPUS_phiDefault' and abs(l1jet_pt - l1jet_pt_JetShapeAndAlgoWise['Default']['RawPUS_phiDefault']) > 2:
                            selectedEvents_list.append( "%d:%d:%d" % (int(Evt_br['run']), int(Evt_br['luminosityBlock']), int(Evt_br['event'])) )
                            print("l1jet_pt_JetShapeAndAlgoWise['Default']['RawPUS_phiDefault'] {}, l1jet_pt_JetShapeAndAlgoWise['9x9']['RawPUS_phiDefault'] {} \t\t\t <<<<< CHECK THIS ".format( \
                                                                                                                                                                                                    l1jet_pt_JetShapeAndAlgoWise['Default']['RawPUS_phiDefault'], l1jet_pt_JetShapeAndAlgoWise['9x9']['RawPUS_phiDefault'])
                                    )
                            
                            
                                
                        # L1JetDefault -------------------------------------------------------
                        if (algo1 == 'L1JDefault') and (jetShape == 'Default'):
                            l1jet_pt = jetEtPUS_L1JetDefault
                        # --------------------------------------------------------------------
                        
                        iL1JetPtCat = getJetPtCategory( l1jet_pt )
                        iL1JetPtCat_woLayer2Calib = getJetPtCategory( l1jet_pt_woLayer2Calib )    
                            
                        '''
                        print "jetShape1: {}, algo1 {}, str(jetIEta) {}, jPt {}, src {}, iEvent {}".format('jet_byHand_den%s' % (jetShape1), algo1, str(jetIEta), jPt, src, iEvent)
                        if 'jet_byHand_den%s' % (jetShape1) not in hist2:
                            print "  {} doesn't exist".format('jet_byHand_den%s' % (jetShape1))
                        if algo1 not in hist2['jet_byHand_den%s' % (jetShape1)]:
                            print "  {} doesn't exist".format(algo1)
                        if str(jetIEta) not in hist2['jet_byHand_den%s' % (jetShape1)][algo1]:
                            print "  {} doesn't exist".format(str(jetIEta))
                        if 'HBEF' not in hist2['jet_byHand_den%s' % (jetShape1)][algo1]:
                            print "  {} doesn't exist".format()
                        if jPt not in hist2['jet_byHand_den%s' % (jetShape1)][algo1][str(jetIEta)]:
                            print "  {} doesn't exist".format(jPt)
                        if src not in hist2['jet_byHand_den%s' % (jetShape1)][algo1][str(jetIEta)][jPt]:
                            print "  {} doesn't exist".format(src)
                        if iEvent not in hist2['jet_byHand_den%s' % (jetShape1)][algo1][str(jetIEta)][jPt][src]:
                            print "  {} doesn't exist".format(iEvent)
                        print "hist2['jet_byHand_den%s' % (jetShape1)][algo1][str(jetIEta)][jPt][src]: {}".format(hist2['jet_byHand_den%s' % (jetShape1)][algo1][str(jetIEta)][jPt][src])
                        '''
                        
                        #l1JetCollection[src][jetShape][algo1].append( [l1jet_pt, sL1JetEtaCat] )
                        #l1JetCollection[src][jetShape][algo1][sL1JetEtaCat].append( [l1jet_pt, sL1JetEtaCat] ) 
                        l1JetCollection[src][jetShape][algo1][sEtaCat_PFJet].append( [l1jet_pt, sEtaCat_PFJet] )  # use EtaCat of PF jet for rate plots
                            
                        for jPt in PT_CAT.keys():
                            if 'jet_byHand_den' in dists2:
                                hist2['jet_byHand_den%s' % (jetShape1)][algo1][sjetIEta_toUse][jPt][src][iEvent].Fill( vOff.Pt() )
                                hist2['jet_byHand_den%s' % (jetShape1)][algo1]['HBEF'        ][jPt][src][iEvent].Fill( vOff.Pt() )

                            if 'jet_byHand_eff_den_vs_PU' in dists3:    
                                hist3['jet_byHand_eff_den_vs_PU%s' % (jetShape1)][algo1][sEtaCat_PFJet][jPt][src][0].Fill( vOff.Pt(), nVtx, puWeight )       #revisit
                                hist3['jet_byHand_eff_den_vs_PU%s' % (jetShape1)][algo1]['HBEF'       ][jPt][src][0].Fill( vOff.Pt(), nVtx, puWeight )
                            
                            if l1jet_pt > PT_CAT[jPt][1]:
                                if 'jet_byHand_num' in dists2:
                                    hist2['jet_byHand_num%s' % (jetShape1)][algo1][sjetIEta_toUse][jPt][src][iEvent].Fill( vOff.Pt() )
                                    hist2['jet_byHand_num%s' % (jetShape1)][algo1]['HBEF'        ][jPt][src][iEvent].Fill( vOff.Pt() )

                                # if 'jet_byHand_eff_num_vs_PU' in dists3:
                                #     hist3['jet_byHand_eff_num_vs_PU%s' % (jetShape1)][algo1][sEtaCat_PFJet][jPt][src][iEvent].Fill( vOff.Pt(), nVtx, puWeight )
                                #     hist3['jet_byHand_eff_num_vs_PU%s' % (jetShape1)][algo1]['HBEF'       ][jPt][src][iEvent].Fill( vOff.Pt(), nVtx, puWeight )

                        
                                
                        res = (l1jet_pt - vOff.Pt()) / vOff.Pt()
                        res_woLayer2Calib = (l1jet_pt_woLayer2Calib - vOff.Pt()) / vOff.Pt()

                        if PrintLevel >= 3 :
                            print("%4s%8s, %24s: vOff %7.2f, l1j: %7.2f -> %7.2f,  res: %4.2f %4.2f" % \
                                    (' ',jetShape,algo1, vOff.Pt(), l1jet_pt_woLayer2Calib, l1jet_pt, res, res_woLayer2Calib)
                                    )

                            
                        #hist2['jet_byHand_res_vs_iEta%s' % (jetShape1)][algo1]['HBEF'        ][iPFJetPtCat        ][src][iEvent].Fill(jetIEta_toUse, res, puWeight)
                        #hist2['jet_byHand_res_vs_iEta%s' % (jetShape1)][algo1]['HBEF'        ]['PtAllBins'][src][iEvent].Fill(jetIEta_toUse, res, puWeight)

                        #hist2['jet_byHand_PU_vs_iEta%s' % (jetShape1)][algo1]['HBEF'        ][iPFJetPtCat        ][src][iEvent].Fill(jetIEta_toUse, l1jet_PU_pt, puWeight)
                        #hist2['jet_byHand_PU_vs_iEta%s' % (jetShape1)][algo1]['HBEF'        ]['PtAllBins'][src][iEvent].Fill(jetIEta_toUse, l1jet_PU_pt, puWeight)

                        #hist2['jet_byHand_PUByRawPt_vs_iEta%s' % (jetShape1)][algo1]['HBEF'        ][iPFJetPtCat        ][src][iEvent].Fill(jetIEta_toUse, l1jet_PU_pt/Raw_TT_iET_tmp, puWeight)
                        #hist2['jet_byHand_PUByRawPt_vs_iEta%s' % (jetShape1)][algo1]['HBEF'        ]['PtAllBins'][src][iEvent].Fill(jetIEta_toUse, l1jet_PU_pt/Raw_TT_iET_tmp, puWeight)
                        
                        ## Calibration plots with fixBinWidthPFJetPt
                        #hist2['jet_byHand_L1JetPt_vs_PFJetPt%s' % (jetShape1)][algo1][sjetIEta_toUse]['PtAllBins'][src][iEvent].Fill(l1jet_pt, vOff.Pt(), puWeight)
                        #hist2['jet_byHand_L1JetPt_vs_PFJetPt%s' % (jetShape1)][algo1]['HBEF'        ]['PtAllBins'][src][iEvent].Fill(l1jet_pt, vOff.Pt(), puWeight)
                        '''
                        if runMode in ['CalibJetByHand'] and jetShape == 'Default' and algo1 == 'RawPUS':
                            # validate Layer2Calibration by hand
                            #print "jetShape1: {}, algo1 {}, str(jetIEta) {}, jPt {}, src {}, iEvent {}".format('jet_byHand_den%s' % (jetShape1), algo1, str(jetIEta), jPt, src, iEvent)
                            #hist2['jet_byHand_L1JetPt_vs_DefaultL1JetPt%s' % (jetShape1)][algo1][sjetIEta_toUse]['PtAllBins'][src][iEvent].Fill(l1jet_pt_woLayer2Calib, (l1jet_pt - jetEt) / jetEt, puWeight)
                            #hist2['jet_byHand_L1JetPt_vs_DefaultL1JetPt%s' % (jetShape1)][algo1]['HBEF'        ]['PtAllBins'][src][iEvent].Fill(l1jet_pt_woLayer2Calib, (l1jet_pt - jetEt) / jetEt, puWeight)
                            pass
                        '''
                        
                        '''
                        # w.r.t. nVts
                        if iL1JetPtCat != 'None':
                            if 'jet_byHand_res_vs_iEta_vs_nVtx' in dists2:
                                hist2['jet_byHand_res_vs_iEta_vs_nVtx%s' % (jetShape1)][algo1]['HBEF'        ][iL1JetPtCat        ][src][iEvent].Fill(jetIEta_toUse, nVtx, res, puWeight)
                        if 'jet_byHand_res_vs_iEta_vs_nVtx' in dists2:
                            hist2['jet_byHand_res_vs_iEta_vs_nVtx%s' % (jetShape1)][algo1]['HBEF'        ]['PtAllBins'][src][iEvent].Fill(jetIEta_toUse, nVtx, res, puWeight)
                        '''
                        
                        if 'jet_byHand_res_vs_iEta_vs_nVtx' in dists2:        #revisit
                            hist2['jet_byHand_res_vs_iEta_vs_nVtx%s' % (jetShape1)][algo1]['HBEF'        ][iPFJetPtCat][src][0].Fill(jetIEta_toUse, nVtx, res, puWeight)
                            hist2['jet_byHand_res_vs_iEta_vs_nVtx%s' % (jetShape1)][algo1]['HBEF'        ]['PtAllBins'][src][0].Fill(jetIEta_toUse, nVtx, res, puWeight)
                        if 'jet_byHand_res_woLayer2Calib_vs_iEta_vs_nVtx' in dists2:
                            hist2['jet_byHand_res_woLayer2Calib_vs_iEta_vs_nVtx%s' % (jetShape1)][algo1]['HBEF'        ][iPFJetPtCat][src][0].Fill(jetIEta_toUse, nVtx, res_woLayer2Calib, puWeight)
                            hist2['jet_byHand_res_woLayer2Calib_vs_iEta_vs_nVtx%s' % (jetShape1)][algo1]['HBEF'        ]['PtAllBins'][src][0].Fill(jetIEta_toUse, nVtx, res_woLayer2Calib, puWeight)
                            
                            

                #if runMode in ['makeInputForML'] and data_dict['L1JetDefault_RawEtPUS'] < 0.:
                #    print "-ve RawEtPUS: data_dict: {}".format(data_dict)
                if runMode in ['makeInputForML'] and isFirstEntry_WriteInputForML:
                    fOut_MLInputs_writer = csv.DictWriter(fOut_MLInputs, fieldnames=data_dict.keys())
                    fOut_MLInputs_writer.writeheader()
                    isFirstEntry_WriteInputForML = False
                    if PrintLevel >= 14:
                        print("WriteInputForML: data_dict.keys(): {}".format(data_dict.keys()) )
                    
                if runMode in ['makeInputForML']: fOut_MLInputs_writer.writerow( data_dict )
                if PrintLevel >= 14:
                    print("WriteInputForML: data_dict: {}".format(data_dict))


        # end loop: for iOff in range(nOffJets):
        if sFInEventsToRun and hCaloTowers_iEta_vs_iPhi and hCaloTTs_iEta_vs_iPhi:
            hCaloTowers_iEta_vs_iPhi_list.append( hCaloTowers_iEta_vs_iPhi )
            hCaloTTs_iEta_vs_iPhi_list.append( hCaloTTs_iEta_vs_iPhi ) 
        # ## End loop: for iOff in range(nOffJets):


        # ### Fill trigger rates plots per event level --------------------------------------------------------
        for src in ['unp','emu']:

            #for jetShape in ['Default'] + JetShapes:
            for jetShape in JetShapes + JetShapesType2:
                # JetShape = "" plots are with the first version of code for 9x9 jets
                jetShape1 = jetShape
                if jetShape == 'Default':  jetShape1 = ""
                else:                      jetShape1 = "_%s" % (jetShape)
                
                for algo1 in PUSAlgosAll + PUSAlgosAllType2:
                    # read proper jetShape and PUSAlgo conbination
                    if (jetShape in JetShapes      and algo1 not in PUSAlgosAll) or \
                        (jetShape in JetShapesType2 and algo1 not in PUSAlgosAllType2 ):
                        continue
                    
                    if (algo1 == 'L1JDefault') and (jetShape != 'Default'): continue
                    
                    for ieta_cat in IETA_CAT.keys():
                        if ieta_cat == 'HBEF': continue
                        
                        l1JetsInEvent_sortedByL1JPt = sorted(l1JetCollection[src][jetShape][algo1][ieta_cat], key=lambda x: x[0], reverse=True)
                        if len(l1JetsInEvent_sortedByL1JPt) == 0: continue
                        
                        if PrintLevel >= 5:
                            print("    l1JetCollection[{}][{}][{}][{}]: {} \t\t sorted {}".format(src, jetShape, algo1,  ieta_cat,  l1JetCollection[src][jetShape][algo1][ieta_cat], l1JetsInEvent_sortedByL1JPt))
                            
                        # leading jet
                        l1JetPt_toConsider      = l1JetsInEvent_sortedByL1JPt[0][0]
                        l1JetIEtaCat_toConsider = l1JetsInEvent_sortedByL1JPt[0][1]
                        for jetBin in range(jetRate_bins[0]):
                            pT_thrsh = jetRate_bins[1] + (jetBin * jetRate_binWidth)
                            if l1JetPt_toConsider <= pT_thrsh: continue
                            # if 'jet_byHand_rates_singleJet' in dists4:
                                # hist4['jet_byHand_rates_singleJet%s' % (jetShape1)][algo1][l1JetIEtaCat_toConsider][src][iEvent].Fill( pT_thrsh, nVtx, puWeight )
                                # hist4['jet_byHand_rates_singleJet%s' % (jetShape1)][algo1]['HBEF'                 ][src][iEvent].Fill( pT_thrsh, nVtx, puWeight )
                        
                        if len(l1JetsInEvent_sortedByL1JPt) <= 1: continue
                        # 2nd leading jet
                        l1JetPt_toConsider      = l1JetsInEvent_sortedByL1JPt[1][0]
                        l1JetIEtaCat_toConsider = l1JetsInEvent_sortedByL1JPt[1][1]
                        for jetBin in range(jetRate_bins[0]):
                            pT_thrsh = jetRate_bins[1] + (jetBin * jetRate_binWidth)
                            if l1JetPt_toConsider <= pT_thrsh: continue
                            # if 'jet_byHand_rates_doubleJet' in dists4:
                                # hist4['jet_byHand_rates_doubleJet%s' % (jetShape1)][algo1][l1JetIEtaCat_toConsider][src][iEvent].Fill( pT_thrsh, nVtx, puWeight )
                                # hist4['jet_byHand_rates_doubleJet%s' % (jetShape1)][algo1]['HBEF'                 ][src][iEvent].Fill( pT_thrsh, nVtx, puWeight )
                        
                        if len(l1JetsInEvent_sortedByL1JPt) <= 2: continue
                        # 3rd leading jet
                        l1JetPt_toConsider      = l1JetsInEvent_sortedByL1JPt[2][0]
                        l1JetIEtaCat_toConsider = l1JetsInEvent_sortedByL1JPt[2][1]
                        for jetBin in range(jetRate_bins[0]):
                            pT_thrsh = jetRate_bins[1] + (jetBin * jetRate_binWidth)
                            if l1JetPt_toConsider <= pT_thrsh: continue
                            if 'jet_byHand_rates_trippleJet' in dists4:       #revisit
                                hist4['jet_byHand_rates_trippleJet%s' % (jetShape1)][algo1][l1JetIEtaCat_toConsider][src][0].Fill( pT_thrsh, nVtx, puWeight )
                                hist4['jet_byHand_rates_trippleJet%s' % (jetShape1)][algo1]['HBEF'                 ][src][0].Fill( pT_thrsh, nVtx, puWeight )
                        
                        if len(l1JetsInEvent_sortedByL1JPt) <= 3: continue
                        # 4th leading jet
                        l1JetPt_toConsider      = l1JetsInEvent_sortedByL1JPt[3][0]
                        l1JetIEtaCat_toConsider = l1JetsInEvent_sortedByL1JPt[3][1]
                        for jetBin in range(jetRate_bins[0]):
                            pT_thrsh = jetRate_bins[1] + (jetBin * jetRate_binWidth)
                            if l1JetPt_toConsider <= pT_thrsh: continue
                            if 'jet_byHand_rates_quadJet' in dists4:      #revisit
                                hist4['jet_byHand_rates_quadJet%s' % (jetShape1)][algo1][l1JetIEtaCat_toConsider][src][0].Fill( pT_thrsh, nVtx, puWeight )
                                hist4['jet_byHand_rates_quadJet%s' % (jetShape1)][algo1]['HBEF'                 ][src][0].Fill( pT_thrsh, nVtx, puWeight )
            ###  trigger rates plots per event level --------------------------------------------------------                


            ## End loop: for jEvt in range(chains['Unp'][iEvent].GetEntries()):

        print("\n\n iEvent: {} ".format( nTotalEvents_byChains[0]))
        hnTotalEvents.SetBinContent(1, nTotalEvents_byChains[0])
        # ## End loop: for iEvent in range(len(chains['Unp'])):
    # print(f"Total Events: {NTot}")
    print('\nFinished loop over Events')
    print(f"Processing completed in {time.time() - start_time:.2f} seconds")
    if runMode in ['makeInputForML']: 
        fOut_MLInputs.close()

    out_file = R.TFile(out_file_str,'recreate')
    out_file.cd()

    colors = [R.kGray, R.kViolet, R.kBlue, R.kCyan+3, R.kTeal, R.kSpring, R.kRed, R.kOrange, R.kMagenta+3, R.kMagenta]
    while len(colors) < len(in_file_names):
        colors += colors

    if usePUReweighting:
        dir1 = out_file.mkdir("PUWeights")
        dir1.cd()
        hPUWt.Write();
        hnVtxData.Write();
        hnVtxMC.Write();
        hnVtxMCWtd.Write();

    out_file.cd()
    hnVtx.Write();
    hnVtx_ReWtd.Write();
    hnTotalEvents.Write();


    '''
    hist_L1Jet_unp_TowerIEta_vs_IEta.Write()
    hist_L1Jet_emu_TowerIEta_vs_IEta.Write()
    hist_L1Jet_unp_TowerIPhi_vs_IPhi.Write()
    hist_L1Jet_emu_TowerIPhi_vs_IPhi.Write()

    for iEta in ETA_Bins:
        hist_PFJetPt_iEtawise[iEta].Write()

    for src in ['unp','emu']:    
        for iEta in ETA_Bins:
            hist_nPV_vs_L1JetDefaultRAW_SF[src][iEta].Write()
            hist_nPV_vs_L1JetDefaultPUS_SF[src][iEta].Write()
    '''
    l1MatchOffline = True
    if   l1MatchOffline:
        hdR_OffJet_OffMu_min.Write();
        hdR_OffJet_OffEle_min.Write();
        hdR_OffJet_OffMu_min_vs_vOffMuPtByvOffJetPt.Write();
        hdR_OffJet_OffEle_min_vs_vOffElePtByvOffJetPt.Write();
        hdR_OffJet_OffMu_min_forPtFracLt0p5.Write();
        hdR_OffJet_OffEle_min_forPtFracLt0p5.Write();

        for src in ['unp','emu']:
            hdR_OffJet_L1Jet[src].Write();
        

    hRefJet_pt_0.Write();
    hRefJet_eta_0.Write();
    hRefJet_phi_0.Write(); #
    hRefJet_pt_0_1.Write();
    hRefJet_eta_0_1.Write();
    hRefJet_phi_0_1.Write(); #
    for idx_ in range( nJetFilters + 1 ):
        hRefJet_pt_0_1_test[idx_].Write()
        hRefJet_eta_0_1_test[idx_].Write()
        hRefJet_phi_0_1_test[idx_].Write()
    hRefJet_pt_0_2.Write();
    hRefJet_eta_0_2.Write();
    hRefJet_phi_0_2.Write(); #   
    hRefJet_pt_0_3.Write();
    hRefJet_eta_0_3.Write();
    hRefJet_phi_0_3.Write(); #
    hRefJet_pt_0_4.Write();
    hRefJet_eta_0_4.Write();
    hRefJet_phi_0_4.Write(); #    
    hRefJet_pt_1.Write();
    hRefJet_eta_1.Write();
    hRefJet_phi_1.Write(); #
    hRefJet_pt_2.Write();
    hRefJet_eta_2.Write();
    hRefJet_phi_2.Write(); #
    hRefJet_pt_2p01.Write();
    hRefJet_eta_2p01.Write();
    hRefJet_phi_2p01.Write(); #
    hRefJet_pt_2p02.Write();
    hRefJet_eta_2p02.Write();
    hRefJet_phi_2p02.Write(); #
    hRefJet_pt_2p03.Write();
    hRefJet_eta_2p03.Write();
    hRefJet_phi_2p03.Write(); #

    hL1JetUnp_Pt_0.Write()
    hL1JetUnp_Eta_0.Write()
    hL1JetUnp_Phi_0.Write() #
    hL1JetEmu_Pt_0.Write()
    hL1JetEmu_Eta_0.Write()
    hL1JetEmu_Phi_0.Write() #
    hOfflineJet_Pt_0.Write()
    hOfflineJet_Eta_0.Write()
    hOfflineJet_Phi_0.Write() #

    hL1JetUnp_Pt_1.Write()
    hL1JetUnp_Eta_1.Write()
    hL1JetUnp_Phi_1.Write() #
    hL1JetEmu_Pt_1.Write()
    hL1JetEmu_Eta_1.Write()
    hL1JetEmu_Phi_1.Write() #
    hL1JetUnp_Pt_2.Write()
    hL1JetUnp_Eta_2.Write()
    hL1JetUnp_Phi_2.Write() #
    hL1JetEmu_Pt_2.Write()
    hL1JetEmu_Eta_2.Write()
    hL1JetEmu_Phi_2.Write() #        

    ## Loop over all histograms
    if runMode not in ['CalCalibSF', 'CalibJetByHand', 'makeInputForML'] and len(hist) > 0:
        out_file.cd()
        for algo in ['PUS','noPUS','Raw','RawPUS']:
            for iEta in ETA_CAT.keys():
                for iPt in PT_CAT.keys():
                    for src in ['unp','emu']:
                        # for iTP in range(len(in_file_names)):

                        hist['jet_eff'][algo][iEta][iPt][src] = R.TEfficiency( hist['jet_num'][algo][iEta][iPt][src],
                                                                                    hist['jet_den'][algo][iEta][iPt][src] )
                        hist['jet_eff'][algo][iEta][iPt][src].SetName( hist['jet_num'][algo][iEta][iPt][src].GetName().replace('num','eff') )

                        for dist in dists:

                            hist[dist][algo][iEta][iPt][src].SetLineWidth(2)
                            if src == 'unp': hist[dist][algo][iEta][iPt][src].SetLineColor(R.kBlack)
                            if src == 'emu': hist[dist][algo][iEta][iPt][src].SetLineColor(colors)
                            hist[dist][algo][iEta][iPt][src].Write()

                            ## End loop: for dist in dists+['jet_eff']
                        ## End loop: for iTP in range(len(in_file_names))
                    ## End loop: for src in ['unp','emu']
                ## End loop: for iPt in PT_CAT.keys()
            ## End loop: for iEta in ETA_CAT.keys()
        ## End loop: for algo in ['PUS','noPUS','Raw','RawPUS']
        
    if runMode not in ['CalCalibSF', 'CalibJetByHand', 'makeInputForML'] and len(hist1) > 0:
        out_file.cd()
        for dist in dists1:
            for iEta in ETA_Bins:
                for src in ['unp','emu']:
                    # for iTP in range(len(in_file_names)):
                    hist1[dist][iEta][src].SetLineWidth(2)
                    if src == 'unp': hist1[dist][iEta][src].SetLineColor(R.kBlack)
                    if src == 'emu': hist1[dist][iEta][src].SetLineColor(colors)
                    hist1[dist][iEta][src].Write()


    ## Loop over all histograms
    # hist2
    #if runMode not in ['CalibJetByHand', 'makeInputForML']:
    if runMode not in ['makeInputForML']:
        out_file.cd()        
        for algo in PUSAlgosAll: # ['Raw', 'RawPUS', 'RawPUS_phiRingMin4', 'RawPUS_phiRingSide4', 'RawPUS_phiRingAdjacent']: # ['PUS','noPUS','Raw','RawPUS']:        
            for iEta in ETA_Bins:
                for iPt in list(PT_CAT.keys()) + ['PtAllBins']:
                    for src in ['unp','emu']: #['unp','emu']:
                        # for iTP in range(len(in_file_names)):
                            #for jetShape in ['Default'] + JetShapes:
                        for jetShape in JetShapes:
                            if (algo == 'L1JDefault') and (jetShape != 'Default'): continue
                            if (not l1nanoChunkyDonut) and (jetShape == 'Default') and (algo == 'RawPUS'):            continue
                            if (not l1nanoPhiRing)     and (jetShape == 'Default') and (algo == 'RawPUS_phiDefault'): continue

                            # JetShape = "" plots are with the first version of code for 9x9 jets
                            jetShape1 = jetShape
                            if jetShape == 'Default':  jetShape1 = ""
                            else:                      jetShape1 = "_%s" % (jetShape)
                                # JetShape = "" plots are with the first version of code for 9x9 jets                        

                            if 'jet_byHand_eff' in dists2:
                                hist2['jet_byHand_eff%s' % (jetShape1)][algo][iEta][iPt][src] = R.TEfficiency( hist2['jet_byHand_num%s' % (jetShape1)][algo][iEta][iPt][src],
                                                                                                    hist2['jet_byHand_den%s' % (jetShape1)][algo][iEta][iPt][src] )
                                hist2['jet_byHand_eff%s' % (jetShape1)][algo][iEta][iPt][src].SetName( hist2['jet_byHand_num%s' % (jetShape1)][algo][iEta][iPt][src].GetName().replace('num','eff') )

                            for dist_1 in dists2:
                                dist = '%s%s' % (dist_1, jetShape1)
                                if '_vs_iEta' in dist and iEta != 'HBEF': continue
                                if dist_1 in ['jet_byHand_L1JetPt_vs_DefaultL1JetPt'] and algo != 'RawPUS': continue
                                if dist_1 in ['jet_byHand_L1JetPt_vs_PFJetPt', 'jet_byHand_L1JetPt_vs_DefaultL1JetPt'] and iPt != 'PtAllBins': continue
                                if dist_1 in ['jet_byHand_L1JetPt_vs_DefaultL1JetPt'] and jetShape != 'Default': continue

                                hist2[dist][algo][iEta][iPt][src].SetLineWidth(2)
                                if src == 'unp': hist2[dist][algo][iEta][iPt][src].SetLineColor(R.kBlack)
                                if src == 'emu': hist2[dist][algo][iEta][iPt][src].SetLineColor(colors)
                                hist2[dist][algo][iEta][iPt][src].Write()



    if runMode in ['CalibJetByHand']:
        # hist3
        out_file.cd()        
        for algo in PUSAlgosAll + PUSAlgosAllType2: # ['Raw', 'RawPUS', 'RawPUS_phiRingMin4', 'RawPUS_phiRingSide4', 'RawPUS_phiRingAdjacent']: # ['PUS','noPUS','Raw','RawPUS']:
            for iEta in ETA_CAT.keys():
                for iPt in PT_CAT.keys():
                    for src in ['unp','emu']: #['unp','emu']:
                        # for iTP in range(len(in_file_names)):
                            #for jetShape in ['Default'] + JetShapes:
                        for jetShape in JetShapes + JetShapesType2:
                            # read proper jetShape and PUSAlgo conbination
                            if (jetShape in JetShapes      and algo not in PUSAlgosAll) or \
                                (jetShape in JetShapesType2 and algo not in PUSAlgosAllType2 ):
                                continue
            
                            if (algo == 'L1JDefault') and (jetShape != 'Default'): continue

                            # JetShape = "" plots are with the first version of code for 9x9 jets
                            jetShape1 = jetShape
                            if jetShape == 'Default':  jetShape1 = ""
                            else:                      jetShape1 = "_%s" % (jetShape)
                                # JetShape = "" plots are with the first version of code for 9x9 jets   

                            for dist_1 in dists3:
                                dist = '%s%s' % (dist_1, jetShape1)

                                hist3[dist][algo][iEta][iPt][src].SetLineWidth(2)
                                if src == 'unp': hist3[dist][algo][iEta][iPt][src].SetLineColor(R.kBlack)
                                if src == 'emu': hist3[dist][algo][iEta][iPt][src].SetLineColor(colors)
                                hist3[dist][algo][iEta][iPt][src].Write()

        # hist4
        out_file.cd()        
        for algo in PUSAlgosAll + PUSAlgosAllType2: # ['Raw', 'RawPUS', 'RawPUS_phiRingMin4', 'RawPUS_phiRingSide4', 'RawPUS_phiRingAdjacent']: # ['PUS','noPUS','Raw','RawPUS']:
            for iEta in ETA_CAT.keys():
                for src in ['unp','emu']: #['unp','emu']:
                    # for iTP in range(len(in_file_names)):
                    #for jetShape in ['Default'] + JetShapes:
                    for jetShape in JetShapes + JetShapesType2:
                        # read proper jetShape and PUSAlgo conbination
                        if (jetShape in JetShapes      and algo not in PUSAlgosAll) or \
                            (jetShape in JetShapesType2 and algo not in PUSAlgosAllType2 ):
                            continue
                        
                        if (algo == 'L1JDefault') and (jetShape != 'Default'): continue
                        
                        # JetShape = "" plots are with the first version of code for 9x9 jets
                        jetShape1 = jetShape
                        if jetShape == 'Default':  jetShape1 = ""
                        else:                      jetShape1 = "_%s" % (jetShape)
                            # JetShape = "" plots are with the first version of code for 9x9 jets   

                        for dist_1 in dists4:
                            dist = '%s%s' % (dist_1, jetShape1)

                            hist4[dist][algo][iEta][src].SetLineWidth(2)
                            if src == 'unp': hist4[dist][algo][iEta][src].SetLineColor(R.kBlack)
                            if src == 'emu': hist4[dist][algo][iEta][src].SetLineColor(colors)
                            hist4[dist][algo][iEta][src].Write()


    # JetShapesType2, PUSAlgosAllType2
    # hist5
    if runMode not in ['makeInputForML']:
        out_file.cd()
        hStat.Fill(50)
        
        for iEta in ETA_CAT.keys():
            for iPt in list(PT_CAT.keys()) + ['PtAllBins']:
                for src in ['emu']: # ['unp', 'emu']:
                    # for iTP in range(len(in_file_names)):
                    for jetShape in JetShapesType2:
                        jetShape1 = "_%s" % (jetShape)      
                        
                        for algo in PUSAlgosAllType2: # ['Et', 'RawEt']
                        
                            for dist_1 in dists5:
                                dist = '%s%s' % (dist_1, jetShape1)
                                

                                if dist_1 in ['jet_byHand_res_vs_iEta_vs_nVtx'] and \
                                    iEta != 'HBEF':
                                    continue

                                #print "dist5 {} algo {}, iEta {}', iPt {}, src {}, iTP {}".format(dist, algo,iEta, iPt,src,iTP)
                                hist5[dist][algo][iEta][iPt][src].Write()


    # hist5
    if runMode not in ['makeInputForML']:
        out_file.cd()

        for src in ['unp','emu']:
            # for iTP in range(len(in_file_names)):
            for dist in dists6:
                
                hist6[dist][src].Write()



    # hist7
    if runMode in ['']:
        out_file.cd()

        for src in ['emu']: # ['unp','emu']:
            # for iTP in range(len(in_file_names)):
            for dist in dists7:
                
                hists7[dist][src].Write()




    if runMode in ['trbshtPhiRingPUS']:
        out_file.cd()
        outDir1 = out_file.mkdir( "caloTTs" )
        outDir1.cd()
        
        if len(hCaloTowers_iEta_vs_iPhi_list) != len(hCaloTTs_iEta_vs_iPhi_list):
            print("len(hCaloTowers_iEta_vs_iPhi_list) != len(hCaloTTs_iEta_vs_iPhi_list):: something went wrong \t\t\t **** ERROR ****")
        else:
            for iH in range(len(hCaloTowers_iEta_vs_iPhi_list)):
                hCaloTowers_iEta_vs_iPhi_list[iH].Write()
                hCaloTTs_iEta_vs_iPhi_list[iH].Write()
                

    if runMode in ['trbshtPhiRingPUS']:
        out_file.cd()
        outDir1 = out_file.mkdir( "caloTTs_1" )
        outDir1.cd()

        hTTEtMax_forL1JetPUEt0.Write()
        hL1JetRawEt_vs_L1JetEt_forL1JetPUEt0.Write()
        

    # hist8
    if runMode in ['', 'makeInputForML']:
        out_file.cd()
        outDir1 = out_file.mkdir( "L1JetPt" )
        outDir1.cd()

        for src in ['unp','emu']:
            for iEta in ETA_Bins:
                for dist in dists8:
                    hist8[dist][src][iEta].Write()

        

    out_file.cd()
    hStat.Write()

    hnVts_vs_nTT_unp.Write()
    hnVts_vs_nTT_emu.Write()

    hJEC_iEta_vs_Pt.Write()

    out_file.Close()

    #print '\nWrote out file:  plots/'+out_file_str+'.root'
    print('\nWrote out file:  '+out_file_str)


    print("selectedEvents_list: {}".format(selectedEvents_list))
    if sFOutSelectedEvents and len(selectedEvents_list) > 0: 
        with open(sFOutSelectedEvents, 'w') as fOutSelectedEvents:
            for evt in selectedEvents_list:
                fOutSelectedEvents.write( '%s\n' % evt )

    print(f"Processing completed in {time.time() - start_time:.2f} seconds")


if __name__ == '__main__':
    run()
