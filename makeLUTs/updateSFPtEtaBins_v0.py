#!/bin/env python

'''
Update Pt-Eta bins for SFs to be compatible with firware requirements
'''


import os
import sys
import math 
from array import array
import ROOT as R
R.gROOT.SetBatch(False)  ## Don't print histograms to screen while processing

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

PrintLevel       = 2
JetShapes        = ['Default' ]
PUSAlgosAll      = ['RawPUS']  #['RawPUS', 'RawPUS_phiDefault']

nLinesToRead = -1 # Max line to read from SF.csv file

# Read L1Jet CalibLayer2 SFs from csv files provided by Syed
useAdditionalSFForLUT = True # 8/7 factor is needed for PhiRing as different PU estimation considered in CMSSW and Andrew's computation
icalibSF = 0 # 0, 1
calibSFLable = ['SF'][icalibSF]  
sipFileCalibSF = {
    'Default': {
       'RawPUS': { # Chunky donut
           'fileName': '../data/L1T_Jet_SFs_Run3_QCD_Pt15to7000_PFA1p_CMSSW12_6_0_pre1_nVtxAll_20220925_v3_HyperOpt_test_L1JetEt_PUS_ChunkyDonut.csv', 
           'SFLabel': ['SF'][icalibSF],
           'L1JetPtVarName':'L1JetEt_PUS_ChunkyDonut',
           'additionalFactorForLUT': 1.0,
       },
       
       'RawPUS_phiDefault': {
           'fileName': '../data/L1T_Jet_SFs_Run3_QCD_Pt15to7000_PFA1p_CMSSW12_6_0_pre1_nVtxAll_20220925_v3_HyperOpt_test_L1JetEt_PUS_ChunkyDonut.csv',
           'SFLabel': ['ScaleFactor_Log(GenJet)', 'ScaleFactor_Log(PhiRing)/Log(GenJet)'][icalibSF],
           'L1JetPtVarName':'PhiRingEnergy',
           'additionalFactorForLUT': 8.0/7.0,
       },
    },


}

'''
#calibSF_additionalCorr = 8./7
#updateCalibSF_DefaultRaw = True # update calibSFs_additionalCorr['Default']['Raw'] depending on 
calibSFs_additionalCorr = {
    'Default': { # trigger tower Et was scalled by 0.5 while calculating l1Jet_pT by-hand and to calculate layer2 calibSFs
        'Raw': {'l1NtupleChunkyDonut': 1. * 0.5,  'l1NtuplePhiRing': 8./7 * 0.5 },
        'RawPUS': 1. * 0.5,
        'RawPUS_phiDefault': 8./7 * 0.5,
    },
    
}
'''

sLUTVersion = '2022v4'
nBinsMaxForEtaCompressionLUT = 64 # no. of lines in eta compression LUT
makeLUTForIEta29 = [False, 1.0]
makeLUTForIEta41 = [False] # SFs for iEta=41 are missing in SFv6. Copy SFs from IEta=40
separatePtBinForLowPt = True # True: Use separate pT bin quantile for 0 < pT <= 15 GeV where SFs are truckated
IEtaBinOffsetForEtaCompressedLUT = 0 # 0: IEtaBin_forLUT = IEtaBin = [1, 41];  -1: IEtaBin_forLUT = IEtaBin - 1 = [0, 40]. # 0 give correct calibration.
sFOut_LUT_pt_compress    = 'lut_pt_compress_%s.txt' % (sLUTVersion)
sFOut_LUT_eta_compress   = 'lut_eta_compress_%s.txt' % (sLUTVersion)
sFOut_LUT_calib_compress = 'lut_calib_%s_ECALZS_decimal.txt' % (sLUTVersion)

makeLUTsInUncompressedBins = False # make LUTs without compressing Pt and Eta bins
sFOut_LUT_pt_uncompress  = 'lut_pt_uncompress_%s.txt' % (sLUTVersion)
sFOut_LUT_eta_uncompress = 'lut_eta_uncompress_%s.txt' % (sLUTVersion)
sFOut_LUT_calib_uncompress = 'lut_calib_%s_ECALZS_decimal_uncompress.txt' % (sLUTVersion)


NCompPtBins = 16 # 16 # No. of compressed pT bins
calibSF_L1JetPtRange = [15., 255., 1.] #[15., 250., 1.] # [<lowest pT>,  <hightest pT>,  <pT bin width>] # pT range for SFs to read from Syed's SF.csv file
LUT_PtRange = [0., 255., 1.] # pT range for SFs for LUT
SF_forZeroPt = 1.

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


IEtaBinsMerge_forPlots = [
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




useAbsEtaBins = True
ETA_Bins = []
for iEta in range(-41,42):
    if iEta in [-29, 0, 29]:        continue;
    if useAbsEtaBins and iEta < 0:  continue;
    ETA_Bins.append(str(iEta))




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
        print("idx_pTBinQuantiles final ({}): {}".format(
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

    
    
    

if __name__ == '__main__':


    calibSF_L1JetPtRangeMin = calibSF_L1JetPtRange[0]
    calibSF_L1JetPtRangeMax = calibSF_L1JetPtRange[1]
    calibSF_L1JetPtBinWidth = calibSF_L1JetPtRange[2]

    
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
            

            data_calibSFs = OrderedDict()
            for Eta_bin in ETA_Bins:
                data_calibSFs[int(Eta_bin)] = OrderedDict()

            with open(sipFileCalibSF_toRead, mode='r') as fipFileCalibSF_toRead:
                calibSF_csv_reader = csv.DictReader(fipFileCalibSF_toRead)

                nLines = 0
                
                #print("\n\ncalibSF_csv_reader: {} {}\n\n\n".format(type(calibSF_csv_reader), calibSF_csv_reader))
                for calibSF_csv_row in calibSF_csv_reader:
                    iEta_tmp = int(calibSF_csv_row['L1JetTowerIEtaAbs']) # str
                    #l1JetPt_tmp = float(calibSF_csv_row['PhiRingEnergy'])
                    l1JetPt_tmp = float(calibSF_csv_row[ sipFileCalibSF[jetShape][PUSAlgo]['L1JetPtVarName'] ])
                    SF_tmp = float(calibSF_csv_row[ sipFileCalibSF[jetShape][PUSAlgo]['SFLabel'] ])
                    if useAdditionalSFForLUT:
                        SF_tmp = SF_tmp * sipFileCalibSF[jetShape][PUSAlgo]['additionalFactorForLUT']
                        
                    if PrintLevel >= 11:
                        print("calibSF_csv_row: {}, iEta {} {}, l1JetPt {} {}, SF {} {}={}".format(
                            calibSF_csv_row, type(iEta_tmp),iEta_tmp, type(l1JetPt_tmp),l1JetPt_tmp, type(SF_tmp),SF_tmp, sipFileCalibSF[jetShape][PUSAlgo]['SFLabel']))

                        
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
                        #data_calibSFs[Eta_bin][Pt] = 1.0
                        # Add Pt=0 at the beining in OrderedDict
                        data_calibSFs[Eta_bin].update( {Pt : 1.0} )
                        data_calibSFs[Eta_bin].move_to_end(Pt, last=False)
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
            print("\nmakeLUTsInUncompressedBins:: \nEtaBins_sorted ({}): {} \nPtBins_forLUT ({}): {}".format(
                len(EtaBins_sorted), EtaBins_sorted,
                nPtBins_forLUT, PtBins_forLUT
            ))

            
            # write LUTs with uncompressed Pt and eta -------------------------------------------------------------
            if makeLUTsInUncompressedBins:
                sDir_LUT = "LUTs/%s_%s_%s" % (jetShape, PUSAlgo, calibSFLable)
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
                fOut_LUT_eta_uncompress.write('# Converts abs(MP ieta) (6 bits) into 4 bit index\n')
                fOut_LUT_eta_uncompress.write('# This is NOT calo ieta\n')
                fOut_LUT_eta_uncompress.write('# anything after # is ignored with the exception of the header\n')
                fOut_LUT_eta_uncompress.write('# the header is first valid line starting with #<header> versionStr nrBitsAddress nrBitsData </header>\n')
                fOut_LUT_eta_uncompress.write('#<header> v1 6 4 </header>\n')

                fOut_LUT_pt_uncompress.write('# PT compression LUT\n')
                fOut_LUT_pt_uncompress.write('# maps 8 bits to 4 bits\n')
                fOut_LUT_pt_uncompress.write('# the 1st column is the integer value after selecting bits 1:8\n')
                fOut_LUT_pt_uncompress.write('# anything after # is ignored with the exception of the header\n')
                fOut_LUT_pt_uncompress.write('# the header is first valid line starting with #<header> versionStr nrBitsAddress nrBitsData </header>\n')
                fOut_LUT_pt_uncompress.write('#<header> v1 8 4 </header>\n')
                
                

                for IEta in EtaBins_sorted:
                    IEtaBin_forLUT = IEta - 1
                    fOut_LUT_eta_uncompress.write('%d %d\n' % (int(IEtaBin_forLUT), int(IEtaBin_forLUT)))
                    
                    for Pt in PtBins_forLUT:
                        PtBin_forLUT = Pt

                        
                        if IEta == EtaBins_sorted[0]:
                            fOut_LUT_pt_uncompress.write('%d %d\n' % (int(Pt), int(PtBin_forLUT)))
                            
                        sComments = ""
                        if Pt == PtBins_forLUT[0]:
                            sComments = "   # eta_bin %d, pt %d" % (IEtaBin_forLUT, PtBin_forLUT)

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
                    print("data_calibSFs[Eta_bin] ({}): {} ".format(len(list(data_calibSFs[Eta_bin].keys())), data_calibSFs[Eta_bin].keys()))
                if len(list(data_calibSFs[Eta_bin].keys())) == 0:
                    print("Eta_bin {}: len(list(data_calibSFs[Eta_bin].keys())) == 0".format(Eta_bin))
                    continue

                if PrintLevel >= 5:
                    print("Eta_Bins_Read.index(Eta_bin)): {} ".format(Eta_Bins_Read.index(Eta_bin)))
                    
                idx_Eta_Bin_Read = Eta_Bins_Read.index(Eta_bin)
                
                SFs_inEtaBin = np.array( list(data_calibSFs[Eta_bin].values()) )
                #SFs_array[Eta_bin] = SFs_inEtaBin
                SFs_array[idx_Eta_Bin_Read] = SFs_inEtaBin

                if PrintLevel >= 5:
                    print("\nEta_bin {}: \nSF ({}): {} \nSF_all ({}): {}".format(
                        Eta_bin,
                        SFs_inEtaBin.shape, SFs_inEtaBin,
                        SFs_array.shape,SFs_array
                    ))

            print("\n\nSF_all final ({}): \n{}\n\n".format(SFs_array.shape,SFs_array))
            # ----------------------------------------------------------------------------------------------------------
            
            ## average SF over all eta bins
            SFs_AvgOverEtaBins = np.average(SFs_array, axis=0)
            print("\nSFs_AvgOverEtaBins ({}): {}".format(len(SFs_AvgOverEtaBins), SFs_AvgOverEtaBins))


            
            ### Divide pT in quantiles with constant DeltaSF over the quatiles  ---------------------------
            idx_quantiles, SFs_AvgOverEtaBins_splitInQuantile = splitInQuantile_ConstDeltaSF(Pt_Bins_Read, SFs_AvgOverEtaBins, nQuantiles=NCompPtBins)
            

            Pt_quantiles = []
            for idx in idx_quantiles:
                Pt_quantiles.append( (Pt_Bins_Read[idx-1] + Pt_Bins_Read[idx]) / 2. )
               
            if PrintLevel >= 2:
                print("idx_quantiles ({}): {}, \nPt_quantiles ({}): {} \nSFs_AvgOverEtaBins_splitInQuantile: {}".format(
                    len(idx_quantiles), idx_quantiles,
                    len(Pt_quantiles), Pt_quantiles,
                    len(SFs_AvgOverEtaBins_splitInQuantile), SFs_AvgOverEtaBins_splitInQuantile)
                      )

             
            SFs_AvgInPtQuant_AvgOverEtaBins = SFs_AvgOverEtaBins.copy()
            Pt_AvgInPtQuant = []
            SFs_AvgInPtQuant_AvgOverEtaBins_Coarse = []
            for iQuant in range(NCompPtBins):
                idxStart                  = 0             if iQuant == 0               else idx_quantiles[(iQuant-1)]
                idxEnd                    = nPt_Bins_Read if iQuant == (NCompPtBins-1) else idx_quantiles[(iQuant)]
                Pt_Bins_Read_iQuant       = Pt_Bins_Read[idxStart : idxEnd]
                SFs_AvgOverEtaBins_iQuant = SFs_AvgOverEtaBins[idxStart : idxEnd]
                avgPt_iQuant              = statistics.mean(Pt_Bins_Read_iQuant)
                avgSF_iQuant              = statistics.mean(SFs_AvgOverEtaBins_iQuant)
                if separatePtBinForLowPt and iQuant == 0:
                    avgSF_iQuant          = SFs_AvgOverEtaBins_iQuant[1]
                Pt_AvgInPtQuant.append( avgPt_iQuant )
                SFs_AvgInPtQuant_AvgOverEtaBins_Coarse.append( avgSF_iQuant )
                if PrintLevel >= 2:
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


                    
            sDir1 = "plots_%s_%s_%s" % (jetShape, PUSAlgo, calibSFLable)
            sDir1 = sDir1.replace("(", "_")
            sDir1 = sDir1.replace(")", "_")
            sDir1 = sDir1.replace("__", "_")
            if not os.path.exists(sDir1):
                os.makedirs(sDir1)


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





            sDir_LUT = "LUTs/%s_%s_%s" % (jetShape, PUSAlgo, calibSFLable)
            sDir_LUT = sDir_LUT.replace("(","_")
            sDir_LUT = sDir_LUT.replace(")","_")
            sDir_LUT = sDir_LUT.replace("__","_")
            if not os.path.exists(sDir_LUT):
                os.makedirs(sDir_LUT)

            fOut_LUT_pt_compress    = open("%s/%s" % (sDir_LUT, sFOut_LUT_pt_compress),    'w')
            fOut_LUT_eta_compress   = open("%s/%s" % (sDir_LUT, sFOut_LUT_eta_compress),   'w')
            fOut_LUT_calib_compress = open("%s/%s" % (sDir_LUT, sFOut_LUT_calib_compress), 'w')

            # write headers
            fOut_LUT_eta_compress.write('# MP ieta compression LUT\n')
            fOut_LUT_eta_compress.write('# Converts abs(MP ieta) (6 bits) into 4 bit index\n')
            fOut_LUT_eta_compress.write('# This is NOT calo ieta\n')
            fOut_LUT_eta_compress.write('# anything after # is ignored with the exception of the header\n')
            fOut_LUT_eta_compress.write('# the header is first valid line starting with #<header> versionStr nrBitsAddress nrBitsData </header>\n')
            #fOut_LUT_eta_compress.write('#<header> v1 6 4 </header>\n')
            fOut_LUT_eta_compress.write('#<header> v1 6 6 </header>\n')

            fOut_LUT_pt_compress.write('# PT compression LUT\n')
            fOut_LUT_pt_compress.write('# maps 8 bits to 4 bits\n')
            fOut_LUT_pt_compress.write('# the 1st column is the integer value after selecting bits 1:8\n')
            fOut_LUT_pt_compress.write('# anything after # is ignored with the exception of the header\n')
            fOut_LUT_pt_compress.write('# the header is first valid line starting with #<header> versionStr nrBitsAddress nrBitsData </header>\n')
            fOut_LUT_pt_compress.write('#<header> v1 8 4 </header>\n')


            ## Update SFs in each IEta bin to have in pT_quantiles --------------------------------
            if PrintLevel >= 2:
                print("Update SFs in each IEta bin to have in pT_quantiles".format())
                #print("\n\n data_calibSFs @10 : {} \n\n".format(data_calibSFs))
                
            data_calibSFs_compressed       = OrderedDict()
            #data_calibSFs_compressed_check = data_calibSFs.copy()
            data_calibSFs_compressed_check = copy.deepcopy(data_calibSFs)

            IEtaBin_forLUT_col0 = None
            for idx_EtaBins_sorted, IEtaBin in enumerate(EtaBins_sorted):               
                data_calibSFs_compressed[IEtaBin] = OrderedDict()

                if IEtaBin == 29: continue

                # CaloTools::mpEta(jet->hwEta())
                CaloTool_mpEta = map_CaloIEta_to_CaloTool_mpEta[IEtaBin]
                
                IEtaBin_forLUT_col0 = CaloTool_mpEta  # IEtaBin if IEtaBinOffsetForEtaCompressedLUT == 0 else IEtaBin - 1;
                IEtaBin_forLUT_col1 = CaloTool_mpEta - 1 #if IEtaBin <= 28 else CaloTool_mpEta - 2
                #fOut_LUT_eta_compress.write('%d %d\n' % (int(IEtaBin_forLUT), int(IEtaBin_forLUT)))
                if idx_EtaBins_sorted == 0:
                    fOut_LUT_eta_compress.write('0 0  # dummy ieta_bin\n')
                fOut_LUT_eta_compress.write('%d %d  # ieta %d\n' % (int(IEtaBin_forLUT_col0), int(IEtaBin_forLUT_col1), IEtaBin))



                for iPtQuant in range(NCompPtBins):                  
                    idxStart                  = 0             if iPtQuant == 0               else idx_quantiles[(iPtQuant-1)]
                    idxEnd                    = nPt_Bins_Read if iPtQuant == (NCompPtBins-1) else idx_quantiles[(iPtQuant)]
                    Pt_Bins_Read_iPtQuant     = Pt_Bins_Read[idxStart : idxEnd]

                    #SFs_iPtQuant              = list(data_calibSFs[IEtaBin].values())[idxStart : idxEnd]
                    SFs_iPtQuant              = list(data_calibSFs[IEtaBin].values()).copy()[idxStart : idxEnd]

                    PtMin_iPtQuant            = Pt_Bins_Read_iPtQuant[0]
                    PtMax_iPtQuant            = Pt_Bins_Read_iPtQuant[-1]
                    avgPt_iPtQuant            = statistics.mean( Pt_Bins_Read_iPtQuant )
                    avgSF_iPtQuant            = statistics.mean( SFs_iPtQuant )
                    if separatePtBinForLowPt and iPtQuant == 0:
                        # ignore SF(pT=0)=1 and set SF(pT<=15 GeV) = SF(pT=15GeV)
                        avgSF_iPtQuant        = SFs_iPtQuant[1]
                        avgPt_iPtQuant        = Pt_Bins_Read_iPtQuant[-1] # calculate SF_inBits at bin higher edge (pT=15 GeV)

                    #data_calibSFs_compressed[IEtaBin][iPtQuant] = avgSF_iPtQuant
                    data_calibSFs_compressed[IEtaBin][avgPt_iPtQuant] = avgSF_iPtQuant
                    
                    
                
                    #data_calibSFs_compressed_check[IEtaBin][idxStart : idxEnd] = [avgSF_iPtQuant] * (idxEnd - idxStart)
                    for idx_PtUncompressesBin in range(idxStart, idxEnd):
                        data_calibSFs_compressed_check[IEtaBin][idx_PtUncompressesBin] = avgSF_iPtQuant

                    
                    PtBin_forLUT = iPtQuant                    
                    if idx_EtaBins_sorted == 0:
                        for Pt in Pt_Bins_Read_iPtQuant:
                            fOut_LUT_pt_compress.write('%d %d\n' % (int(Pt), int(PtBin_forLUT)))
                            
                    sComments = ""
                    if iPtQuant == 0:
                        sComments = "   # ieta %d, pt %d" % (IEtaBin, PtBin_forLUT)
                    fOut_LUT_calib_compress.write('%d %d %d %d %d %d %f%s\n' % (int(IEtaBin), int(PtMin_iPtQuant + 0.5), int(PtMax_iPtQuant + 0.5), int(avgPt_iPtQuant + 0.5), int(IEtaBin_forLUT_col1), int(PtBin_forLUT), avgSF_iPtQuant, sComments ))

                if PrintLevel >= 2:
                    print("\ndata_calibSFs_compressed[{}]: ({}) {} \ndata_calibSFs_compressed_check[{}]: ({}) {}".format(
                        IEtaBin, len(data_calibSFs_compressed[IEtaBin].keys()), data_calibSFs_compressed[IEtaBin],
                        IEtaBin, len(data_calibSFs_compressed_check[IEtaBin].keys()), data_calibSFs_compressed_check[IEtaBin]
                    ))

                if PrintLevel >= 20:
                    print("\n\n data_calibSFs[{}] @11  : {} \n\n".format(IEtaBin, data_calibSFs[IEtaBin]))

            

            # dummy bins
            # nBinsMaxForEtaCompressionLUT
            #for IEtaBin_forLUT in range(len(EtaBins_sorted), nBinsMaxForEtaCompressionLUT):
            #    fOut_LUT_eta_compress.write('%d %d\n' % (int(IEtaBin_forLUT), int(IEtaBin_forLUT)))
            IEtaBin_forLUT_col0 += 1 # the next bin = 42
            while IEtaBin_forLUT_col0 < nBinsMaxForEtaCompressionLUT:
                IEtaBin_forLUT_col1 = IEtaBin_forLUT_col0 - 1 
                fOut_LUT_eta_compress.write('%d %d\n' % (int(IEtaBin_forLUT_col0), int(IEtaBin_forLUT_col1)))
                IEtaBin_forLUT_col0 += 1

            
            fOut_LUT_pt_compress.close()
            fOut_LUT_eta_compress.close()
            fOut_LUT_calib_compress.close()

            print("Wrote %s/%s, %s/%s, %s/%s" % (
                sDir_LUT, sFOut_LUT_pt_compress,
                sDir_LUT, sFOut_LUT_eta_compress,
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
                if IEtaBin == 29: continue
                
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

            print("Done !!!")

            
