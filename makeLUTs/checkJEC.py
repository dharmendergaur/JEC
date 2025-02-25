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


sInFileDir = "/home/siddhesh/Work/CMS/L1_Trigger_Work/L1T_ServiceTask/hcalPUsub/myAna/hcalPUsub_v5_20220311/run_1/makeLUTs/LUTs/Default_RawPUS_SF_RegressedTo_log_GenJetPt_v6";
sInFIle_lut_JECInDecimal = "lut_calib_2022v1_ECALZS_decimal.txt"
sInFIle_lut_PtCompression = "lut_pt_compress_2022v1.txt"

sInFIle_MonitorJetCalib = "L1T_HCALL2Calib_stage1_l1NtupleChunkyDonut_PFA1p_nVtxAll_part0_of_1.root"
sHistoName_JEC_iEta_vs_Pt = "hJEC_iEta_vs_Pt"


if __name__ == '__main__':

    sInFIle_lut_JECInDecimal_fullName  = '%s/%s' % (sInFileDir,sInFIle_lut_JECInDecimal)
    sInFIle_lut_PtCompression_fullName = '%s/%s' % (sInFileDir,sInFIle_lut_PtCompression)

    lut_pTCompression = OrderedDict()
    with open(sInFIle_lut_PtCompression_fullName, mode='r') as fIn_lut_PtCompression:
        for line0 in fIn_lut_PtCompression.readlines():
            line = line0.replace('\n','')
            colmns = line.split(" ")
            if '#' in colmns[0]: continue

            if PrintLevel >= 15:
                print("line {}, colmns {}".format(
                    line, colmns
                ))
            Pt = float(colmns[0])
            PtCompressionBin = int(colmns[1])
            lut_pTCompression[Pt] = PtCompressionBin

    if PrintLevel >= 2:
        print("lut_pTCompression: {}".format(lut_pTCompression))

        
    lut_JEC           = OrderedDict()
    with open(sInFIle_lut_JECInDecimal_fullName, mode='r') as fIn_lut_JECInDecimal:
        for line0 in fIn_lut_JECInDecimal.readlines():
            line = line0.replace('\n','')
            colmns = line.split(" ")
            if '#' in colmns[0]: continue

            if PrintLevel >= 15:
                print("line {}, colmns {}".format(
                    line, colmns
                ))
            IEta             = int(colmns[0])
            PtCompressionBin = int(colmns[5])
            SF               = float(colmns[6])

            if IEta not in lut_JEC.keys():
                lut_JEC[IEta] = OrderedDict()
            lut_JEC[IEta][PtCompressionBin] = SF

    if PrintLevel >= 2:
        print("lut_JEC: {}".format(lut_JEC))

        

    EtaBins_sorted = list(lut_JEC.keys())
    EtaBins_sorted.sort()
    PtBins = list(lut_pTCompression.keys()) # [float for Pt in list(lut_pTCompression.keys()) ]


    fIn_MonitorJetCalib = R.TFile(sInFIle_MonitorJetCalib)
    hJEC_iEta_vs_Pt = fIn_MonitorJetCalib.Get(sHistoName_JEC_iEta_vs_Pt)
    print("hJEC_iEta_vs_Pt.GetNbinsX() {}, hJEC_iEta_vs_Pt.GetNbinsY() {}".format(hJEC_iEta_vs_Pt.GetNbinsX(), hJEC_iEta_vs_Pt.GetNbinsY()))

    sDir1 = "plots_checkJEC"
    if not os.path.exists(sDir1):
        os.makedirs(sDir1)

    
    for IEta in EtaBins_sorted:
        fig, ax = plt.subplots()  # Create a figure containing a single axes.


        # JEC applied in L1Ntuples --------------------------------------------
        Pt_appliedJEC = []
        AppiedJEC = []
        IEta_bin_tmp = hJEC_iEta_vs_Pt.GetXaxis().FindBin( IEta )
        hJEC_forIEta = hJEC_iEta_vs_Pt.ProjectionY("%s_ProjY_iEta%d" % (hJEC_iEta_vs_Pt.GetName(), IEta), IEta_bin_tmp, IEta_bin_tmp)
        if PrintLevel >= 20:
            print("hJEC_forIEta.GetNbinsX(): {}".format(hJEC_forIEta.GetNbinsX()))

        for iPtBin in range(1, hJEC_forIEta.GetNbinsX()+1):
            Pt_appliedJEC.append( hJEC_forIEta.GetXaxis().GetBinCenter(iPtBin) )
            AppiedJEC.append( hJEC_forIEta.GetBinContent(iPtBin) )

        if PrintLevel >= 20:
            sTmp = [ "%f: %f" % (Pt_tmp, AppiedJEC[idx_tmp])  for idx_tmp, Pt_tmp in enumerate(Pt_appliedJEC) ]
            print("\n\nsTmp {}: {}".format(IEta, sTmp))

        ax.plot( np.array(Pt_appliedJEC),  np.array(AppiedJEC), label='L1Ntuples pTCorrected/(pTRaw - PU)', color='None', marker='o', markerfacecolor='b', markeredgecolor='b', markersize=1.5)
        # ---------------------------------------------------------------------
        
        
        # Syed's SFsInDecimal number
        SFs = [ float(lut_JEC[IEta][lut_pTCompression[Pt]])  for Pt in PtBins]

        sTmp =""
        for idx, Pt in enumerate(PtBins):
            sTmp += "%s: %f,  " % (Pt, SFs[idx])
        if PrintLevel >= 20:
            print("SFs: {}".format(sTmp))
        ax.plot( np.array(PtBins),  np.array(SFs), label='SFs stored in JEC LUT', color='None', marker='o', markerfacecolor='None', markeredgecolor='r', markersize=1.9, alpha=0.4)


        ax.set_xlabel('L1 jet pT [GeV]')
        ax.set_ylabel('SF')
        #ax.set_title("Simple Plot")  # Add a title to the axes.
        ax.legend()  # Add a legend.

        ax.set_xscale('log', base=10)
        #plt.xscale('log', base=2)

        # Set axes limit
        yMax = max(SFs) * 1.10
        yMin = min(SFs) * 0.90
        '''
        if yMax > hJEC_forIEta.GetMaximum():
            yMax = hJEC_forIEta.GetMaximum() * 1.05
        for iPtBin in range(1, hJEC_forIEta.GetNbinsX()+1):
            if
        '''
        npAppiedJEC = np.array(AppiedJEC)
        npAppiedJEC = npAppiedJEC[ np.where(npAppiedJEC > 0.1) ]
        #print("IEta {}, npAppiedJEC {}".format(IEta, npAppiedJEC))
        if npAppiedJEC.size > 0:
            yMax = max(yMax, np.amax(npAppiedJEC) * 1.10 )
            yMin = min(yMin, np.amin(npAppiedJEC) * 0.90 )
        if yMax > 3: yMax = 3
        plt.ylim(yMin, yMax)
        plt.xlim(1, 300)

        #plt.show()
        plt.savefig('%s/checkJEC_iEta%d.png' % (sDir1, int(IEta)))
    # -------------------------------------------------------------------------------------
