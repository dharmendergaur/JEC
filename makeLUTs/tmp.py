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


for iEta_category, iEtaBinRange in IEta_Cat_forML.items():
    iEtaBins_i = range(iEtaBinRange[0], iEtaBinRange[-1]+1)

    for Pt_category, PtRange in Pt_Cat_forML.items():
        PtRangeMin = PtRange[0]
        PtRangeMax = PtRange[1]
        
        data_all_iEtaBins = data_all[
            (data_all[sL1JetTowerIEtaAbs] >= iEtaBinRange[0]) & 
            (data_all[sL1JetTowerIEtaAbs] <= iEtaBinRange[-1]) &
            (data_all[sL1JetEt] >= PtRangeMin) &
            (data_all[sL1JetEt] <  PtRangeMax)
        ][varsOfInterest]
        if printLevel >= 0:
            print("\niEta_category {}, iEtaBinRange {}, data_all_iEtaBins.describe(): \n{}".format(
                iEta_category, iEtaBinRange, data_all_iEtaBins.describe()))


        X = data_all_iEtaBins[train_vars]
        y = data_all_iEtaBins[target_var]

        xgb_rg = train_MLModel_wHyperopt(X, y)


        data_SFs_i = prepareDataframeForSFs(iEtaBins_i)
        data_SFs_i[sL1JetEt_forML_predict] = xgb_rg.predict(data_SFs_i[train_vars])
        data_SFs_i[sL1JetEt_predict]       = transform_back_JetEt_fromML( data_SFs_i[sL1JetEt_forML_predict] )
        data_SFs_i[sSF]                    = data_SFs_i[sL1JetEt_predict] / data_SFs_i[sL1JetEt]
        if printLevel >= 11:
            print("iEtaBins_i: {}".format(iEtaBins_i))
            print("data_SFs_i: {}".format(data_SFs_i.describe()))

        if data_SFs is None:
            data_SFs = data_SFs_i
        else:
            data_SFs = pd.concat([data_SFs, data_SFs_i])    

