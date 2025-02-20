def run():
    start_time = time.time()
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
    Emu_br = BranchCollection(branch_data["emu"], "L1EmulJet_", "nL1EmulJet")
    Unp_br = BranchCollection(branch_data["unp"], "L1Jet_", "nL1Jet")
    Vtx_br = BranchCollection(branch_data["vtx"], "PV_", "PV_npvsGood")
    Evt_br = BranchCollection(branch_data["evt"], "", "")

for iEvent in range(Ntot):
    # for iEvent in range(500):

        l1JetRef_br = None
        nRefJets    = 0
        et_RefJets  = None
        eta_RefJets = None
        phi_RefJets = None
        NJets= Jet_br['nObjects'][iEvent]
        nOffJets= Jet_br['nObjects'][iEvent]
        et_RefJets=Jet_br["pt"][iEvent]
        eta_RefJets=Jet_br["eta"][iEvent]
        phi_RefJets=Jet_br["phi"][iEvent]

        # print("Event: %d, NJets: %d" % (iEvent, NJets))

        nTotalEvents_byChains[0] += 1
        hStat.Fill(0)

        #Analyze (GEN.nVtx == 0) events from SinglePhoton_EpsilonPU sample to trouble-shoot high SFs in iEta 28 ----
        # if isMC and useCutGenNVtxEq0:
        #     if Gen_br.nVtx > 0: continue
        # ----------------------------------------------------------------------------------------------------------
            

    
        # if not isMC and len(GoldenJSONForData_list) > 0:
        #     if not passGoldenJSON(goldenJSON, int(Evt_br['run'][iEvent]), int(Evt_br['luminosityBlock'][iEvent])):
        #         #print(f"Run:LS:Event:  %d:%d:%d   fails GoldenJSON " %(int(Evt_br['run']), int(Evt_br['luminosityBlock']), int(Evt_br['event']))); sys.stdout.flush();
        #         continue


        dataEra = ''
        if not isMC:
            for Era_, eraRunRange_ in dataErasRunRange.items():
                if int(Evt_br['run'][iEvent]) >= eraRunRange_[0] and int(Evt_br['run'][iEvent]) <= eraRunRange_[1]:
                    dataEra = Era_
                    break
        
        # print("dataEra: %s" % (dataEra))

        hStat.Fill(1)

        
        hCaloTowers_iEta_vs_iPhi = None
        hCaloTTs_iEta_vs_iPhi    = None

        if sFInEventsToRun: # and runMode in ['trbshtPhiRingPUS']:
            sRunLSEvent = "%s:%s:%s" % (int(Evt_br['run'][iEvent]), int(Evt_br['luminosityBlock'][iEvent]), int(Evt_br['event'][iEvent]))
            sRunLSEvent_toUse = sRunLSEvent.replace(":", "_") 
            if sRunLSEvent not in eventsToRun_list: continue
            print("Ruuning on selected event %s" % (sRunLSEvent))

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
        l1MatchOffline= True
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
                # PF CHS jets
                nRefJets     = Jet_br['nObjects'][iEvent]
                et_RefJets   = Jet_br['pt'][iEvent]
                eta_RefJets  = Jet_br['eta'][iEvent]
                phi_RefJets  = Jet_br['phi'][iEvent]                    
        elif l1MatchGen:
            nRefJets     = Gen_br.nJet[iEvent]
            et_RefJets   = Gen_br.jetPt[iEvent]
            eta_RefJets  = Gen_br.jetEta[iEvent]
            phi_RefJets  = Gen_br.jetPhi[iEvent]

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

        if not isMC and len(HLT_Triggers_Required) > 0:
            nOffMuons_passingTrigThsh = [0] * len(TrigThshs_OffMuPt)
            if PrintLevel >= 3:
                print(f"nOffMuons_passingTrigThsh_0: {nOffMuons_passingTrigThsh}")
            for iMu in range(nOffMuons):
                if not Muon_br['tightId'][iEvent][iMu]: continue

                for iTrigThsh in range(len(TrigThshs_OffMuPt)):
                    if Muon_br['pt'][iEvent][iMu] > TrigThshs_OffMuPt[iTrigThsh]:
                        nOffMuons_passingTrigThsh[iTrigThsh] += 1

            passingTrigThshs = True
            for iTrigThsh in range(len(TrigThshs_OffMuPt)):
                if nOffMuons_passingTrigThsh[iTrigThsh] == 0:
                    passingTrigThshs = False
                    break
                
            if PrintLevel >= 3:
                print(f"nOffMuons_passingTrigThsh: {nOffMuons_passingTrigThsh},   passingTrigThshs: {passingTrigThshs}")

            if not passingTrigThshs: continue

            hStat.Fill(3)
        
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
        for iOff in range(NJets):

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

                if any(reject_if):
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
            # elif l1MatchGen:
            #     data_dict['GenJetEt']             = vOff.Pt()
            #     data_dict['nVertexGen']           = int(Gen_br.nVtx)
            #     data_dict['nMeanPUGen']           = int(Gen_br.nMeanPU)
            #     data_dict['matchedPFJetEtCorr']   = off_vec_matc
        
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

                    if etaCat[src][algo] == 'None':
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

            if VERBOSE or PrintLevel >= 1:
                sTmp = "matchedEmuIdx: "
                for src in ['unp','emu']:
                    for algo in ['PUS','noPUS','Raw','RawPUS']:
                        sTmp += "  %s_%s %d" % (src,algo,matchedEmuIdx[src][algo])
                print("        {}".format(sTmp))

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
                l1jet_idx = matchedEmuIdx['emu']['PUS']
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
                    jetIEta_tmp    = convert_jetIEta_to_jetTowerIEta(l1jet_br['hwEta'][iEvent][l1jet_idx])         #revisit
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
                    jetIEta = convert_jetIEta_to_jetTowerIEta( 2*l1jet_br['eta'][iEvent][l1jet_idx])
                    jetIPhi = convert_jetIPhi_to_jetTowerIPhi( 2*l1jet_br['phi'][iEvent][l1jet_idx] )
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

        print("\n\n nTotalEvents_byChains[iEvent {}]: {} ".format(0, nTotalEvents_byChains[0]))
        hnTotalEvents.SetBinContent(1, nTotalEvents_byChains[0])
        # ## End loop: for iEvent in range(len(chains['Unp'])):

    print('\nFinished loop over Events')
    print(f"Processing completed in {time.time() - start_time:.2f} seconds")
    if runMode in ['makeInputForML']: 
        fOut_MLInputs.close()

    out_file = R.TFile(out_file_str,'recreate')
    out_file.cd()
