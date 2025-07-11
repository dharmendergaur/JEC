# L1T Jet layer 2 calibration

Script to make L1T JEC LUTs

## Step 1: Produce L1T Nano Files

- Recipe to make L1T Nanos can be found at [here](https://github.com/cms-sw/cmssw/tree/master/DPGAnalysis/L1TNanoAOD). Choose the latest CMSSW release and l1t-integration tag.

- **Dummy `cmsDriver` command** for running over `Muon1/RAW-RECO`(_change the input file if not available_). Make sure you have required caloparams (`L1TSettingsToCaloParams_2025_baseline_v2_iET`) correctly defined in `L1Trigger/Configuration/python/customiseSettings.py`, if not, refer [here](https://indico.cern.ch/event/1504117/) for instructions.

        cmsDriver.py customL1toNANO --conditions 140X_dataRun3_v20 -s RAW2DIGI,L1Reco,NANO:@PHYS+@L1DPG --datatier NANOAOD --eventcontent NANOAOD --data --process customl1nano --scenario pp --era Run3 --customise Configuration/DataProcessing/RecoTLR.customisePostEra_Run3 --customise_unsch L1Trigger/Configuration/customiseReEmul.L1TReEmulFromRAW --customise_unsch L1Trigger/Configuration/customiseSettings.L1TSettingsToCaloParams_2025_baseline_v2_iET --secondfilein /store/data/Run2024I/Muon0/RAW/v1/000/386/505/00000/46111224-ec19-4680-b777-4a175e09659e.root --filein /store/data/Run2024I/Muon1/RAW-RECO/ZMu-PromptReco-v1/000/386/679/00000/429dca6f-1658-46fe-a9f1-b45c16530332.root -n -1 --no_exec

- Sample files are available in `sample_l1nano_script/` with sample crab submit file.

## Step 2: Produce BDT input .csv file

The produced L1Nano files are read and BDT input in `.csv` files are produced by `BDTInputProducer/script.py`.

Input parameters read by default:

- `--l1nano`: ../sample_root_files/Nano.root.root
- `--output`: Nano_out
- `--l1MatchOffline`: True
- `--l1nanoPhiRing`

## Step 3: BDT trainning

1.  Calibration JEC SF computation script is `setup/calculate_L1JetSFs_usingBD.ipynb`.
2.  Run for 4 different setups using the correct --ipFile flag:

        python3 calculate_L1JetSFs_usingBDT.py --PhiRing     --l1MatchGen --MLTarget GenEt --ipFile _csv_file_

        python3 calculate_L1JetSFs_usingBDT.py --PhiRing     --l1MatchGen --MLTarget logGenEt --ipFile _csv_file_

        python3 calculate_L1JetSFs_usingBDT.py --PhiRing     --l1MatchGen --MLTarget GenEtByL1Et --ipFile _csv_file_

        python3 calculate_L1JetSFs_usingBDT.py --PhiRing     --l1MatchGen --MLTarget logGenEtByL1Et --ipFile _csv_file_

    Modify `--fracOfDataToUse` accordingly.

3.  Run `setup/compare_JEC_SFs.ipynb` to get final `JEC_SF.csv` file.

4.  Run `setup/check_L1JetSFs.ipynb` to make quality control plots.

**Note:** Outputs will be written in the `data/` directory.

## Step 4: Make LUTs

1. cd makeLUTs

2. Run `python3 updateSFPtEtaBins.py`

3. Run `check_L1TJetSFs.ipynb` to make performance plots using the new JEC LUTs

# git commit:

git remote add origin git@github.com:dharmendergaur/JEC.git
git push -u origin master
