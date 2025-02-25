# L1T Jet layer 2 calibration

Script to make L1T JEC LUTs

## Step 1: Produce L1T Nano Files

Recipe to make L1T Nanos can be found at [here](https://github.com/cms-sw/cmssw/tree/master/DPGAnalysis/L1TNanoAOD). Choose the latest CMSSW release and l1t-integration tag.

Sample files are available in `sample_l1nano_script` with crab submit file.

## Step 2: Produce BDT input .csv file

The produced L1Nano files are read and BDT input in `.csv` files are produced by `BDTInputProducer/script.py`.

Input parameters read by default:

- `--l1nano`: ../sample_root_files/Nano.root.root
- `--output`: Nano_out
- `--l1MatchOffline`: True
- `--l1nanoPhiRing`

## Step 3: BDT trainning

1.  Calibration JEC SF computation script is `setup/calculate_L1JetSFs_usingBD.ipynb`.
2.  Run for 4 different setups:

        python3 calculate_L1JetSFs_usingBDT.py --PhiRing     --l1MatchGen --MLTarget GenEt

        python3 calculate_L1JetSFs_usingBDT.py --PhiRing     --l1MatchGen --MLTarget logGenEt

        python3 calculate_L1JetSFs_usingBDT.py --PhiRing     --l1MatchGen --MLTarget GenEtByL1Et

        python3 calculate_L1JetSFs_usingBDT.py --PhiRing     --l1MatchGen --MLTarget logGenEtByL1Et

Modify `--fracOfDataToUse` accordingly.

3. Run `setup/compare_JEC_SFs.ipynb` to get final `JEC_SF.csv` file.

4. Run `setup/check_L1JetSFs.ipynb` to make quality control plots.

Notes: Outputs will be written in the `data/` directory.

## Step 4: Make LUTs

1. cd makeLUTs

2. Run `python3 updateSFPtEtaBins.py`

3. Run `check_L1TJetSFs.ipynb` to make performance plots using the new JEC LUTs

# git commit:

git remote add origin git@github.com:dharmendergaur/JEC.git
git push -u origin master
