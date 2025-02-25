from CRABClient.UserUtilities import config

username='ddharmen'
lumi_mask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions24/Cert_Collisions2024_378981_386951_Golden.json'
output_dataset_tag = 'jetMET24I'
input_dataset = '/JetMET0/Run2024I-PromptReco-v1/MINIAOD'
pset_name = 'nano.py'
output_files = ['nano.root']
file_output_path = f'/store/user/{username}/JECs2025/JETMET24I/miniaod/baseline'
config = config()

# General settings
config.General.requestName = output_dataset_tag 
config.General.transferLogs = True
config.General.transferOutputs = True


config.JobType.pluginName = 'Analysis'
config.JobType.psetName = pset_name
config.JobType.outputFiles = output_files
config.JobType.maxMemoryMB = 2500

# Data settings
config.Data.inputDataset = input_dataset
config.Data.inputDBS = 'global'

config.Data.useParent = True
config.Data.partialDataset = True

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1

config.Data.allowNonValidInputDataset = True

config.Data.outputDatasetTag = output_dataset_tag
config.Data.outLFNDirBase = file_output_path

config.Data.lumiMask = lumi_mask
config.Site.storageSite = 'T2_CH_CERN'
