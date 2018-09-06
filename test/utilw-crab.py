from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'testBTAG-17_TT4'
config.General.workArea = '/scratch/lgiannini/bTAG_retrain/CRAB/'
config.General.transferOutputs = True
config.General.transferLogs = True


config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'miniAOD17CRAB_config.py'


config.Data.inputDataset = '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'
#config.Data.inputDataset = '/QCD_Pt_120to170_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/%s/testBTAG-17' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'CRAB3_tutorial_May2015_MC_analysis'
config.Data.totalUnits = config.Data.unitsPerJob * 10


config.Site.storageSite = "T2_IT_Pisa"


##o
#from CRABAPI.RawCommand import crabCommand
#crabCommand('submit', config = config)

##o
#crab submit -c crabConfig_data_slimMiniAOD.py


##check what you do
#crab submit -c crabConfig_tutorial_Data_analysis.py --dryrun
