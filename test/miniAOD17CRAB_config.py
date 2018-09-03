import FWCore.ParameterSet.Config as cms
process = cms.Process("Pippo")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.load('Configuration/StandardSequences/GeometryDB_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/EventContent/EventContent_cff')
process.GlobalTag.globaltag = '94X_mc2017_realistic_v14'
process.options = cms.untracked.PSet()

#process.options.numberOfThreads = cms.untracked.uint32(1)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(

        "root://xrootd-cms.infn.it//store/mc/RunIIFall17MiniAODv2/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/30000/0055C65C-E558-E811-AB0E-008CFA582BF4.root"
       
    )#,
    #skipEvents=cms.untracked.uint32(36590)
)


process.analyzer1 = cms.EDAnalyzer('Analyzer_MINIAOD_new')

outfile='example_x.root'

process.TFileService = cms.Service("TFileService",
   fileName = cms.string(outfile) )


process.p = cms.Path(process.analyzer1)
 
