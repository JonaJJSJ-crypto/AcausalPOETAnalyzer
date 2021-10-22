import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
import os
import sys

relBase = os.environ['CMSSW_BASE']

process = cms.Process('MERGE')

#process.options = cms.untracked.PSet(
#    wantSummary = cms.untracked.bool(True)
#)
#process.load('FWCore.MessageService.MessageLogger_cfi')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#    reportEvery = cms.untracked.int32(1), # every!
#    limit = cms.untracked.int32(-1)     # no limit!
#    )
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000 # only report every 10th event start
#process.MessageLogger.cerr_stats.threshold = 'WARNING' 

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = "WARNING"
process.MessageLogger.categories.append("MERGE")
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit=cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))


files = FileUtils.loadListFromFile("data/CMSPrivate_MonteCarlo_LWSM200DnR_unmerged_file_index.txt")
process.source = cms.Source(
    "PoolSource", fileNames=cms.untracked.vstring(*files))


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


# define the PoolOutputModule
process.poolOutput = cms.OutputModule('PoolOutputModule',
                                      fileName = cms.untracked.string('mymerged.root'),
                                      outputCommands = cms.untracked.vstring('keep *'),
                                      maxSize = cms.untracked.int32(2500000)
)

#process.poolOutput.outputCommands.append('keep FEDRawDataCollection_source_*_*')
#process.poolOutput.outputCommands.append('keep L1GlobalTriggerObjectMapRecord_hltL1GtObjectMap_*_*')
#process.poolOutput.outputCommands.append('keep triggerTriggerEvent_hltTriggerSummaryAOD_*_*')
#process.poolOutput.outputCommands.append('keep edmTriggerResults_TriggerResults_*_*')
process.output = cms.EndPath(process.poolOutput)
