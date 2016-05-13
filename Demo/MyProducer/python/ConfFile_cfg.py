import FWCore.ParameterSet.Config as cms

process = cms.Process("Producer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/cms/Tutorials/TWIKI_DATA/TTJets_8TeV_53X.root'
    )
)

process.load("Demo.MyProducer.CfiFile_cfi")
process.MyTracks.generalTracks='generalTracks'

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root')
)

  
process.p = cms.Path(process.MyTracks)

process.e = cms.EndPath(process.out)
