import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

#initialize MessageLogger and output report
#process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('Demo')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#    limit = cms.untracked.int32(-1)
#)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/cms/Tutorials/TWIKI_DATA/TTJets_8TeV_53X.root'
#         'file:myOutputFile.root'
    )
)

process.load("Demo.DemoAnalyzer.CfiFile_cfi")
process.demo.minTracks=0
process.demo.tracks='generalTracks'
process.demo.dEdxEstimator='dedxHarmonic2'
process.demo.nFoundHitsCut=0
process.demo.etaCut=1.4
process.demo.chi2Cut=3
process.demo.deltaRCut=0.01
process.demo.dEdxCut=0
# process.producer.generalTracks='generalTracks'

#process.dump=cms.EDAnalyzer('EventContentAnalyzer')

#process.Tracer = cms.Service("Tracer")

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('histodemo.root')
                                   )

#process.out = cms.OutputModule("PoolOutputModule",
 #   fileName = cms.untracked.string('myOutputFile.root')
    
#)

process.p = cms.Path(process.demo)

# process.e = cms.EndPath(process.out)
