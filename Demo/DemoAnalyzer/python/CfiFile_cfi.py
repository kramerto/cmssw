import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('DemoAnalyzer',
           minTracks=cms.untracked.uint32(0),
           tracks=cms.untracked.string('generalTracks'),
           dEdxEstimator=cms.untracked.string('dedxHarmonic2')
)

