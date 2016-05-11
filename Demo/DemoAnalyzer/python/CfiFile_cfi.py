import FWCore.ParameterSet.Config as cms
import math

infinity = float("inf")

demo = cms.EDAnalyzer('DemoAnalyzer',
           minTracks=cms.untracked.uint32(0),
           tracks=cms.untracked.string('generalTracks'),
           dEdxEstimator=cms.untracked.string('dedxHarmonic2'),
           nFoundHitsCut=cms.untracked.uint32(0),
           etaCut=cms.untracked.double(infinity),
           chi2Cut=cms.untracked.double(infinity),
           deltaRCut=cms.untracked.double(100),
           dEdxCut=cms.untracked.double(-1)
)

