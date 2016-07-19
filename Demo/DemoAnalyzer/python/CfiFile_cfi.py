import FWCore.ParameterSet.Config as cms
import math

infinity = float("inf")

demo = cms.EDAnalyzer('DemoAnalyzer',
           tracks=cms.untracked.string('generalTracks'),
           dEdxEstimator=cms.untracked.string('dedxHarmonic2'),
           genParticles=cms.untracked.string('genParticles'),
           isData=cms.untracked.bool(False),
           gainsConfig=cms.untracked.string('./Demo/DemoAnalyzer/data/Data13TeVGains_v2.root'),
           crossTalkInvAlgo=cms.untracked.bool(False),
           correctFEDSat=cms.untracked.bool(False),
           clusterCleaning=cms.untracked.bool(False),
           useScalefactors=cms.untracked.bool(False),
           useHSCPCalibration=cms.untracked.bool(False),
           useMyCalibration=cms.untracked.bool(False),
           nFoundHitsCut=cms.untracked.uint32(0),
           etaCut=cms.untracked.double(infinity),
           chi2Cut=cms.untracked.double(infinity),
           deltaRMatching=cms.untracked.double(100),
           dEdxCut=cms.untracked.double(-1),
           dZCut=cms.untracked.double(9999),
           deltaPTMatching=cms.untracked.double(9999),
           deltaPMatching=cms.untracked.double(9999)             
)

