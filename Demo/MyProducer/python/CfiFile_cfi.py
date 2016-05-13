import FWCore.ParameterSet.Config as cms

MyTracks = cms.EDProducer('MyProducer',
            generalTracks=cms.InputTag('generalTracks')
)
