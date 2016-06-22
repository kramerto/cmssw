import FWCore.ParameterSet.Config as cms

process = cms.Process("dEdx")
process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.Services_cff')

process.GlobalTag.globaltag = 'START72_V1::All'



process.load("SUSYBSMAnalysis.HSCP.HSCParticleProducer_cff") 

# 
# process.load("RecoTracker.MeasurementDet.MeasurementTrackerEventProducer_cfi")
# process.load("RecoVertex.BeamSpotProducer.BeamSpot_cff")
# process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
# 
# ####################################################################################                                                                                                                         
# #   HIT-DEDX Information                                                                                                                                                                                     
# ####################################################################################                                                                                                                         
# 
# process.dedxHitInfo               = cms.EDProducer("HSCPDeDxInfoProducer",
#     tracks                     = cms.InputTag("TrackRefitter"),
#     trajectoryTrackAssociation = cms.InputTag("TrackRefitter"),
# 
#     Reccord            = cms.untracked.string("SiStripDeDxMip_3D_Rcd"),
#     Formula            = cms.untracked.uint32(0),
#     ProbabilityMode    = cms.untracked.string("Accumulation"),
# 
#     UseStrip           = cms.bool(True),
#     UsePixel           = cms.bool(True),
#     MeVperADCStrip     = cms.double(3.61e-06*265),
#     MeVperADCPixel     = cms.double(3.61e-06),
# 
#     UseCalibration     = cms.bool(False),
#     calibrationPath    = cms.string("file:Gains.root"),
#     ShapeTest          = cms.bool(True),
# )
# 
# process.load("Demo.DemoAnalyzer.CfiFile_cfi")
# 
# process.maxEvents = cms.untracked.PSet(
#     input = cms.untracked.int32(100)
# )
#   
# # Input source
# process.source = cms.Source("PoolSource",
#     fileNames = cms.untracked.vstring(
# #                                     'file:/nfs/dust/cms/user/vormwald/TOP-RunIIFall15DR76-00002.root'
#                                     'file:/nfs/dust/cms/user/vormwald/RunIIFall15DR76_TT_TuneCUETP8M1_13TeV-powheg-pythia8_AODSIM_25nsFlat10to25TSG_76X_mcRun2_asymptotic_v11_ext3-v1.root'                            
#                                       ),
#     secondaryFileNames = cms.untracked.vstring()
# )
#    
# process.out = cms.OutputModule("PoolOutputModule",
#     fileName = cms.untracked.string('myOutputFile.root')
# )
#   
# # HSCParticleProducerSeq = cms.Sequence(offlineBeamSpot + MeasurementTrackerEvent + TrackRefitter  + dedxHitInfo)
# process.p =cms.Path(process.offlineBeamSpot + process.MeasurementTrackerEvent + process.TrackRefitter  + process.dedxHitInfo)
# process.e = cms.EndPath(process.out)
