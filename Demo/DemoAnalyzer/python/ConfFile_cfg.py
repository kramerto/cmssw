import FWCore.ParameterSet.Config as cms
process = cms.Process("Demo")

# process.load("RecoVertex.BeamSpotProducer.BeamSpot_cff")
# process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
# process.TrackRefitter.src = 'generalTracks'
# 
# process.dedxHitInfo = cms.EDProducer("HSCPDeDxInfoProducer",
#     tracks = cms.InputTag("TrackRefitter"),
#     trajectoryTrackAssociation = cms.InputTag("TrackRefitter"),
# 
#     UseStrip  = cms.bool(True),
#     UsePixel  = cms.bool(True),
#     MeVperADCStrip = cms.double(3.61e-06*265),
#     MeVperADCPixel = cms.double(3.61e-06),
# 
#     UseCalibration = cms.bool(False),
#     calibrationPath = cms.string("file:Gains.root"),
#     ShapeTest = cms.bool(True),
# )


#initialize MessageLogger and output report
#process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('Demo')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#    limit = cms.untracked.int32(-1)
#)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
inputpath='file:/pnfs/desy.de/cms/tier2/store/data/Run2015D/JetHT/AOD/16Dec2015-v1/00000/'
outputpath='file:/nfs/dust/cms/user/kramerto/'
name='FEF28D68-A4AF-E511-A3FD-0025904C63F8'
root='.root'
ntuple='_NTuple'
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(inputpath+name+root),
#                 cms.untracked.vstring(
#             'file:/nfs/dust/cms/user/vormwald/HSCPgluino_M-1000_TuneCUETP8M1_13TeV-pythia8_00001.root'                            
#             'file:/nfs/dust/cms/user/vormwald/RunIIFall15DR76_TT_TuneCUETP8M1_13TeV-powheg-pythia8_AODSIM_25nsFlat10to25TSG_76X_mcRun2_asymptotic_v11_ext3-v1.root'                            
#             'file:/afs/cern.ch/cms/Tutorials/TWIKI_DATA/TTJets_8TeV_53X.root'
#             'file:/afs/desy.de/user/k/kramerto/CMSSW_8_0_6/src/SUSYBSMAnalysis/HSCP/test/MakeEDMtuples/HSCP.root'
#     )
#     firstEvent = cms.untracked.uint32(4553)
)
process.load("Demo.DemoAnalyzer.CfiFile_cfi")
process.demo.tracks='generalTracks'
process.demo.dEdxEstimator='dedxHarmonic2'
process.demo.genParticles='genParticles'
process.demo.isData=True
process.demo.gainsConfig='./Demo/DemoAnalyzer/data/Data13TeVGains_v2.root'
process.demo.crossTalkInvAlgo=False
process.demo.correctFEDSat=True
process.demo.clusterCleaning=True
process.demo.nFoundHitsCut=10
process.demo.etaCut=2.3
process.demo.chi2Cut=3
process.demo.deltaRMatching=0.03
process.demo.dEdxCut=0
process.demo.dZCut=0.1
# process.demo.deltaPTMatching=30
# process.demo.deltaPMatching=9999
# process.producer.generalTracks='generalTracks'

#process.dump=cms.EDAnalyzer('EventContentAnalyzer')

#process.Tracer = cms.Service("Tracer")

process.TFileService = cms.Service("TFileService",
#                                        fileName = cms.string('HSCPgluino_M-1000_TuneCUETP8M1_13TeV-pythia8_00001_NTuple.root')
                                       fileName = cms.string(outputpath+name+ntuple+root)
                                   )

#process.out = cms.OutputModule("PoolOutputModule",
 #   fileName = cms.untracked.string('myOutputFile.root')
    
#)

process.p = cms.Path(process.demo)
# process.p = cms.Path(process.offlineBeamSpot*process.TrackRefitter*process.dedxHitInfo*process.demo)

# process.e = cms.EndPath(process.out)
