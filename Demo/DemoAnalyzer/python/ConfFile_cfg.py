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
inputpath='file:/pnfs/desy.de/cms/tier2/store/data/Run2016B/SingleMuon/AOD/PromptReco-v2/000/274/284/00000/'
outputpath='file:/nfs/dust/cms/user/kramerto/'
name='Run2015D_16Dec2015_00000'
name1='08DF53A6-6729-E611-84D7-02163E0119A6'
# name2='A644739F-C0AF-E511-86AB-C4346BC8C638'
# name3='A6BF018A-AEAF-E511-89DA-0CC47A78A3F8'
# name4='BEF2EA39-A2AF-E511-9AD2-0025907B50FC'
# name5='C0660B15-A9AF-E511-A511-002590D0B01C'
# name6='CC44563A-AEAF-E511-9D89-C4346BC78D10'
# name7='D0B4D81A-AFAF-E511-B30E-002618943939'
# name8='D0CE3CE8-9FAF-E511-BD9B-00259073E466'
# name9='D4DC2895-AEAF-E511-B6DB-00261894390E'
# name10='DAC1908D-AEAF-E511-AB4E-0CC47A4D76AA'
# name11='DE841E8E-AEAF-E511-A7B6-0CC47A4D768C'
# name12='E668AA98-AEAF-E511-B5A9-002590593902'
# name13='EACEA337-A2AF-E511-A0EC-0025907B4FD6'
# name14='F2288E7C-EDAF-E511-B836-0CC47A4D7674'
# name15='F25A3F8D-AEAF-E511-82AF-0CC47A4C8EC8'
# name16='F893AE8C-A6AF-E511-8F57-0025907B5002'
# name17='FEF28D68-A4AF-E511-A3FD-0025904C63F8'
root='.root'
ntuple='_NTuple'
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
                                      inputpath+name1+root,
#                                       inputpath+name2+root,
#                                       inputpath+name3+root,
#                                       inputpath+name4+root,
#                                       inputpath+name5+root,
#                                       inputpath+name6+root,
#                                       inputpath+name7+root,
#                                       inputpath+name8+root,
#                                       inputpath+name9+root,
#                                       inputpath+name10+root,
#                                       inputpath+name11+root,
#                                       inputpath+name12+root,
#                                       inputpath+name13+root,
#                                       inputpath+name14+root,
#                                       inputpath+name15+root,
#                                       inputpath+name16+root,
#                                       inputpath+name17+root
                                      )
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
process.demo.correctFEDSat=False
process.demo.clusterCleaning=False
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
