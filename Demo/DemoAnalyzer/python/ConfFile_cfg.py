import FWCore.ParameterSet.Config as cms
process = cms.Process("Demo")

process.load("RecoVertex.BeamSpotProducer.BeamSpot_cff")
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.TrackRefitter.src = 'generalTracks'

process.dedxHitInfo = cms.EDProducer("HSCPDeDxInfoProducer",
    tracks = cms.InputTag("TrackRefitter"),
    trajectoryTrackAssociation = cms.InputTag("TrackRefitter"),

    UseStrip  = cms.bool(True),
    UsePixel  = cms.bool(True),
    MeVperADCStrip = cms.double(3.61e-06*265),
    MeVperADCPixel = cms.double(3.61e-06),

    UseCalibration = cms.bool(False),
    calibrationPath = cms.string("file:Gains.root"),
    ShapeTest = cms.bool(True),
)


#initialize MessageLogger and output report
#process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('Demo')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#    limit = cms.untracked.int32(-1)
#)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#             'file:/nfs/dust/cms/user/vormwald/RunIIFall15DR76_TT_TuneCUETP8M1_13TeV-powheg-pythia8_GEN-SIM-RAW_25nsFlat10to25TSG_76X_mcRun2_asymptotic_v11_ext3-v1.root '                            
            'file:/nfs/dust/cms/user/vormwald/RunIIFall15DR76_TT_TuneCUETP8M1_13TeV-powheg-pythia8_AODSIM_25nsFlat10to25TSG_76X_mcRun2_asymptotic_v11_ext3-v1.root'                            
#             'file:/afs/cern.ch/cms/Tutorials/TWIKI_DATA/TTJets_8TeV_53X.root'
#             'file:myOutputFile.root'
    )
)

process.load("Demo.DemoAnalyzer.CfiFile_cfi")
process.demo.minTracks=0
process.demo.tracks='generalTracks'
process.demo.dEdxEstimator='dedxHarmonic2'
process.demo.nFoundHitsCut=10
process.demo.etaCut=1.4
process.demo.chi2Cut=3
process.demo.deltaRCut=0.022
process.demo.dEdxCut=0
process.demo.dZCut=0.145
# process.producer.generalTracks='generalTracks'

#process.dump=cms.EDAnalyzer('EventContentAnalyzer')

#process.Tracer = cms.Service("Tracer")

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('histodemo.root')
                                   )

#process.out = cms.OutputModule("PoolOutputModule",
 #   fileName = cms.untracked.string('myOutputFile.root')
    
#)

process.p = cms.Path(process.offlineBeamSpot+process.TrackRefitter+process.dedxHitInfo+process.demo)

# process.e = cms.EndPath(process.out)
