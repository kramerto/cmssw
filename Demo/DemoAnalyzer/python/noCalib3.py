import FWCore.ParameterSet.Config as cms
process = cms.Process("Demo")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2500000) )

import os

basedir = "/pnfs/desy.de/cms/tier2/store/data/Run2016B/SingleMuon/AOD/PromptReco-v2"
content = os.listdir(basedir)
listoffiles = []

for runnumber1 in content:
    basedirrun1 = basedir+'/'+runnumber1
    runnumbers2 = os.listdir(basedirrun1)
#     for runnumber2 in runnumbers2:
#         basedirrun2= basedirrun1+'/'+runnumber2
    basedirrun2= basedirrun1+'/275'
    runnumbers3 = os.listdir(basedirrun2)
    for runnumber3 in runnumbers3:
        basedirrun3= basedirrun2+'/'+runnumber3
        versions= os.listdir(basedirrun3)
        for version in versions:
            basedirversion= basedirrun3+'/'+version
            files = os.listdir(basedirversion)
            for f in files:
                    #print(basedirversion+'/'+f)
                listoffiles.append('file:'+basedirversion+'/'+f)

# print(listoffiles) 

outputpath='file:/nfs/dust/cms/user/kramerto/'
name='Run2016B_SingleMuon_noCalibration_part3'
root='.root'
ntuple='_NTuple'
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(listoffiles),
                            skipEvents=cms.untracked.uint32(5000000))
  
process.load("Demo.DemoAnalyzer.CfiFile_cfi")
process.demo.tracks='generalTracks'
process.demo.dEdxEstimator='dedxHarmonic2'
process.demo.genParticles='genParticles'
process.demo.isData=True
process.demo.gainsConfig='./Demo/DemoAnalyzer/data/Data13TeVGains_v2.root'
process.demo.crossTalkInvAlgo=False
process.demo.correctFEDSat=False
process.demo.clusterCleaning=False
process.demo.useScalefactors=False
process.demo.useHSCPCalibration=False
process.demo.useMyCalibration=False
process.demo.nFoundHitsCut=10
# process.demo.etaCut=1.4
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
                                       fileName = cms.string(outputpath+name+ntuple+root)
                                   )

#process.out = cms.OutputModule("PoolOutputModule",
 #   fileName = cms.untracked.string('myOutputFile.root')
    
#)

process.p = cms.Path(process.demo)
# process.p = cms.Path(process.offlineBeamSpot*process.TrackRefitter*process.dedxHitInfo*process.demo)

# process.e = cms.EndPath(process.out)
