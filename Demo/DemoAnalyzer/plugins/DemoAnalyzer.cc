// -*- C++ -*-
//
// Package:    Demo/DemoAnalyzer
// Class:      DemoAnalyzer
//
/**\class DemoAnalyzer DemoAnalyzer.cc Demo/DemoAnalyzer/plugins/DemoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Tobias Robert Jakob Kramer
//         Created:  Wed, 27 Apr 2016 15:04:51 GMT
//
//
//
//
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "TBranch.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "Geometry/CommonDetUnit/interface/GeomDetType.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TLorentzVector.h"

//
// class declaration
//

class DemoAnalyzer : public edm::EDAnalyzer {
public:
  explicit DemoAnalyzer(const edm::ParameterSet &);
  ~DemoAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event &, const edm::EventSetup &) override;
  virtual void endJob() override;

  // virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  // virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  // virtual void beginLuminosityBlock(edm::LuminosityBlock const&,
  // edm::EventSetup const&) override;
  // virtual void endLuminosityBlock(edm::LuminosityBlock const&,
  // edm::EventSetup const&) override;

  // ----------member data ---------------------------
  unsigned int minTracks_;
  std::string tracks_;
  std::string dEdxEstimator_;
  TTree *tree_;
  int nTracks_;
  std::vector<double> dEdx_;
  std::vector<double> nHits_;
  std::vector<double> charge_;
  std::vector<double> nLostHits_;
  std::vector<double> nFoundHits_;
  std::vector<double> nDEdxMeasurements_;
  std::vector<double> chi2_;
  std::vector<double> ndof_;
  std::vector<double> p_;
  std::vector<int> particleID_;
  std::vector<double> deltaR_;
  //  std::vector<int> particleID_2;
  //  std::vector<double> deltaR_2;
  std::vector<double> deltaEta_;
  std::vector<double> deltaPhi_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet &iConfig)
    : minTracks_(iConfig.getUntrackedParameter<unsigned int>("minTracks", 0)),
      tracks_(iConfig.getUntrackedParameter<std::string>("tracks",
                                                         "generalTracks")),
      dEdxEstimator_(iConfig.getUntrackedParameter<std::string>(
          "dEdxEstimator", "dedxHarmonic2")) {

  // now do what ever initialization is needed
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("Tree", "Tree");
  tree_->Branch("dEdx", &dEdx_);
  tree_->Branch("nTracks", &nTracks_);
  tree_->Branch("nHits", &nHits_);
  tree_->Branch("charge", &charge_);
  tree_->Branch("nLostHits", &nLostHits_);
  tree_->Branch("nFoundHits", &nFoundHits_);
  tree_->Branch("nDEdxMeasurements", &nDEdxMeasurements_);
  tree_->Branch("chi2", &chi2_);
  tree_->Branch("ndof", &ndof_);
  tree_->Branch("p", &p_);
  tree_->Branch("particleID", &particleID_);
  tree_->Branch("deltaR", &deltaR_);
  tree_->Branch("deltaEta", &deltaEta_);
  tree_->Branch("deltaPhi", &deltaPhi_);

  //  tree_->Branch("particleID2", &particleID_2);
  //  tree_->Branch("deltaR2", &deltaR_2);
}

DemoAnalyzer::~DemoAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void DemoAnalyzer::analyze(const edm::Event &iEvent,
                           const edm::EventSetup &iSetup) {
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel("MyTracks", tracks);
  edm::Handle<edm::ValueMap<reco::DeDxData>> dEdxTrackHandle;
  iEvent.getByLabel(dEdxEstimator_, dEdxTrackHandle);
  const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxTrackHandle.product();
  nTracks_ = tracks->size();
  dEdx_.clear();
  nHits_.clear();
  nLostHits_.clear();
  nFoundHits_.clear();
  charge_.clear();
  chi2_.clear();
  nDEdxMeasurements_.clear();
  ndof_.clear();
  p_.clear();
  particleID_.clear();
  deltaR_.clear();
  //  particleID_2.clear();
  //  deltaR_2.clear();
  deltaEta_.clear();
  deltaPhi_.clear();
  for (unsigned int i = 0; i < tracks->size(); i++) {
    nHits_.push_back(tracks->at(i).found() + tracks->at(i).lost());
    nLostHits_.push_back(tracks->at(i).lost());
    nFoundHits_.push_back(tracks->at(i).found());
    charge_.push_back(tracks->at(i).charge());
    chi2_.push_back(tracks->at(i).chi2());
    ndof_.push_back(tracks->at(i).ndof());
    reco::TrackRef track = reco::TrackRef(tracks, i);
    // Track momentum is given by:
    // track->p();
    // You can access dE/dx Estimation of your track with:
    // dEdxTrack[track].numberOfSaturatedMeasurements();
    dEdx_.push_back(dEdxTrack[track].dEdx());
    p_.push_back(track->p());
    nDEdxMeasurements_.push_back(dEdxTrack[track].numberOfMeasurements());
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel("genParticles", genParticles);
    double deltaR = 100;
    double deltaEta = 100;
    double deltaPhi = 100;
    int pID = 0;
    for (unsigned int j = 0; j < genParticles->size(); ++j) {
      const reco::GenParticle &p = (*genParticles)[j];
      if ((p.status() == 1) && (p.charge() != 0)) {
        if (std::sqrt(std::pow((tracks->at(i).eta() - p.eta()), 2) +
                      std::pow(fmod(std::abs(tracks->at(i).phi() - p.phi()),
                                    (4 * atan(1))),
                               2)) < deltaR) {
          deltaR =
              std::sqrt(std::pow((tracks->at(i).eta() - p.eta()), 2) +
                        std::pow(fmod(std::abs(tracks->at(i).phi() - p.phi()),
                                      (4 * atan(1))),
                                 2));
          deltaEta = std::abs(tracks->at(i).eta() - p.eta());
          deltaPhi =
              fmod(std::abs(tracks->at(i).phi() - p.phi()), (4 * atan(1)));

          pID = p.pdgId();
        }
        //    int id = p.pdgId();
        //    int st = p.status();
        //    const reco::Candidate *mom = p.mother();
        //    double pt = p.pt(), eta = p.eta(), phi = p.phi(), mass =
        //    p.mass();
        //    double vx = p.vx(), vy = p.vy(), vz = p.vz();
        //    int charge = p.charge();
        //    int n = p.numberOfDaughters();
        //    for (int j = 0; j < n; ++j) {
        //      const reco::Candidate *d = p.daughter(j);
        //      int dauId = d->pdgId();
        //    }
      }
    }
    deltaR_.push_back(deltaR);
    deltaEta_.push_back(deltaEta);
    deltaPhi_.push_back(deltaPhi);
    if (deltaR < 0.02)
      particleID_.push_back(pID);
    else
      particleID_.push_back(0);
    deltaR = 100;
    pID = 0;

    for (unsigned int j = 0; j < tracks->at(i).found(); j++) {
      if (tracks->at(i).recHit(j)->detUnit()->type().isTrackerPixel()) {
        std::cout << "Pixel" << std::endl;
      }
      if (tracks->at(i).recHit(j)->detUnit()->type().isTrackerStrip()) {
        std::cout << "Strip" << std::endl;
      }
    }
  }

  if (minTracks_ <= tracks->size()) {
    edm::LogInfo("Demo") << "number of tracks " << tracks->size();
  }
  //  // LÖSCHEN
  //  edm::Handle<reco::GenParticleCollection> genParticles;
  //  iEvent.getByLabel("genParticles", genParticles);
  //  for (unsigned int j = 0; j < genParticles->size(); ++j) {
  //    double deltaR2 = 100;
  //    int pID2 = 0;
  //    int index = 80000;
  //    const reco::GenParticle &p = (*genParticles)[j];
  //    if ((p.status() == 1) && (p.charge() != 0)) {
  //      for (unsigned int i = 0; i < tracks->size(); i++) {
  //        if (std::sqrt(std::pow((tracks->at(i).eta() - p.eta()), 2) +
  //                      std::pow(fmod(std::abs(tracks->at(i).phi() - p.phi()),
  //                                    (4 * atan(1))),
  //                               2)) < deltaR2) {
  //          deltaR2 =
  //              std::sqrt(std::pow((tracks->at(i).eta() - p.eta()), 2) +
  //                        std::pow(fmod(std::abs(tracks->at(i).phi() -
  //                        p.phi()),
  //                                      (4 * atan(1))),
  //                                 2));
  //          index = i;
  //          pID2 = p.pdgId();
  //        }
  //      }
  //    }
  //    //    if (deltaR2 < 100)
  //    deltaR_2.push_back(deltaR2);
  //    if (deltaR2 < 0.03)
  //      particleID_2[index] = pID2;
  //    ;
  //    //    else
  //    //      particleID_2.push_back(0);
  //  }
  //  // BIS HIER
  tree_->Fill();

  //#ifdef THIS_IS_AN_EVENT_EXAMPLE
  //  Handle<ExampleData> pIn;
  //  iEvent.getByLabel("example", pIn);
  //#endif
  //
  //#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  //  ESHandle<SetupData> pSetup;
  //  iSetup.get<SetupRecord>().get(pSetup);
  //#endif
}

// ------------ method called once each job just before starting event loop
// ------------
void DemoAnalyzer::beginJob() {}

// ------------ method called once each job just after ending the event loop
// ------------
void DemoAnalyzer::endJob() {}

// ------------ method called when starting to processes a run  ------------
/*
void
DemoAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run
// ------------
/*
void
DemoAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block
// ------------
/*
void
DemoAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&,
edm::EventSetup
const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity
// block
// ------------
/*
void
DemoAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&,
edm::EventSetup
const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for
// the
// module  ------------
void DemoAnalyzer::fillDescriptions(
    edm::ConfigurationDescriptions &descriptions) {
  // The following says we do not know what parameters are allowed so do no
  // validation
  // Please change this to state exactly what you do use, even if it is no
  // parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer);
