// -*- C++ -*-
//
// Package:    Demo/DemoAnalyzer
// Class:      DemoAnalyzer
//
/**\class DemoAnalyzer DemoAnalyzer.cc
/Demo/DemoAnalyzer/plugins/DemoAnalyzer.cc

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
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/DeDxHitInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "TMatrix.h"
//#include "TBranch.h"
//#include "TMath.h"

//#include "Geometry/CommonDetUnit/interface/GeomDetType.h"

#include <unordered_map>

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
  virtual void beginRun(edm::Run const &, edm::EventSetup const &) override;

  // virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  // virtual void beginLuminosityBlock(edm::LuminosityBlock const&,
  // edm::EventSetup const&) override;
  // virtual void endLuminosityBlock(edm::LuminosityBlock const&,
  // edm::EventSetup const&) override;
  std::vector<int> crossTalkInv(const std::vector<int> &Q, const float x1,
                                const float x2, bool way, float threshold,
                                float thresholdSat);
  std::vector<int> convert(const std::vector<unsigned char> &input);
  bool clusterCleaning(const SiStripCluster *cluster, bool crosstalkInv);

  // ----------member data ---------------------------

  // Config
  std::string tracks_;
  std::string dEdxEstimator_;
  std::string genParticles_;
  bool isData_;
  std::string gainsConfig_;
  bool crossTalkInvAlgo_;
  bool correctFEDSat_;
  bool clusterCleaning_;

  // Cuts
  unsigned int nFoundHitsCut_;
  double chi2Cut_;
  double etaCut_;
  double deltaRMatching_;
  double dEdxCut_;
  double dZCut_;
  double deltaPTMatching_;
  double deltaPMatching_;

  // Calibration
  double dEdxSF[2] = {1.0, 1.21836}; // Data
  std::map<unsigned int, std::unordered_map<unsigned int, double>>
      TrackerGainsPerRuns;
  std::unordered_map<unsigned int, double> *TrackerGains;

  TTree *tree_;
  unsigned int nTracks_;
  std::vector<double> recoHarmonic2_;
  std::vector<double> charge_;
  std::vector<double> nLostHits_;
  std::vector<double> nFoundHits_;
  std::vector<double> chi2_;
  std::vector<double> ndof_;
  std::vector<double> p_;
  std::vector<int> matchedID;
  std::vector<double> deltaR_;
  std::vector<double> deltaEta_;
  std::vector<double> deltaPhi_;
  std::vector<double> dZ_;
  std::vector<double> pT_;
  std::vector<double> deltaPT_;
  std::vector<double> myHarmonic2_;
  std::vector<double> myHarmonic2Pixel_;
  std::vector<double> myHarmonic2Strip_;
  std::vector<double> myHarmonic2Amplitudes_;
  std::vector<double> myHarmonic2StripAmplitudes_;
  std::vector<double> myHarmonic2Gains_;
  std::vector<double> myHarmonic2StripGains_;
  std::vector<double> myHarmonic2Crosstalk_;
  std::vector<double> myHarmonic2StripCrosstalk_;
  std::vector<double> myHarmonic2FEDSat_;
  std::vector<double> myHarmonic2StripFEDSat_;
  std::vector<double> myHarmonic2ClusterCleaning_;
  std::vector<double> myHarmonic2StripClusterCleaning_;
  std::vector<double> dEdxPerPixel_;
  std::vector<double> dEdxPerStrip_;
  std::vector<double> eta_;
  std::vector<double> phi_;
  std::vector<unsigned int> runID_;
  unsigned int run_;
};

//
// constants, enums and typedefs
//
double infinity = std::numeric_limits<double>::infinity();
//
// static data member definitions
//

//
// constructors and destructor
//
DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet &iConfig)
    : tracks_(iConfig.getUntrackedParameter<std::string>("tracks",
                                                         "generalTracks")),
      dEdxEstimator_(iConfig.getUntrackedParameter<std::string>(
          "dEdxEstimator", "dedxHarmonic2")),
      genParticles_(iConfig.getUntrackedParameter<std::string>("genParticles",
                                                               "genParticles")),
      isData_(iConfig.getUntrackedParameter<bool>("isData", false)),
      gainsConfig_(iConfig.getUntrackedParameter<std::string>(
          "gainsConfig", "./Demo/DemoAnalyzer/data/Data13TeVGains_v2.root")),
      crossTalkInvAlgo_(
          iConfig.getUntrackedParameter<bool>("crossTalkInvAlgo", false)),
      correctFEDSat_(
          iConfig.getUntrackedParameter<bool>("correctFEDSat", false)),
      clusterCleaning_(
          iConfig.getUntrackedParameter<bool>("clusterCleaning", false)),
      nFoundHitsCut_(
          iConfig.getUntrackedParameter<unsigned int>("nFoundHitsCut", 0)),
      chi2Cut_(iConfig.getUntrackedParameter<double>("chi2Cut", infinity)),
      etaCut_(iConfig.getUntrackedParameter<double>("etaCut", infinity)),
      deltaRMatching_(iConfig.getUntrackedParameter<double>("deltaRCut", 100)),
      dEdxCut_(iConfig.getUntrackedParameter<double>("dEdxCut", -1)),
      dZCut_(iConfig.getUntrackedParameter<double>("dZCut", 9999)),
      deltaPTMatching_(
          iConfig.getUntrackedParameter<double>("deltaPTMatching", 9999)),
      deltaPMatching_(
          iConfig.getUntrackedParameter<double>("deltaPMatching", 9999)) {

  // now do what ever initialization is needed
  TrackerGains = NULL;
  run_ = 0;

  // Consumes

  edm::InputTag vertexTag("offlinePrimaryVertices");
  consumes<reco::VertexCollection>(vertexTag);
  edm::InputTag trackTag("generalTracks");
  consumes<reco::TrackCollection>(trackTag);
  edm::InputTag estimatorTag("dedxHarmonic2");
  consumes<edm::ValueMap<reco::DeDxData>>(estimatorTag);
  edm::InputTag dEdxTag("dedxHitInfo");
  consumes<reco::DeDxHitInfoAss>(dEdxTag);

  if (!isData_) {
    edm::InputTag genParticlesTag("genParticles");
    consumes<reco::GenParticleCollection>(genParticlesTag);
    dEdxSF[0] = 1.09708; // MC
    dEdxSF[1] = 1.01875; // MC
  }

  // Branches

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("Tree", "Tree");
  tree_->Branch("recoHarmonic2", &recoHarmonic2_);
  tree_->Branch("nTracks", &nTracks_);
  tree_->Branch("charge", &charge_);
  tree_->Branch("nLostHits", &nLostHits_);
  tree_->Branch("nFoundHits", &nFoundHits_);
  tree_->Branch("chi2", &chi2_);
  tree_->Branch("ndof", &ndof_);
  tree_->Branch("p", &p_);
  tree_->Branch("matchedID", &matchedID);
  tree_->Branch("deltaR", &deltaR_);
  tree_->Branch("deltaEta", &deltaEta_);
  tree_->Branch("deltaPhi", &deltaPhi_);
  tree_->Branch("dZ", &dZ_);
  tree_->Branch("pT", &pT_);
  tree_->Branch("deltaPT", &deltaPT_);
  tree_->Branch("myHarmonic2", &myHarmonic2_);
  tree_->Branch("myHarmonic2Pixel", &myHarmonic2Pixel_);
  tree_->Branch("myHarmonic2Strip", &myHarmonic2Strip_);
  tree_->Branch("dEdxPerPixel", &dEdxPerPixel_);
  tree_->Branch("dEdxPerStrip", &dEdxPerStrip_);
  tree_->Branch("eta", &eta_);
  tree_->Branch("phi", &phi_);
  tree_->Branch("runID", &runID_);

  tree_->Branch("myHarmonic2Amplitudes", &myHarmonic2Amplitudes_);
  tree_->Branch("myHarmonic2StripAmplitudes", &myHarmonic2StripAmplitudes_);
  tree_->Branch("myHarmonic2Gains", &myHarmonic2Gains_);
  tree_->Branch("myHarmonic2StripGains", &myHarmonic2StripGains_);
  tree_->Branch("myHarmonic2Crosstalk", &myHarmonic2Crosstalk_);
  tree_->Branch("myHarmonic2StripCrosstalk", &myHarmonic2StripCrosstalk_);
  tree_->Branch("myHarmonic2FEDSat", &myHarmonic2FEDSat_);
  tree_->Branch("myHarmonic2StripFEDSat", &myHarmonic2StripFEDSat_);
  tree_->Branch("myHarmonic2ClusterCleaning", &myHarmonic2ClusterCleaning_);
  tree_->Branch("myHarmonic2StripClusterCleaning",
                &myHarmonic2StripClusterCleaning_);
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

  // Get collections from event

  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByLabel("offlinePrimaryVertices", primaryVertices);

  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(tracks_, tracks);

  edm::Handle<reco::GenParticleCollection> genParticles;
  if (!isData_) {
    iEvent.getByLabel(genParticles_, genParticles);
  }

  edm::Handle<edm::ValueMap<reco::DeDxData>> dEdxTrackHandle;
  iEvent.getByLabel(dEdxEstimator_, dEdxTrackHandle);
  const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxTrackHandle.product();

  edm::Handle<reco::DeDxHitInfoAss> dedxCollH;
  iEvent.getByLabel("dedxHitInfo", dedxCollH);

  // clear vectors

  nTracks_ = 0;
  recoHarmonic2_.clear();
  nLostHits_.clear();
  nFoundHits_.clear();
  charge_.clear();
  chi2_.clear();
  ndof_.clear();
  p_.clear();
  matchedID.clear();
  deltaR_.clear();
  deltaEta_.clear();
  deltaPhi_.clear();
  dZ_.clear();
  pT_.clear();
  deltaPT_.clear();
  myHarmonic2_.clear();
  myHarmonic2Pixel_.clear();
  myHarmonic2Strip_.clear();
  dEdxPerPixel_.clear();
  dEdxPerStrip_.clear();
  eta_.clear();
  phi_.clear();
  runID_.clear();

  // loop over tracks

  for (unsigned int i = 0; i < tracks->size(); i++) {
    reco::TrackRef track = reco::TrackRef(tracks, i);
    if ((std::abs(tracks->at(i).eta()) < etaCut_) &&
        (tracks->at(i).found() >= nFoundHitsCut_) &&
        ((tracks->at(i).chi2() / tracks->at(i).ndof()) < chi2Cut_) &&
        (dEdxTrack[reco::TrackRef(tracks, i)].dEdx() > dEdxCut_) &&
        (std::abs(tracks->at(i).dz(primaryVertices->front().position())) <
         dZCut_)) {

      // Matching

      if (!isData_) {
        double deltaR = 9999;
        double deltaEta = 9999;
        double deltaPhi = 9999;
        double deltaPT = 9999;
        double matching = 9999;
        int pID = 0;
        for (unsigned int j = 0; j < genParticles->size(); ++j) {
          const reco::GenParticle &p = (*genParticles)[j];
          if ((p.status() == 1) && ((p.charge() != 0))) {

            if (std::sqrt(std::pow((tracks->at(i).eta() - p.eta()), 2) +
                          std::pow(fmod(std::abs(tracks->at(i).phi() - p.phi()),
                                        (4 * atan(1))),
                                   2)) +
                    std::abs(p.pt() - tracks->at(i).pt()) /
                        std::min(p.pt(), tracks->at(i).pt()) <
                matching) {
              deltaR = std::sqrt(
                  std::pow((tracks->at(i).eta() - p.eta()), 2) +
                  std::pow(fmod(std::abs(tracks->at(i).phi() - p.phi()),
                                (4 * atan(1))),
                           2));
              deltaEta = std::abs(tracks->at(i).eta() - p.eta());
              deltaPhi =
                  fmod(std::abs(tracks->at(i).phi() - p.phi()), (4 * atan(1)));

              pID = p.pdgId();

              deltaPT = std::abs(p.pt() - tracks->at(i).pt());
              matching = (std::abs(p.pt() - tracks->at(i).pt()) /
                          std::min(p.pt(), tracks->at(i).pt())) +
                         deltaR;
            }
          }
        }
        deltaR_.push_back(deltaR);
        deltaEta_.push_back(deltaEta);
        deltaPhi_.push_back(deltaPhi);
        deltaPT_.push_back(deltaPT);
        if (matching < 0.2) {
          matchedID.push_back(pID);
        } else {
          matchedID.push_back(0);
        }
        deltaR = 9999;
        pID = 0;
      } else {
        matchedID.push_back(0);
        deltaR_.push_back(0);
        deltaEta_.push_back(0);
        deltaPhi_.push_back(0);
        deltaPT_.push_back(0);
      }

      // Fill vectors

      nTracks_++;
      nLostHits_.push_back(tracks->at(i).lost());
      nFoundHits_.push_back(tracks->at(i).found());
      charge_.push_back(tracks->at(i).charge());
      chi2_.push_back(tracks->at(i).chi2());
      ndof_.push_back(tracks->at(i).ndof());
      dZ_.push_back(
          std::abs(tracks->at(i).dz(primaryVertices->front().position())));
      recoHarmonic2_.push_back(dEdxTrack[track].dEdx());
      p_.push_back(track->p());
      pT_.push_back(tracks->at(i).pt());
      eta_.push_back(tracks->at(i).eta());
      phi_.push_back(tracks->at(i).phi());
      runID_.push_back(run_);

      // dEdx

      const reco::DeDxHitInfo *dedxHits = NULL;
      reco::DeDxHitInfoRef dedxHitsRef = dedxCollH->get(track.key());
      if (!dedxHitsRef.isNull())
        dedxHits = &(*dedxHitsRef);
      if (dedxHits) {
        double count = 0;
        double pixelCount = 0;
        double stripCount = 0;
        double harmonic2 = 0;
        double harmonic2Pixel = 0;
        double harmonic2Strip = 0;
        for (unsigned int h = 0; h < dedxHits->size(); h++) {
          DetId detid(dedxHits->detId(h));
          double scaleFactor = dEdxSF[0];
          double Norm = (detid.subdetId() < 3) ? 3.61e-06 : 3.61e-06 * 265;
          int charge = dedxHits->charge(h);
          if (charge > 0) {
            if (detid.subdetId() >= 3) {
              if (clusterCleaning_ &&
                  !clusterCleaning(dedxHits->stripCluster(h),
                                   crossTalkInvAlgo_))
                continue;
              const SiStripCluster *cluster = dedxHits->stripCluster(h);
              std::vector<int> amplitudes = convert(cluster->amplitudes());
              if (crossTalkInvAlgo_)
                amplitudes = crossTalkInv(amplitudes, 0.10, 0.04, true, 20, 25);
              int firstStrip = cluster->firstStrip();
              int prevAPV = -1;
              double gain = 1.0;
              //              bool isSatCluster = false;
              charge = 0;
              for (unsigned int s = 0; s < amplitudes.size(); s++) {
                if (TrackerGains != NULL) {
                  int APV = (firstStrip + s) / 128;
                  if (APV != prevAPV) {
                    gain = TrackerGains->at(detid.rawId() << 3 | APV);
                    prevAPV = APV;
                  }
                }

                int StripCharge = amplitudes[s];
                if (StripCharge < 254) {
                  StripCharge = (int)(StripCharge / gain);
                  if (StripCharge >= 1023) {
                    StripCharge = 255;
                  } else if (StripCharge >= 254) {
                    StripCharge = 254;
                  }
                }

                if (StripCharge >= 254) {
                  //                  isSatCluster = true;
                }
                if ((StripCharge >= 255) && correctFEDSat_) {
                  StripCharge = 512;
                }
                charge += StripCharge;
              }
              //              if (isSatCluster) {
              //                std::cout << "SatCluster" <<std::endl;
              dEdxPerStrip_.push_back(scaleFactor * Norm * charge /
                                      dedxHits->pathlength(h));
              harmonic2Strip +=
                  pow((scaleFactor * Norm * charge / dedxHits->pathlength(h)),
                      -2.0);
              stripCount++;
            } else if (detid.subdetId() < 3) {
              scaleFactor *= dEdxSF[1];
              dEdxPerPixel_.push_back(scaleFactor * Norm * charge /
                                      dedxHits->pathlength(h));
              harmonic2Pixel +=
                  pow((scaleFactor * Norm * charge / dedxHits->pathlength(h)),
                      -2.0);
              pixelCount++;
            }
            harmonic2 += pow(
                (scaleFactor * Norm * charge / dedxHits->pathlength(h)), -2.0);
            count++;
          }
        }

        if (count != 0) {
          myHarmonic2_.push_back(pow((harmonic2 / count), -0.5));
        } else {
          myHarmonic2_.push_back(-1);
        }
        if (pixelCount != 0)
          myHarmonic2Pixel_.push_back(pow((harmonic2Pixel / pixelCount), -0.5));
        else
          myHarmonic2Pixel_.push_back(-1);
        if (stripCount != 0)
          myHarmonic2Strip_.push_back(pow((harmonic2Strip / stripCount), -0.5));
        else
          myHarmonic2Strip_.push_back(-1);
      } else {
        myHarmonic2Strip_.push_back(-1);
        myHarmonic2Pixel_.push_back(-1);
        myHarmonic2_.push_back(-1);
      }
    }
  }

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
void DemoAnalyzer::beginJob() {
  // Pixel & Strip Calibration

  // Load Calibration File
  if (isData_) {
    TrackerGainsPerRuns.clear();
    TrackerGains = NULL;

    TFile *InputFile = new TFile(gainsConfig_.c_str(), "r");
    TList *ObjList = InputFile->GetListOfKeys();
    for (int i = 0; i < ObjList->GetSize(); i++) {
      std::string Path = ObjList->At(i)->GetName();
      TObject *tmp = InputFile->Get(Path.c_str());
      if (tmp->InheritsFrom("TTree")) {
        std::string dirName = ObjList->At(i)->GetName();
        unsigned int FirstRun, LastRun;
        sscanf(dirName.c_str(), "Gains_%d_to_%d", &FirstRun, &LastRun);
        printf("Add a new gain starting at run %d\n", FirstRun);

        TTree *t1 = (TTree *)tmp;
        unsigned int tree_DetId;
        t1->SetBranchAddress("DetId", &tree_DetId);
        unsigned char tree_APVId;
        t1->SetBranchAddress("APVId", &tree_APVId);
        double tree_Gain;
        t1->SetBranchAddress("Gain", &tree_Gain);

        TrackerGains = &TrackerGainsPerRuns[FirstRun];
        for (unsigned int ientry = 0; ientry < t1->GetEntries(); ientry++) {
          t1->GetEntry(ientry);
          (*TrackerGains)[tree_DetId << 3 | tree_APVId] = tree_Gain;
        }
      }
    }
    InputFile->Close();
    delete InputFile;
  }
}

// ------------ method called once each job just after ending the event loop
// ------------
void DemoAnalyzer::endJob() {}

// ------------ method called when starting to processes a run  ------------

void DemoAnalyzer::beginRun(edm::Run const &iRun,
                            edm::EventSetup const &iSetup) {
  run_ = iRun.id().run();
  std::map<unsigned int, std::unordered_map<unsigned int, double>>::iterator it,
      itPrev = TrackerGainsPerRuns.begin();
  for (it = TrackerGainsPerRuns.begin(); it != TrackerGainsPerRuns.end();
       it++) {
    if (it->first > run_) {
      TrackerGains = &(itPrev->second);
      return;
    } // runs are ordered, so the previous iterator correspond to our run
    itPrev = it;
  }
  TrackerGains = &(itPrev->second); // just in case we go beyond the list of
                                    // run for which we have a correction
}

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

std::vector<int> DemoAnalyzer::crossTalkInv(const std::vector<int> &Q,
                                            const float x1, const float x2,
                                            bool way, float threshold,
                                            float thresholdSat) {
  const unsigned N = Q.size();
  std::vector<int> QII;
  std::vector<float> QI(N, 0);
  Double_t a = 1 - 2 * x1 - 2 * x2;
  TMatrix A(N, N);

  if (Q.size() < 2 || Q.size() > 8) {
    for (unsigned int i = 0; i < Q.size(); i++) {
      QII.push_back((int)Q[i]);
    }
    return (QII);
  }
  if (way) {
    std::vector<int>::const_iterator mQ = max_element(Q.begin(), Q.end());
    if (*mQ > 253) {
      if (*mQ == 255 && *(mQ - 1) > 253 && *(mQ + 1) > 253)
        return (Q);
      if (*(mQ - 1) > thresholdSat && *(mQ + 1) > thresholdSat &&
          *(mQ - 1) < 254 && *(mQ + 1) < 254 &&
          abs(*(mQ - 1) - *(mQ + 1)) < 40) {
        QII.push_back((10 * (*(mQ - 1)) + 10 * (*(mQ + 1))) / 2);
        return (QII);
      }
    }
  }

  for (unsigned int i = 0; i < N; i++) {
    A(i, i) = a;
    if (i < N - 1) {
      A(i + 1, i) = x1;
      A(i, i + 1) = x1;
    } else
      continue;
    if (i < N - 2) {
      A(i + 2, i) = x2;
      A(i, i + 2) = x2;
    }
  }

  if (N == 1)
    A(0, 0) = 1 / a;
  else
    A.InvertFast();

  for (unsigned int i = 0; i < N; i++) {
    for (unsigned int j = 0; j < N; j++) {
      QI[i] += A(i, j) * (float)Q[j];
    }
  }

  for (unsigned int i = 0; i < QI.size(); i++) {
    if (QI[i] < threshold)
      QI[i] = 0;
    QII.push_back((int)QI[i]);
  }

  return (QII);
}

std::vector<int>
DemoAnalyzer::convert(const std::vector<unsigned char> &input) {
  std::vector<int> output;
  for (unsigned int i = 0; i < input.size(); i++) {
    output.push_back((int)input[i]);
  }
  return (output);
}

bool DemoAnalyzer::clusterCleaning(const SiStripCluster *cluster,
                                   bool crosstalkInv) {
  if (!cluster)
    return (true);
  std::vector<int> ampls = convert(cluster->amplitudes());
  if (crosstalkInv)
    ampls = crossTalkInv(ampls, 0.10, 0.04, true, 20, 25);

  Int_t NofMax = 0;
  Int_t recur255 = 1;
  Int_t recur254 = 1;
  bool MaxOnStart = false;
  bool MaxInMiddle = false, MaxOnEnd = false;
  Int_t MaxPos = 0;

  if (ampls.size() != 1 &&
      ((ampls[0] > ampls[1]) ||
       (ampls.size() > 2 && ampls[0] == ampls[1] && ampls[1] > ampls[2] &&
        ampls[0] != 254 && ampls[0] != 255) ||
       (ampls.size() == 2 && ampls[0] == ampls[1] && ampls[0] != 254 &&
        ampls[0] != 255))) {
    NofMax = NofMax + 1;
    MaxOnStart = true;
  }

  if (ampls.size() > 2) {
    for (unsigned int i = 1; i < ampls.size() - 1; i++) {
      if ((ampls[i] > ampls[i - 1] && ampls[i] > ampls[i + 1]) ||
          (ampls.size() > 3 && i > 0 && i < ampls.size() - 2 &&
           ampls[i] == ampls[i + 1] && ampls[i] > ampls[i - 1] &&
           ampls[i] > ampls[i + 2] && ampls[i] != 254 && ampls[i] != 255)) {
        NofMax = NofMax + 1;
        MaxInMiddle = true;
        MaxPos = i;
      }
      if (ampls[i] == 255 && ampls[i] == ampls[i - 1]) {
        recur255 = recur255 + 1;
        MaxPos = i - (recur255 / 2);
        if (ampls[i] > ampls[i + 1]) {
          NofMax = NofMax + 1;
          MaxInMiddle = true;
        }
      }
      if (ampls[i] == 254 && ampls[i] == ampls[i - 1]) {
        recur254 = recur254 + 1;
        MaxPos = i - (recur254 / 2);
        if (ampls[i] > ampls[i + 1]) {
          NofMax = NofMax + 1;
          MaxInMiddle = true;
        }
      }
    }
  }

  if (ampls.size() > 1) {
    if (ampls[ampls.size() - 1] > ampls[ampls.size() - 2] ||
        (ampls.size() > 2 &&
         ampls[ampls.size() - 1] == ampls[ampls.size() - 2] &&
         ampls[ampls.size() - 2] > ampls[ampls.size() - 3]) ||
        ampls[ampls.size() - 1] == 255) {
      NofMax = NofMax + 1;
      MaxOnEnd = true;
    }
  }

  if (ampls.size() == 1) {
    NofMax = 1;
  }

  //               ____
  //              |    |____
  //          ____|    |    |
  //         |    |    |    |____
  //     ____|    |    |    |    |
  //    |    |    |    |    |    |____
  //  __|____|____|____|____|____|____|__
  //    C_Mnn C_Mn C_M  C_D  C_Dn C_Dnn
  //
  bool shapecdtn = false;

  if (crosstalkInv) {
    if (NofMax == 1) {
      shapecdtn = true;
    }
    return (shapecdtn);
  }

  Float_t C_M = 0.0;
  Float_t C_D = 0.0;
  Float_t C_Mn = 10000;
  Float_t C_Dn = 10000;
  Float_t C_Mnn = 10000;
  Float_t C_Dnn = 10000;
  Int_t CDPos;
  Float_t coeff1 = 1.7;
  Float_t coeff2 = 2.0;
  Float_t coeffn = 0.10;
  Float_t coeffnn = 0.02;
  Float_t noise = 4.0;

  if (NofMax == 1) {

    if (MaxOnStart == true) {
      C_M = (Float_t)ampls[0];
      C_D = (Float_t)ampls[1];
      if (ampls.size() < 3)
        shapecdtn = true;
      else if (ampls.size() == 3) {
        C_Dn = (Float_t)ampls[2];
        if (C_Dn <=
                coeff1 * coeffn * C_D + coeff2 * coeffnn * C_M + 2 * noise ||
            C_D == 255)
          shapecdtn = true;
      } else if (ampls.size() > 3) {
        C_Dn = (Float_t)ampls[2];
        C_Dnn = (Float_t)ampls[3];
        if ((C_Dn <=
                 coeff1 * coeffn * C_D + coeff2 * coeffnn * C_M + 2 * noise ||
             C_D == 255) &&
            C_Dnn <=
                coeff1 * coeffn * C_Dn + coeff2 * coeffnn * C_D + 2 * noise) {
          shapecdtn = true;
        }
      }
    }

    if (MaxOnEnd == true) {
      C_M = (Float_t)ampls[ampls.size() - 1];
      C_D = (Float_t)ampls[ampls.size() - 2];
      if (ampls.size() < 3)
        shapecdtn = true;
      else if (ampls.size() == 3) {
        C_Dn = (Float_t)ampls[0];
        if (C_Dn <=
                coeff1 * coeffn * C_D + coeff2 * coeffnn * C_M + 2 * noise ||
            C_D == 255)
          shapecdtn = true;

      } else if (ampls.size() > 3) {
        C_Dn = (Float_t)ampls[ampls.size() - 3];
        C_Dnn = (Float_t)ampls[ampls.size() - 4];
        if ((C_Dn <=
                 coeff1 * coeffn * C_D + coeff2 * coeffnn * C_M + 2 * noise ||
             C_D == 255) &&
            C_Dnn <=
                coeff1 * coeffn * C_Dn + coeff2 * coeffnn * C_D + 2 * noise) {
          shapecdtn = true;
        }
      }
    }

    if (MaxInMiddle == true) {
      C_M = (Float_t)ampls[MaxPos];
      int LeftOfMaxPos = MaxPos - 1;
      if (LeftOfMaxPos <= 0)
        LeftOfMaxPos = 0;
      int RightOfMaxPos = MaxPos + 1;
      if (RightOfMaxPos >= (int)ampls.size())
        RightOfMaxPos = ampls.size() - 1;

      if (ampls[LeftOfMaxPos] < ampls[RightOfMaxPos]) {
        C_D = (Float_t)ampls[RightOfMaxPos];
        C_Mn = (Float_t)ampls[LeftOfMaxPos];
        CDPos = RightOfMaxPos;
      } else {
        C_D = (Float_t)ampls[LeftOfMaxPos];
        C_Mn = (Float_t)ampls[RightOfMaxPos];
        CDPos = LeftOfMaxPos;
      }
      if (C_Mn < coeff1 * coeffn * C_M + coeff2 * coeffnn * C_D + 2 * noise ||
          C_M == 255) {
        if (ampls.size() == 3)
          shapecdtn = true;
        else if (ampls.size() > 3) {
          if (CDPos > MaxPos) {
            if (ampls.size() - CDPos - 1 == 0) {
              C_Dn = 0;
              C_Dnn = 0;
            }
            if (ampls.size() - CDPos - 1 == 1) {
              C_Dn = (Float_t)ampls[CDPos + 1];
              C_Dnn = 0;
            }
            if (ampls.size() - CDPos - 1 > 1) {
              C_Dn = (Float_t)ampls[CDPos + 1];
              C_Dnn = (Float_t)ampls[CDPos + 2];
            }
            if (MaxPos >= 2) {
              C_Mnn = (Float_t)ampls[MaxPos - 2];
            } else if (MaxPos < 2)
              C_Mnn = 0;
          }
          if (CDPos < MaxPos) {
            if (CDPos == 0) {
              C_Dn = 0;
              C_Dnn = 0;
            }
            if (CDPos == 1) {
              C_Dn = (Float_t)ampls[0];
              C_Dnn = 0;
            }
            if (CDPos > 1) {
              C_Dn = (Float_t)ampls[CDPos - 1];
              C_Dnn = (Float_t)ampls[CDPos - 2];
            }
            if (ampls.size() - LeftOfMaxPos > 1 &&
                MaxPos + 2 < (int)(ampls.size()) - 1) {
              C_Mnn = (Float_t)ampls[MaxPos + 2];
            } else
              C_Mnn = 0;
          }
          if ((C_Dn <=
                   coeff1 * coeffn * C_D + coeff2 * coeffnn * C_M + 2 * noise ||
               C_D == 255) &&
              C_Mnn <=
                  coeff1 * coeffn * C_Mn + coeff2 * coeffnn * C_M + 2 * noise &&
              C_Dnn <=
                  coeff1 * coeffn * C_Dn + coeff2 * coeffnn * C_D + 2 * noise) {
            shapecdtn = true;
          }
        }
      }
    }
  }
  if (ampls.size() == 1) {
    shapecdtn = true;
  }
  return (shapecdtn);
}

// define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer);
