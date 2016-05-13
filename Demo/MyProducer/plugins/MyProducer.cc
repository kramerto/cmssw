// -*- C++ -*-
//
// Package:    Demo/MyProducer
// Class:      MyProducer
//
/**\class MyProducer MyProducer.cc Demo/MyProducer/plugins/MyProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Tobias Robert Jakob Kramer
//         Created:  Mon, 09 May 2016 13:57:51 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <vector>
#include <iostream>
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

//
// class declaration
//

class MyProducer : public edm::EDProducer {
public:
  explicit MyProducer(const edm::ParameterSet &);
  ~MyProducer();

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  virtual void beginJob();
  virtual void produce(edm::Event &, const edm::EventSetup &);
  virtual void endJob();

  // virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  // virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  // virtual void beginLuminosityBlock(edm::LuminosityBlock const&,
  // edm::EventSetup const&) override;
  // virtual void endLuminosityBlock(edm::LuminosityBlock const&,
  // edm::EventSetup const&) override;

  // ----------member data ---------------------------
  edm::InputTag generalTracks_;
  typedef std::vector<reco::Track> TrackCollection;
  //  typedef reco::Track *TrackPointer;
  //  typedef std::vector<TrackPointer> TrackPointerCollection;
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
MyProducer::MyProducer(const edm::ParameterSet &iConfig) {

  // register your products
  produces<TrackCollection>();
  //  produces<TrackPointerCollection>();
  generalTracks_ = iConfig.getParameter<edm::InputTag>("generalTracks");

  /*
    // if do put with a label
    produces<ExampleData2>("label");

    // if you want to put into the Run
    produces<ExampleData2, InRun>();
    */
  // now do what ever other initialization is needed
}

MyProducer::~MyProducer() {

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void MyProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {

  edm::Handle<TrackCollection> generalTracks;
  iEvent.getByLabel(generalTracks_, generalTracks);
  std::auto_ptr<TrackCollection> MyTracks(new TrackCollection);

  for (unsigned int i = 0; i < generalTracks->size(); i++) {
    if ((std::abs(generalTracks->at(i).eta()) < 1.4) &&
        (generalTracks->at(i).found() > 9)) {
      MyTracks->push_back(generalTracks->at(i));
    }
  }
  iEvent.put(MyTracks);
  //  edm::Handle<TrackPointerCollection> generalTracks;
  //  iEvent.getByLabel(generalTracks_, generalTracks);
  //  std::auto_ptr<TrackPointerCollection> MyTracks(new
  //  TrackPointerCollection);
  //
  //  for (unsigned int i = 0; i < generalTracks->size(); i++) {
  //    if ((std::abs(generalTracks->at(i)->eta()) < 1.4) &&
  //        (generalTracks->at(i)->found() > 9)) {
  //      TrackPointer track = (generalTracks->at(i));
  //      MyTracks->push_back(track);
  //    }
  //  }
  //  iEvent.put(MyTracks);

  /* this is an EventSetup example
     //Read SetupData from the SetupRecord in the EventSetup
     ESHandle<SetupData> pSetup;
     iSetup.get<SetupRecord>().get(pSetup);
  */
}

// ------------ method called once each job just before starting event loop
// ------------
void MyProducer::beginJob() {}

// ------------ method called once each job just after ending the event loop
// ------------
void MyProducer::endJob() {}

// ------------ method called when starting to processes a run  ------------
/*
void
MyProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
MyProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block
// ------------
/*
void
MyProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup
const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block
// ------------
/*
void
MyProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup
const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the
// module  ------------
void MyProducer::fillDescriptions(
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
DEFINE_FWK_MODULE(MyProducer);
