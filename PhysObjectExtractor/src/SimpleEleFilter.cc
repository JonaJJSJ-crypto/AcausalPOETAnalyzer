// -*- C++ -*-
//
// Package:    SimpleFilter
// Class:      SimpleFilter
//
/**\class SimpleEleFilter SimpleEleFilter.cc PhysObjectExtractorTool/SimpleEleFilter/src/SimpleEleFilter.cc

 Description: [one line class summary]

This is a simple filter example to filter on at least a muon and at least
a tau of certain characteristics, which are mostly configurable.

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:
//         Created:  Sat Jul 17 22:23:23 CEST 2021
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//Classes to extract ctron information
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

/*//classes to extract Muon information
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//classes to extract tau information
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"*/

//
// class declaration
//

class SimpleEleFilter : public edm::EDFilter {
   public:
      explicit SimpleEleFilter(const edm::ParameterSet&);
      ~SimpleEleFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
  edm::InputTag electronInput;
  double ele_minpt_;
  //double ele_etacut_;


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
SimpleEleFilter::SimpleEleFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  electronInput = iConfig.getParameter<edm::InputTag>("InputCollectionElectrons");
  ele_minpt_ = iConfig.getParameter<double>("ele_minpt");
  //ele_etacut_ = iConfig.getParameter<double>("ele_etacut");
  //ele_mindxy_ = iConfig.getParameter<double>("ele_mindxy")
}


SimpleEleFilter::~SimpleEleFilter()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
SimpleEleFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;
  using namespace std;


 //Filter on at least one good muon
  Handle<reco::GsfElectronCollection> myelectrons;
  iEvent.getByLabel(electronInput, myelectrons);

  Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(InputTag("offlinePrimaryVertices"), vertices);

  bool isGoodElectron =false;

  if(myelectrons.isValid()){
    //math::XYZPoint pv(vertices->begin()->position());
    for (reco::GsfElectronCollection::const_iterator itElec=myelectrons->begin(); itElec!=myelectrons->end(); ++itElec){
      auto trk = itElec->gsfTrack();
      if(itElec->pt()>ele_minpt_){
	isGoodElectron = true;
      }
    }
  }

   return isGoodElectron;
}

// ------------ method called once each job just before starting event loop  ------------
void
SimpleEleFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
SimpleEleFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool
SimpleEleFilter::beginRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool
SimpleEleFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool
SimpleEleFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool
SimpleEleFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SimpleEleFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(SimpleEleFilter);
