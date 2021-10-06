// -*- C++ -*-
//
// Package:    GenParticleAnalyzer
// Class:      GenParticleAnalyzer
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

//classes to extract GenParticle information
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

//classes to save data
#include "TTree.h"
#include "TFile.h"
#include<vector>

//
// class declaration
//

class GenParticleAnalyzer : public edm::EDAnalyzer {
   public:
      explicit GenParticleAnalyzer(const edm::ParameterSet&);
      ~GenParticleAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

     std::vector<std::string>  particle;

      // ----------member data ---------------------------

      TTree *mtree;

      int numGenPart;
      std::vector<int> GenPart_status;
      std::vector<float> GenPart_pt;
      std::vector<float> GenPart_eta;
      std::vector<float> GenPart_mass;
      std::vector<int> GenPart_pdgId;
      std::vector<float> GenPart_phi;
      std::vector<float> GenPart_px;
      std::vector<float> GenPart_py;
      std::vector<float> GenPart_pz;
      std::vector<int> GenPart_mompdgId;
//Daugther store
      int numGenDau;
      std::vector<int> GenDau_status;
      std::vector<float> GenDau_pt;
      std::vector<float> GenDau_eta;
      std::vector<float> GenDau_mass;
      std::vector<int> GenDau_pdgId;
      std::vector<float> GenDau_phi;
      std::vector<float> GenDau_px;
      std::vector<float> GenDau_py;
      std::vector<float> GenDau_pz;
      std::vector<int> GenDau_mompdgId;

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

GenParticleAnalyzer::GenParticleAnalyzer(const edm::ParameterSet& iConfig):
particle(iConfig.getParameter<std::vector<std::string> >("input_particle"))

{
//now do what ever initialization is needed
	edm::Service<TFileService> fs;
	mtree = fs->make<TTree>("Events", "Events");

    mtree->Branch("numGenPart",&numGenPart);
    mtree->GetBranch("numGenPart")->SetTitle("number of generator particles");
    mtree->Branch("GenPart_pt",&GenPart_pt);
    mtree->GetBranch("GenPart_pt")->SetTitle("generator particle transverse momentum");
    mtree->Branch("GenPart_eta",&GenPart_eta);
    mtree->GetBranch("GenPart_eta")->SetTitle("generator particle pseudorapidity");
    mtree->Branch("GenPart_mass",&GenPart_mass);
    mtree->GetBranch("GenPart_mass")->SetTitle("generator particle mass");
    mtree->Branch("GenPart_pdgId",&GenPart_pdgId);
    mtree->GetBranch("GenPart_pdgId")->SetTitle("generator particle PDG id");
    mtree->Branch("GenPart_mompdgId",&GenPart_mompdgId);
    mtree->GetBranch("GenPart_mompdgId")->SetTitle("generator particle mother PDG id");
    mtree->Branch("GenPart_phi",&GenPart_phi);
    mtree->GetBranch("GenPart_phi")->SetTitle("generator particle azimuthal angle of momentum vector");
    mtree->Branch("GenPart_px",&GenPart_px);
    mtree->GetBranch("GenPart_px")->SetTitle("generator particle x coordinate of momentum vector");
    mtree->Branch("GenPart_py",&GenPart_py);
    mtree->GetBranch("GenPart_py")->SetTitle("generator particle y coordinate of momentum vector");
    mtree->Branch("GenPart_pz",&GenPart_pz);
    mtree->GetBranch("GenPart_pz")->SetTitle("generator particle z coordinate of momentum vector");
    mtree->Branch("GenPart_status",&GenPart_status);
    mtree->GetBranch("GenPart_status")->SetTitle("Particle status. 1=stable");

    mtree->Branch("numGenDau",&numGenDau);
    mtree->GetBranch("numGenDau")->SetTitle("number of generator particles");
    mtree->Branch("GenDau_pt",&GenDau_pt);
    mtree->GetBranch("GenDau_pt")->SetTitle("generator particle transverse momentum");
    mtree->Branch("GenDau_eta",&GenDau_eta);
    mtree->GetBranch("GenDau_eta")->SetTitle("generator particle pseudorapidity");
    mtree->Branch("GenDau_mass",&GenDau_mass);
    mtree->GetBranch("GenDau_mass")->SetTitle("generator particle mass");
    mtree->Branch("GenDau_pdgId",&GenDau_pdgId);
    mtree->GetBranch("GenDau_pdgId")->SetTitle("generator particle PDG id");
    mtree->Branch("GenDau_mompdgId",&GenDau_mompdgId);
    mtree->GetBranch("GenDau_mompdgId")->SetTitle("generator particle mother PDG id");
    mtree->Branch("GenDau_phi",&GenDau_phi);
    mtree->GetBranch("GenDau_phi")->SetTitle("generator particle azimuthal angle of momentum vector");
    mtree->Branch("GenDau_px",&GenDau_px);
    mtree->GetBranch("GenDau_px")->SetTitle("generator particle x coordinate of momentum vector");
    mtree->Branch("GenDau_py",&GenDau_py);
    mtree->GetBranch("GenDau_py")->SetTitle("generator particle y coordinate of momentum vector");
    mtree->Branch("GenDau_pz",&GenDau_pz);
    mtree->GetBranch("GenDau_pz")->SetTitle("generator particle z coordinate of momentum vector");
    mtree->Branch("GenDau_status",&GenDau_status);
    mtree->GetBranch("GenDau_status")->SetTitle("Particle status. 1=stable");

}

GenParticleAnalyzer::~GenParticleAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
GenParticleAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   numGenPart=0;
   GenPart_pt.clear();
   GenPart_eta.clear();
   GenPart_mass.clear();
   GenPart_pdgId.clear();
   GenPart_mompdgId.clear();
   GenPart_phi.clear();
   GenPart_px.clear();
   GenPart_py.clear();
   GenPart_pz.clear();
   GenPart_status.clear();

   numGenDau=0;
   GenDau_pt.clear();
   GenDau_eta.clear();
   GenDau_mass.clear();
   GenDau_pdgId.clear();
   GenDau_mompdgId.clear();
   GenDau_phi.clear();
   GenDau_px.clear();
   GenDau_py.clear();
   GenDau_pz.clear();
   GenDau_status.clear();

   Handle<reco::GenParticleCollection> gens;
   iEvent.getByLabel("genParticles", gens);

   unsigned int i;
   string s1,s2;
   std::vector<int> status_parsed;
   std::vector<int> pdgId_parsed;
   std::string delimiter = ":";

   for(i=0;i<particle.size();i++)
   {
       //get status and pgdId from configuration
       s1=particle[i].substr(0,particle[i].find(delimiter));
       s2=particle[i].substr(particle[i].find(delimiter)+1,particle[i].size());
       //parse string to int
       status_parsed.push_back(stoi(s1));
       pdgId_parsed.push_back(stoi(s2));
   }

  if(gens.isValid())
  {
        for (reco::GenParticleCollection::const_iterator itGenPart=gens->begin(); itGenPart!=gens->end(); ++itGenPart)
        {
               //loop trough all particles selected in configuration
               for(i=0;i<particle.size();i++)
               {
                  if((status_parsed[i]==itGenPart->status() && pdgId_parsed[i]==itGenPart->pdgId())||(status_parsed[i]==0 && pdgId_parsed[i]==0)||(status_parsed[i]==0 && pdgId_parsed[i]==itGenPart->pdgId()))
                {
		  numGenPart++;//arreglado conteo de particulas generadas
                  GenPart_pt.push_back(itGenPart->pt());
                  GenPart_eta.push_back(itGenPart->eta());
                  GenPart_mass.push_back(itGenPart->mass());
                  GenPart_pdgId.push_back(itGenPart->pdgId());
                  GenPart_phi.push_back(itGenPart->phi());
                  GenPart_status.push_back(itGenPart->status());
                  GenPart_px.push_back(itGenPart->px());
                  GenPart_py.push_back(itGenPart->py());
                  GenPart_pz.push_back(itGenPart->pz());
		  if(abs(itGenPart->pdgId())!=556)GenPart_mompdgId.push_back(itGenPart->mother()->pdgId());
      else GenPart_mompdgId.push_back(0);
		  //Daugther store
		  int n = itGenPart->numberOfDaughters();
		  numGenDau=n;
      		  for(int j = 0; j < n; ++ j) {
       			const reco::Candidate * d = itGenPart->daughter( j );
			GenDau_pt.push_back(d->pt());
                	GenDau_eta.push_back(d->eta());
                  	GenDau_mass.push_back(d->mass());
                  	GenDau_pdgId.push_back(d->pdgId());
                  	GenDau_phi.push_back(d->phi());
                  	GenDau_status.push_back(d->status());
                  	GenDau_px.push_back(d->px());
                  	GenDau_py.push_back(d->py());
                  	GenDau_pz.push_back(d->pz());
                  	GenDau_mompdgId.push_back(itGenPart->pdgId());
			}
                }
               }
        }
  }

  mtree->Fill();
  return;

}

// ------------ method called once each job just before starting event loop  ------------
void
GenParticleAnalyzer::beginJob()
{}

// ------------ method called once each job just after ending the event loop  ------------
void
GenParticleAnalyzer::endJob()
{}

// ------------ method called when starting to processes a run  ------------
void
GenParticleAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a run  ------------
void
GenParticleAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{}
// ------------ method called when starting to processes a luminosity block  ------------
void
GenParticleAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a luminosity block  ------------
void
GenParticleAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenParticleAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
 //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenParticleAnalyzer);
