//
// Package:    ElectronAnalyzer
// Class:      ElectronAnalyzer
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

//classes to extract electron information
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//Secondary Vertex Clases
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

//Drawer
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"

//jetCollectio
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

//classes to save data
#include "TTree.h"
#include "TFile.h"
#include<vector>

//
// class declaration
//

class ElectronAnalyzer : public edm::EDAnalyzer {
public:
  explicit ElectronAnalyzer(const edm::ParameterSet&);
  ~ElectronAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  virtual float effectiveArea0p3cone(float eta);

  //declare the input tag for GsfElectronCollection
  edm::InputTag electronInput;

  // ----------member data ---------------------------
  bool isData;

  TTree *mtree;
  int numelectron; //number of reco electrons in the event
  int numgenelec; //# of gen electrons in the event
  int numsecvec;  //number of secondary vertex disp in the event
  int numdisp;  //tool for finding best match for electrons
  std::vector<float> electron_e;
  std::vector<float> electron_pt;
  std::vector<float> electron_px;
  std::vector<float> electron_py;
  std::vector<float> electron_pz;
  std::vector<float> electron_eta;
  std::vector<float> electron_phi;
  std::vector<float> electron_ch;
  std::vector<float> genelec_e;
  std::vector<float> genelec_pt;
  std::vector<float> genelec_px;
  std::vector<float> genelec_py;
  std::vector<float> genelec_pz;
  std::vector<float> genelec_eta;
  std::vector<float> genelec_phi;
  std::vector<float> genelec_ch;
  std::vector<float> electron_iso;
  std::vector<bool> electron_isLoose;
  std::vector<bool> electron_isMedium;
  std::vector<bool> electron_isTight;
  std::vector<float> electron_dxy;
  std::vector<float> electron_dz;
  std::vector<float> electron_dxyError;
  std::vector<float> electron_dzError;
  std::vector<float> secvec_posx;
  std::vector<float> secvec_posy;
  std::vector<float> secvec_posz;
  std::vector<float> secvec_poserrorx;
  std::vector<float> secvec_poserrory;
  std::vector<float> secvec_poserrorz;
  std::vector<float> secvec_disp;
  std::vector<float> secvec_dispR;
  std::vector<float> secvec_deltaR;
  std::vector<float> secvec_deltaR1;
  std::vector<float> secvec_deltaR2;
  std::vector<int> electron_Bsecvec;//best match for electron sec vec
  std::vector<int> electron_secN;//# of printed electrons in Secvert
  std::vector<float> electron_deltaRsim;
  std::vector<float> electron_deltaR1sim;
  std::vector<float> electron_deltaR2sim;
  std::vector<float> electron_deltaRtrue;
  std::vector<float> electron_deltaR1true;
  std::vector<float> electron_deltaR2true;
  std::vector<float> Zjet_pt;
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

ElectronAnalyzer::ElectronAnalyzer(const edm::ParameterSet& iConfig)
{
//now do what ever initialization is needed
  isData = iConfig.getParameter<bool>("isData");
  electronInput = iConfig.getParameter<edm::InputTag>("InputCollection");
  edm::Service<TFileService> fs;
  mtree = fs->make<TTree>("Events", "Events");


  mtree->Branch("numberelectron",&numelectron);
  mtree->GetBranch("numberelectron")->SetTitle("number of reco electrons");
  mtree->Branch("numbergenelec",&numgenelec);
  mtree->GetBranch("numbergenelec")->SetTitle("number of gen electrons");
  mtree->Branch("numbersecvec",&numsecvec);
  mtree->GetBranch("numbersecvec")->SetTitle("number of secondary vertex displascements");
  mtree->Branch("electron_e",&electron_e);
  mtree->GetBranch("electron_e")->SetTitle("electron energy");
  mtree->Branch("electron_pt",&electron_pt);
  mtree->GetBranch("electron_pt")->SetTitle("electron transverse momentum");
  mtree->Branch("electron_px",&electron_px);
  mtree->GetBranch("electron_px")->SetTitle("electron momentum x-component");
  mtree->Branch("electron_py",&electron_py);
  mtree->GetBranch("electron_py")->SetTitle("electron momentum y-component");
  mtree->Branch("electron_pz",&electron_pz);
  mtree->GetBranch("electron_pz")->SetTitle("electron momentum z-component");
  mtree->Branch("electron_eta",&electron_eta);
  mtree->GetBranch("electron_eta")->SetTitle("electron pseudorapidity");
  mtree->Branch("electron_phi",&electron_phi);
  mtree->GetBranch("electron_phi")->SetTitle("electron polar angle");
  mtree->Branch("electron_ch",&electron_ch);
  mtree->GetBranch("electron_ch")->SetTitle("electron charge");
  mtree->Branch("genelec_e",&genelec_e);
  mtree->GetBranch("genelec_e")->SetTitle("electron gen energy");
  mtree->Branch("genelec_pt",&genelec_pt);
  mtree->GetBranch("genelec_pt")->SetTitle("electron gen transverse momentum");
  mtree->Branch("genelec_px",&genelec_px);
  mtree->GetBranch("genelec_px")->SetTitle("electron gen momentum x-component");
  mtree->Branch("genelec_py",&genelec_py);
  mtree->GetBranch("genelec_py")->SetTitle("electron gen momentum y-component");
  mtree->Branch("genelec_pz",&genelec_pz);
  mtree->GetBranch("genelec_pz")->SetTitle("electron gen momentum z-component");
  mtree->Branch("genelec_eta",&genelec_eta);
  mtree->GetBranch("genelec_eta")->SetTitle("electron gen pseudorapidity");
  mtree->Branch("genelec_phi",&genelec_phi);
  mtree->GetBranch("genelec_phi")->SetTitle("electron gen polar angle");
  mtree->Branch("genelec_ch",&genelec_ch);
  mtree->GetBranch("genelec_ch")->SetTitle("electron gen charge");
  mtree->Branch("electron_iso",&electron_iso);
  mtree->GetBranch("electron_iso")->SetTitle("electron isolation");
  mtree->Branch("electron_isLoose",&electron_isLoose);
  mtree->GetBranch("electron_isLoose")->SetTitle("electron tagged loose");
  mtree->Branch("electron_isMedium",&electron_isMedium);
  mtree->GetBranch("electron_isMedium")->SetTitle("electron tagged medium");
  mtree->Branch("electron_isTight",&electron_isTight);
  mtree->GetBranch("electron_isTight")->SetTitle("electron tagged tight");
  mtree->Branch("electron_dxy",&electron_dxy);
  mtree->GetBranch("electron_dxy")->SetTitle("electron transverse plane impact parameter (mm)");
  mtree->Branch("electron_dz",&electron_dz);
  mtree->GetBranch("electron_dz")->SetTitle("electron longitudinal impact parameter (mm)");
  mtree->Branch("electron_dxyError",&electron_dxyError);
  mtree->GetBranch("electron_dxyError")->SetTitle("electron transverse impact parameter uncertainty (mm)");
  mtree->Branch("electron_dzError",&electron_dzError);
  mtree->GetBranch("electron_dzError")->SetTitle("electron longitudinal impact parameter uncertainty (mm)");
  mtree->Branch("secvec_posx",&secvec_posx);
  mtree->GetBranch("secvec_posx")->SetTitle("secvec position x (mm)");
  mtree->Branch("secvec_posy",&secvec_posy);
  mtree->GetBranch("secvec_posy")->SetTitle("secvec position y (mm)");
  mtree->Branch("secvec_posz",&secvec_posz);
  mtree->GetBranch("secvec_posz")->SetTitle("secvec position z (mm)");
  mtree->Branch("secvec_poserrorx",&secvec_poserrorx);
  mtree->GetBranch("secvec_poserrorx")->SetTitle("secvec position x error (mm)");
  mtree->Branch("secvec_poserrory",&secvec_poserrory);
  mtree->GetBranch("secvec_poserrory")->SetTitle("secvec position y error (mm)");
  mtree->Branch("secvec_poserrorz",&secvec_poserrorz);
  mtree->GetBranch("secvec_poserrorz")->SetTitle("secvec position z error (mm)");
  mtree->Branch("secvec_disp",&secvec_disp);
  mtree->GetBranch("secvec_disp")->SetTitle("secvec displacement from primary vertex (mm)");
  mtree->Branch("secvec_dispR",&secvec_dispR);
  mtree->GetBranch("secvec_dispR")->SetTitle("secvec weigthed displacement from primary vertex");
  mtree->Branch("secvec_deltaR",&secvec_deltaR);
  mtree->GetBranch("secvec_deltaR")->SetTitle("1-2tracks deltaR");
  mtree->Branch("secvec_deltaR1",&secvec_deltaR1);
  mtree->GetBranch("secvec_deltaR1")->SetTitle("1st track deltaR");
  mtree->Branch("secvec_deltaR2",&secvec_deltaR2);
  mtree->GetBranch("secvec_deltaR2")->SetTitle("2nd track deltaR");
  mtree->Branch("electron_Bsecvec",&electron_Bsecvec);
  mtree->GetBranch("electron_Bsecvec")->SetTitle("best match electrons for secondary vertex");
  mtree->Branch("electron_secN",&electron_secN);
  mtree->GetBranch("electron_secN")->SetTitle("# electrons for secondary vertex");
  mtree->Branch("electron_deltaRsim",&electron_deltaRsim);
  mtree->GetBranch("electron_deltaRsim")->SetTitle("1-2Electron deltaRsim");
  mtree->Branch("electron_deltaR1sim",&electron_deltaR1sim);
  mtree->GetBranch("electron_deltaR1sim")->SetTitle("1 Electron deltaRsim");
  mtree->Branch("electron_deltaR2sim",&electron_deltaR2sim);
  mtree->GetBranch("electron_deltaR2sim")->SetTitle("2 Electron deltaRsim");
  mtree->Branch("electron_deltaRtrue",&electron_deltaRtrue);
  mtree->GetBranch("electron_deltaRtrue")->SetTitle("1-2Electron deltaRtrue");
  mtree->Branch("electron_deltaR1true",&electron_deltaR1true);
  mtree->GetBranch("electron_deltaR1true")->SetTitle("1 Electron deltaRtrue");
  mtree->Branch("electron_deltaR2true",&electron_deltaR2true);
  mtree->GetBranch("electron_deltaR2true")->SetTitle("2 Electron deltaRtrue");
  mtree->Branch("Zjet_pt",&Zjet_pt);
  mtree->GetBranch("Zjet_pt")->SetTitle("Z daugthers pt");

}

ElectronAnalyzer::~ElectronAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

float
ElectronAnalyzer::effectiveArea0p3cone(float eta)
{
  if(fabs(eta) < 1.0) return 0.13;
  else if(fabs(eta) < 1.479) return 0.14;
  else if(fabs(eta) < 2.0) return 0.07;
  else if(fabs(eta) < 2.2) return 0.09;
  else if(fabs(eta) < 2.3) return 0.11;
  else if(fabs(eta) < 2.4) return 0.11;
  else return 0.14;
}

// ------------ method called for each event  ------------
void
ElectronAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;

  //cout<<"\n\n\n==========================================ESTA=============================================\n\n\n"<<endl;

  Handle<reco::GsfElectronCollection> myelectrons;
  iEvent.getByLabel(electronInput, myelectrons);
  Handle<reco::ConversionCollection> hConversions;
  iEvent.getByLabel("allConversions", hConversions);
  Handle<reco::BeamSpot> bsHandle;
  iEvent.getByLabel("offlineBeamSpot", bsHandle);
  const reco::BeamSpot &beamspot = *bsHandle.product();
  Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(InputTag("offlinePrimaryVertices"), vertices);
  math::XYZPoint pv(vertices->begin()->position());
  Handle<double> rhoHandle;
  iEvent.getByLabel(InputTag("fixedGridRhoAll"), rhoHandle);

  numelectron = 0;
  numgenelec = 0;
  numsecvec = 0;
  numdisp = 0;
  electron_e.clear();
  electron_pt.clear();
  electron_px.clear();
  electron_py.clear();
  electron_pz.clear();
  electron_eta.clear();
  electron_phi.clear();
  electron_ch.clear();
  genelec_e.clear();
  genelec_pt.clear();
  genelec_px.clear();
  genelec_py.clear();
  genelec_pz.clear();
  genelec_eta.clear();
  genelec_phi.clear();
  genelec_ch.clear();
  electron_iso.clear();
  electron_isLoose.clear();
  electron_isMedium.clear();
  electron_isTight.clear();
  electron_dxy.clear();
  electron_dz.clear();
  electron_dxyError.clear();
  electron_dzError.clear();
  secvec_posx.clear();
  secvec_posy.clear();
  secvec_posz.clear();
  secvec_poserrorx.clear();
  secvec_poserrory.clear();
  secvec_poserrorz.clear();
  secvec_disp.clear();
  secvec_dispR.clear();
  secvec_deltaR.clear();
  secvec_deltaR1.clear();
  secvec_deltaR2.clear();
  electron_secN.clear();
  electron_deltaRsim.clear();
  electron_deltaR1sim.clear();
  electron_deltaR2sim.clear();
  electron_deltaRtrue.clear();
  electron_deltaR1true.clear();
  electron_deltaR2true.clear();
  Zjet_pt.clear();

//drawer
bool Draw=true;//por terminar
if(!Draw){
   Handle<GenParticleCollection> genParticles;
   iEvent.getByLabel("genParticles", genParticles);
   const GenParticle & t = (*genParticles)[0];
   auto ptmp1=t.p4(), ptmp2=t.p4(), ptmpz1=t.p4(), ptmpz2=t.p4();//almacen de hijas
   int counter=0;//cuenta las hijas almacenadas
   for(size_t i = 0; i < genParticles->size(); ++ i) {

     const GenParticle & p = (*genParticles)[i];
     int id = p.pdgId();
     int st = p.status();
     const Candidate * mom = p.mother();
     if(id==556 || id==-556 || id==23){
     cout<<"\n\n\nPDG: "<<id<<" Status: "<< st<<" mother "<< mom->pdgId() <<endl;
     double charge = p.charge();
     double pt = p.pt(), eta = p.eta(), phi = p.phi(), mass = p.mass();
     double vx = p.vx(), vy = p.vy(), vz = p.vz();
     cout <<"pt: "<<pt<<" eta: "<<eta<<" phi: "<<phi<<" mass: "<<mass<<" vx: "<<vx<<" vy: "<<vy<<" vz: "<<vz<<" c: "<<charge<<endl;
     int n = p.numberOfDaughters();
      for(int j = 0; j < n; ++ j) {
       const Candidate * d = p.daughter( j );
       int dauId = d->pdgId();
       cout<<"  daugther: "<<dauId<<" Dpt: "<<d->pt()<< endl;
       if(id==23){
	if(mom->pdgId()==556 || mom->pdgId()==-556)
	{float Zdpt=d->pt();
	Zjet_pt.push_back(Zdpt);}
       }
       //auto ptmp1=p.p4(), ptmp2=p.p4(), ptmpz1=p.p4(), ptmpz2=p.p4();

       if(dauId==11){
	 counter++;
         auto p1=d->p4();
         ptmp1=d->p4();
        //float phi1=d->phi(), eta1=d->eta();
        cout << "   p4_1: "<<p1<<" pt: "<<d->pt()<< endl;
        for(size_t k = 0; k < genParticles->size(); ++ k) {
          const GenParticle & q = (*genParticles)[k];
          auto q4=q.p4();
          float DeltaRMC1=deltaR(p1,q4);
          //cout<<"deltaRsim: "<<DeltaRMC1<<endl;
          if(DeltaRMC1>10){DeltaRMC1=9;}
	  electron_deltaRsim.push_back(DeltaRMC1);
          electron_deltaR1sim.push_back(DeltaRMC1);
	}
       }
       else if(dauId==-11){
	counter++;
        auto p2=d->p4();
	ptmp2=p2;
        //float phi2=d.phi(), eta2=d.eta();
        cout << "   p4_2: "<<p2<<" pt: "<<d->pt()<<endl;
        for(size_t k = 0; k < genParticles->size(); ++ k) {
          const GenParticle & q = (*genParticles)[k];
          auto q4=q.p4();
          float DeltaRMC2=deltaR(p2,q4);
          ////cout<<"deltaRsim: "<<DeltaRMC2<<endl;
	  if(DeltaRMC2>10){DeltaRMC2=9;}
          electron_deltaRsim.push_back(DeltaRMC2);
          electron_deltaR2sim.push_back(DeltaRMC2);
        }
       }
       else if(dauId==23 && id==556){
	counter++;
        auto pz1=d->p4();
	ptmpz1=pz1;
	cout << "   p4_Z1: "<<pz1<<" pt: "<<d->pt()<<endl;
	//int nd = d->numberOfDaughters();
        ////cout<<"\n\nZd#: "<<nd<<"\n\n"<<endl;
        //const GenParticle * dit = d;
	/*for(size_t k = 0; k < genParticles->size(); ++ k) {

	  const GenParticle & dit = (*genParticles)[k];
	  if (dit.pdgId()==23 && dit.pt()==d->pt()){
	     cout<<"genpt: "<<dit.pt()<<" Zd#: "<<dit.numberOfDaughters()<<endl;
	  }

          if(dit->pdgId()!=23){
                cout<<"pdg hija Z: "<<dit->pdgId()<<" ptd: "<<dit->pt()<<endl;

          }
          else{
               	Candidate * dit = dit.daugther(k);
          }
	}*/

       }
       else if(dauId==23 && id==-556){
	counter++;
        auto pz2=d->p4();
	ptmpz2=pz2;
        //cout << "   p4_Z2: "<<pz2<<" pt: "<<d->pt()<<endl;
       }

       if(counter==4){//almacenado de hijas
       cout<<"\nDeltaR1: "<<deltaR(ptmp1,ptmp2)<<' '<<deltaR(ptmp1,ptmpz1)<<' '<<deltaR(ptmp1,ptmpz2)<<endl;
       cout<<"DeltaR2: "<<deltaR(ptmp2,ptmp1)<<' '<<deltaR(ptmp2,ptmpz1)<<' '<<deltaR(ptmp2,ptmpz2)<<"\n"<<endl;

       electron_deltaRtrue.push_back(deltaR(ptmp1,ptmp2));
       electron_deltaRtrue.push_back(deltaR(ptmp1,ptmpz1));
       electron_deltaRtrue.push_back(deltaR(ptmp1,ptmpz2));
       electron_deltaRtrue.push_back(deltaR(ptmp2,ptmp1));
       electron_deltaRtrue.push_back(deltaR(ptmp2,ptmpz1));
       electron_deltaRtrue.push_back(deltaR(ptmp2,ptmpz2));

       electron_deltaR1true.push_back(deltaR(ptmp1,ptmp2));
       electron_deltaR1true.push_back(deltaR(ptmp1,ptmpz1));
       electron_deltaR1true.push_back(deltaR(ptmp1,ptmpz2));

       electron_deltaR2true.push_back(deltaR(ptmp1,ptmp2));
       electron_deltaR2true.push_back(deltaR(ptmp2,ptmpz1));
       electron_deltaR2true.push_back(deltaR(ptmp2,ptmpz2));

       }

      }
     }
   }
}
///Drawer

  if(myelectrons.isValid()){
    // get the number of electrons in the event
    numelectron=myelectrons->size();
    //cout<<"No.Ele: "<<numelectron<<endl;
    for (reco::GsfElectronCollection::const_iterator itElec=myelectrons->begin(); itElec!=myelectrons->end(); ++itElec){

      int missing_hits = itElec->gsfTrack()->trackerExpectedHitsInner().numberOfHits()-itElec->gsfTrack()->hitPattern().numberOfHits();
      bool passelectronveto = !ConversionTools::hasMatchedConversion(*itElec, hConversions, beamspot.position());

      float el_pfIso = 999;
      if (itElec->passingPflowPreselection()) {
	double rho = 0;
	if(rhoHandle.isValid()) rho = *(rhoHandle.product());
	double Aeff = effectiveArea0p3cone(itElec->eta());
	auto iso03 = itElec->pfIsolationVariables();
	el_pfIso = (iso03.chargedHadronIso + std::max(0.0,iso03.neutralHadronIso + iso03.photonIso - rho*Aeff))/itElec->pt();
      }
      auto trk = itElec->gsfTrack();
      bool isLoose = false, isMedium = false, isTight = false;
      if ( abs(itElec->eta()) <= 1.479 ) {
	if ( abs(itElec->deltaEtaSuperClusterTrackAtVtx())<.007 && abs(itElec->deltaPhiSuperClusterTrackAtVtx())<.15 &&
	     itElec->sigmaIetaIeta()<.01 && itElec->hadronicOverEm()<.12 &&
	     abs(trk->dxy(pv))<.02 && abs(trk->dz(pv))<.2 &&
	     missing_hits<=1 && passelectronveto==true &&
	     abs(1/itElec->ecalEnergy()-1/(itElec->ecalEnergy()/itElec->eSuperClusterOverP()))<.05 &&
	     el_pfIso<.15){

	  isLoose = true;

	  if ( abs(itElec->deltaEtaSuperClusterTrackAtVtx())<.004 && abs(itElec->deltaPhiSuperClusterTrackAtVtx())<.06 && abs(trk->dz(pv))<.1 ){
	    isMedium = true;

	    if (abs(itElec->deltaPhiSuperClusterTrackAtVtx())<.03 && missing_hits<=0 && el_pfIso<.10 ){
	      isTight = true;
	    }
	  }
	}
      }
      else if ( abs(itElec->eta()) > 1.479 && abs(itElec->eta()) < 2.5 ) {
	if ( abs(itElec->deltaEtaSuperClusterTrackAtVtx())<.009 && abs(itElec->deltaPhiSuperClusterTrackAtVtx())<.1 &&
	     itElec->sigmaIetaIeta()<.03 && itElec->hadronicOverEm()<.1 &&
	     abs(trk->dxy(pv))<.02 && abs(trk->dz(pv))<.2 &&
	     missing_hits<=1 && el_pfIso<.15 && passelectronveto==true &&
	     abs(1/itElec->ecalEnergy()-1/(itElec->ecalEnergy()/itElec->eSuperClusterOverP()))<.05) {

	  isLoose = true;

	  if ( abs(itElec->deltaEtaSuperClusterTrackAtVtx())<.007 && abs(itElec->deltaPhiSuperClusterTrackAtVtx())<.03 && abs(trk->dz(pv))<.1 ){
	    isMedium = true;

	    if ( abs(itElec->deltaEtaSuperClusterTrackAtVtx())<.005 && abs(itElec->deltaPhiSuperClusterTrackAtVtx())<.02 && missing_hits<=0 && el_pfIso<.10 ){
	      isTight = true;
	    }
	  }
	}
      }

      electron_e.push_back(itElec->energy());
      electron_pt.push_back(itElec->pt());
      electron_px.push_back(itElec->px());
      electron_py.push_back(itElec->py());
      electron_pz.push_back(itElec->pz());
      electron_eta.push_back(itElec->eta());
      electron_phi.push_back(itElec->phi());
      electron_ch.push_back(itElec->charge());
      electron_iso.push_back(el_pfIso);
      electron_isLoose.push_back(isLoose);
      electron_isMedium.push_back(isMedium);
      electron_isTight.push_back(isTight);
      electron_dxy.push_back(trk->dxy(pv));
      electron_dz.push_back(trk->dz(pv));
      electron_dxyError.push_back(trk->d0Error());
      electron_dzError.push_back(trk->dzError());
    }
/////////////////////////////////////Best Gen particle match//////////////////////////////////
if(!isData){
   Handle<GenParticleCollection> genParticles;
   iEvent.getByLabel("genParticles", genParticles);
   vector<GenParticle> genElec;

   //storing selected electrons
   for(size_t i = 0; i < genParticles->size(); ++ i) {

     const GenParticle & p = (*genParticles)[i];
     int id = abs(p.pdgId());
     int st = p.status();
     if(id==11 && st==1){genElec.emplace_back(p);}
   }
   numgenelec=genElec.size();
   //Comparing deltaR GenvsReco
     //cout<<"\n\n No.Reco: "<<myelectrons->size()<<" No.Gen: "<<genElec.size()<<"\n\n"<<endl;

   for (GsfElectronCollection::const_iterator itElec1=myelectrons->begin(); itElec1!=myelectrons->end(); ++itElec1)
   {
     float saveDR=100;
     int idg=-1; //identity gen particle
     for(auto g=genElec.begin(); g!=genElec.end(); g++)
     {
        if(deltaR(g->p4(),itElec1->p4())<saveDR){
       	   saveDR=deltaR(g->p4(),itElec1->p4());
           idg=g-genElec.begin();
        }
     }
     //cout<<"Idex: "<<idg<<" Recop4: "<< itElec1->p4()<<endl;
     if(idg!=-1){
	auto g=genElec.begin()+idg;
	//cout<<"Genp4: "<<g->p4()<<" deltaR: "<<deltaR(g->p4(),itElec1->p4())<<"\n"<<endl;
      	genelec_e.push_back(g->energy());
      	genelec_pt.push_back(g->pt());
      	genelec_px.push_back(g->px());
      	genelec_py.push_back(g->py());
      	genelec_pz.push_back(g->pz());
      	genelec_eta.push_back(g->eta());
      	genelec_phi.push_back(g->phi());
      	genelec_ch.push_back(g->charge());

     }
   }

}

/////////////////////////////////////////End best Gen match//////////////////////////////////

///////////////////////////////////////Secondary Vertices//////////////////////////////////////////////

   Handle<TrackCollection> tracks;
   Handle<reco::TrackCollection> tks;
   ESHandle<TransientTrackBuilder> theB;

   iEvent.getByLabel("generalTracks", tks);
   iEvent.getByLabel("generalTracks", tracks);
/////////////////////TransientTrackBuilding////////////////

      iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
      vector<TransientTrack> t_tks = (*theB).build(tks);
//////////////////////secondary reco////////////////
   Handle<PFJetCollection> myjets;
   iEvent.getByLabel("ak5PFJets", myjets);
//displacement variables
  float dispx;
  float dispy;
  float dispz;
  float disp;

  float dispR;
  float xerr;
  float yerr;
  float zerr;

//iterative identification between track and electron using Delta R, only using tracks pt>20 and electrons from gsf
float saveDR=100;
int identyTrack[myelectrons->size()];
for(size_t x=0; x!=myelectrons->size();x++){identyTrack[x]=-1;}

int k=0;
for (GsfElectronCollection::const_iterator itElec1=myelectrons->begin(); itElec1!=myelectrons->end(); ++itElec1)
{
 if(itElec1->pt()>20)
 {
  int j=0;
  for(TrackCollection::const_iterator itTrack1 = tracks->begin(); itTrack1 != tracks->end(); ++itTrack1)
  {
     if(itTrack1->pt()>20)
     {
	float DR_ET = deltaR(itElec1->eta(),itElec1->phi(),itTrack1->eta(),itTrack1->phi());
	if(DR_ET<saveDR){
	  saveDR=DR_ET;
	  identyTrack[k]=j;
        }
        j++;
     }
  }
 }
 k++;
 electron_Bsecvec.push_back(-1);//initialyzing best secondary vertex vector
}
//fin identifier

vector<float> savedisp (myelectrons->size());
//float savedisp[myelectrons.size()];
int  i=0;
for(TrackCollection::const_iterator itTrack1 = tracks->begin();
       itTrack1 != tracks->end();
       ++itTrack1)
{
  if( itTrack1->pt()>20 ){
    int	j=0;
    for(TrackCollection::const_iterator itTrack2 = tracks->begin();
       itTrack2 != tracks->end();
       ++itTrack2)
       {

       if( itTrack2->pt()>10 ){
	int j3=0;
	for(TrackCollection::const_iterator itTrack3 = tracks->begin(); itTrack3 != tracks->end(); ++itTrack3)
  	{
	if( itTrack3->pt()>10 ){

         KalmanVertexFitter fitter;
         vector<TransientTrack> trackVec;

         if(t_tks.size()>2 && itTrack1->pt()!=itTrack2->pt()){
	  if(itTrack1->pt()!=itTrack3->pt() && itTrack2->pt()!=itTrack3->pt()){

	   /*cout<<"\npt1: "<<itTrack1->pt()<<" pt2: "<<itTrack2->pt()<<" pt3: "<<itTrack3->pt();
	   cout<<" deltaR1-2: "<<deltaR(itTrack1->phi(),itTrack1->eta(),itTrack2->phi(),itTrack2->eta());
	   cout<<" deltaR1-3: "<<deltaR(itTrack3->phi(),itTrack3->eta(),itTrack1->phi(),itTrack1->eta())<<endl;*/

           //auto trk1 = itElec1->gsfTrack();
           TransientTrack t_trk1 = (* theB).build(* itTrack1);

           //auto trk2 = itElec2->gsfTrack();
           TransientTrack t_trk2 = (* theB).build(* itTrack2);
	   TransientTrack t_trk3 = (* theB).build(* itTrack3);

           trackVec.push_back(t_trk1);
           trackVec.push_back(t_trk2);
	   trackVec.push_back(t_trk3);
           TransientVertex myVertex = fitter.vertex(trackVec);//reconstruction of secondary vertex Sometimes
	   trackVec.clear();

        int k=0;
       if(myVertex.isValid()){
	numsecvec++;
	secvec_posx.push_back(myVertex.position().x());
  	secvec_posy.push_back(myVertex.position().y());
  	secvec_posz.push_back(myVertex.position().z());
  	secvec_poserrorx.push_back(myVertex.positionError().cxx());
  	secvec_poserrory.push_back(myVertex.positionError().cyy());
  	secvec_poserrorz.push_back(myVertex.positionError().czz());
	for (GsfElectronCollection::const_iterator itElec1=myelectrons->begin(); itElec1!=myelectrons->end(); ++itElec1)
	{
             numdisp++;//total of secondary vertex disp per event
	     k++;//# of saved disp per SecondaryVertex

	     //primaryvertex
	     //cout<<"ptE: "<<itElec1->pt()<<endl;
             //cout<<"pvx: "<< itElec1->vx()<<" pvy: "<<itElec1->vy()<<" pvz: "<< itElec1->vz()<<endl;
             //cout<<"phiE: "<<itElec1->phi()<<" etaE: "<<itElec1->eta()<<" phiT: "<<itTrack1->phi()<<" etaT: "<<itTrack1->eta()<<endl;
	     //float deltaR_ET=deltaR(itElec1->phi(),itElec1->eta(),itTrack1->phi(),itTrack1->eta());
	     //cout<<"deltaR-ET: "<<deltaR_ET<<endl;

             //secondaryvertex
             //cout<<"PosX: "<< myVertex.position().x()<<" PosY: "<<myVertex.position().y()<<" PosZ: "<<myVertex.position().z()<<endl;
             dispx= itElec1->vx()-myVertex.position().x();
             dispy= itElec1->vy()-myVertex.position().y();
             dispz= itElec1->vz()-myVertex.position().z();

             //error calculation + displacement
             xerr=myVertex.positionError().cxx();
             yerr=myVertex.positionError().cyy();
             zerr=myVertex.positionError().czz();

             float difx = (myVertex.position().x())/(sqrt((xerr*xerr)+(yerr*yerr)+(zerr*zerr)));
             float dify = (myVertex.position().y())/(sqrt((xerr*xerr)+(yerr*yerr)+(zerr*zerr)));
             float difz = (myVertex.position().z())/(sqrt((xerr*xerr)+(yerr*yerr)+(zerr*zerr)));
             float err  = sqrt(difx*difx*xerr+dify*dify*yerr+difz*difz*zerr);

             disp= sqrt(dispx*dispx + dispy*dispy + dispz*dispz);
             //cout<<"Displacement: "<<disp<<endl;
             dispR= sqrt(dispx*dispx + dispy*dispy + dispz*dispz)/err;
	     //cout<<"Weigthed Displacement: "<<dispR<<'\n'<<endl;


             //displacements store
	     secvec_disp.push_back(disp);
             float totDR=1.5;
             if(dispR<totDR){secvec_dispR.push_back(dispR);}
             else{secvec_dispR.push_back(totDR);}

             //DeltaR identification
             float DeltaRPrima = deltaR(itElec1->eta(),itElec1->phi(),itTrack2->eta(),itTrack2->phi());
	     //cout<<"DeltaR: "<<DeltaRPrima<<'\n'<<endl;

             secvec_deltaR.push_back(DeltaRPrima);
             if(i==0){secvec_deltaR1.push_back(DeltaRPrima);}
             else if(i==1){secvec_deltaR2.push_back(DeltaRPrima);}

	     //find best secvec
             if(identyTrack[k-1]==i)
	     {
		if(dispR > savedisp.at(k-1))
		{
		  savedisp.at(k-1)=dispR;
		  electron_Bsecvec.at(k-1)=numdisp;
		}
	     }

	 }// for itElec

        }
	}//nonequal Tracks condition 3-2 3-1
       }//Nonequal Tracks condition 1-2
      }
      j3++;
      }//for Track3
     }
     j++;
    }
    i++;
 }//finalde if pt track
    electron_secN.push_back(i);
}//fin for Track1
}// Final de Is ValidElectron

  mtree->Fill();
  return;

}

// ------------ method called once each job just before starting event loop  ------------
void
ElectronAnalyzer::beginJob()
{}

// ------------ method called once each job just after ending the event loop  ------------
void
ElectronAnalyzer::endJob()
{}

// ------------ method called when starting to processes a run  ------------
void
ElectronAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a run  ------------
void
ElectronAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{}
// ------------ method called when starting to processes a luminosity block  ------------
void
ElectronAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a luminosity block  ------------
void
ElectronAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ElectronAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
 //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronAnalyzer);
