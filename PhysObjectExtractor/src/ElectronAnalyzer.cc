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

//Trigg object match
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

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
  std::string filterName_;

  // ----------member data ---------------------------
  bool isData;

  TTree *mtree;
  int numelectron; //number of reco electrons in the event
  int numgenelec; //# of gen electrons in the event
  int numsecvec;  //number of secondary vertex disp in the event
  int numEsecvec;  //number of secondary vertex disp in the event
  int numdisp;  //tool for finding best match for electrons
  std::vector<float> electron_e;
  std::vector<float> electron_pt;
  std::vector<float> electron_px;
  std::vector<float> electron_py;
  std::vector<float> electron_pz;
  std::vector<float> electron_vx;
  std::vector<float> electron_vy;
  std::vector<float> electron_vz;
  //std::vector<float> electron_vxerr;
  //std::vector<float> electron_vyerr;
  //std::vector<float> electron_vzerr;
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
  std::vector<float> genelec_DRscore;
  std::vector<float> electron_iso;
  std::vector<bool> electron_isLoose;
  std::vector<bool> electron_isMedium;
  std::vector<bool> electron_isTight;
  std::vector<int> trk_identity;
  std::vector<int> trk_isHQ;
  std::vector<float> trk_NChi2;
  std::vector<float> trk_d0;
  std::vector<float> trk_d0E;
  std::vector<float> trk_dZ;
  std::vector<float> trk_dZE;
  std::vector<int> trk_Nlayer;
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
  std::vector<int> secvec_eleTag;
  std::vector<float> secvec_chi2;
  std::vector<float> secvec_nodf;
  std::vector<float> secvec_normchi2;

  std::vector<float> Esecvec_posx;
  std::vector<float> Esecvec_posy;
  std::vector<float> Esecvec_posz;
  std::vector<float> Esecvec_poserrorx;
  std::vector<float> Esecvec_poserrory;
  std::vector<float> Esecvec_poserrorz;
  std::vector<int> Esecvec_eleTag;
  std::vector<float> Esecvec_chi2;
  std::vector<float> Esecvec_nodf;
  std::vector<float> Esecvec_normchi2;
  //std::vector<int> electron_Bsecvec;//best match for electron sec vec
  std::vector<int> electron_BdR;//best match delta R for electron sec vec
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
  //filterName_ = iConfig.getParameter<std::string>("filterName");
  edm::Service<TFileService> fs;
  mtree = fs->make<TTree>("Events", "Events");


  mtree->Branch("numberelectron",&numelectron);
  mtree->GetBranch("numberelectron")->SetTitle("number of reco electrons");
  mtree->Branch("numbergenelec",&numgenelec);
  mtree->GetBranch("numbergenelec")->SetTitle("number of gen electrons");
  mtree->Branch("numbersecvec",&numsecvec);
  mtree->GetBranch("numbersecvec")->SetTitle("number of secondary vertex displascements");
  mtree->Branch("numberEsecvec",&numEsecvec);
  mtree->GetBranch("numberEsecvec")->SetTitle("number of secondary vertex displascements");
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
  mtree->Branch("electron_vx",&electron_vx);
  mtree->GetBranch("electron_vx")->SetTitle("electron vertex x-component");
  mtree->Branch("electron_vy",&electron_vy);
  mtree->GetBranch("electron_vy")->SetTitle("electron vertex y-component");
  mtree->Branch("electron_vz",&electron_vz);
  mtree->GetBranch("electron_vz")->SetTitle("electron vertex z-component");
  //mtree->Branch("electron_vxerr",&electron_vxerr);
  //mtree->GetBranch("electron_vxerr")->SetTitle("electron vertex x-component error");
  //mtree->Branch("electron_vyerr",&electron_vyerr);
  //mtree->GetBranch("electron_vyerr")->SetTitle("electron vertex y-component error");
  //mtree->Branch("electron_vzerr",&electron_vzerr);
  //mtree->GetBranch("electron_vzerr")->SetTitle("electron vertex z-component error");
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
  mtree->Branch("genelec_DRscore",&genelec_DRscore);
  mtree->GetBranch("genelec_DRscore")->SetTitle("electron delta R score obtained during fitting");
  mtree->Branch("electron_iso",&electron_iso);
  mtree->GetBranch("electron_iso")->SetTitle("electron isolation");
  mtree->Branch("electron_isLoose",&electron_isLoose);
  mtree->GetBranch("electron_isLoose")->SetTitle("electron tagged loose");
  mtree->Branch("electron_isMedium",&electron_isMedium);
  mtree->GetBranch("electron_isMedium")->SetTitle("electron tagged medium");
  mtree->Branch("electron_isTight",&electron_isTight);
  mtree->GetBranch("electron_isTight")->SetTitle("electron tagged tight");
  mtree->Branch("trk_identity",&trk_identity);
  mtree->GetBranch("trk_identity")->SetTitle("electron track identity tag");
  mtree->Branch("trk_isHQ",&trk_isHQ);
  mtree->GetBranch("trk_isHQ")->SetTitle("electron track tagged HighQuality");
  mtree->Branch("trk_NChi2",&trk_NChi2);
  mtree->GetBranch("trk_NChi2")->SetTitle("electron track Normalized Chi squared (mm)");
  mtree->Branch("trk_d0",&trk_d0);
  mtree->GetBranch("trk_d0")->SetTitle("electron track transverse distance to beamspot (mm)");
  mtree->Branch("trk_d0E",&trk_d0E);
  mtree->GetBranch("trk_d0E")->SetTitle("electron track transverse distance to beamspot error (mm)");
  mtree->Branch("trk_dZ",&trk_dZ);
  mtree->GetBranch("trk_dZ")->SetTitle("electron track longitudinal distance to beamspot (mm)");
  mtree->Branch("trk_dZE",&trk_dZE);
  mtree->GetBranch("trk_dZE")->SetTitle("electron track longitudinal distance to beamspot error (mm)");
  mtree->Branch("trk_Nlayer",&trk_Nlayer);
  mtree->GetBranch("trk_Nlayer")->SetTitle("electron track number of hit layers (mm)");
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
  mtree->Branch("secvec_eleTag",&secvec_eleTag);
  mtree->GetBranch("secvec_eleTag")->SetTitle("secvec electron Tag");
  mtree->Branch("secvec_chi2",&secvec_chi2);
  mtree->GetBranch("secvec_chi2")->SetTitle("Vertex chi squared");
  mtree->Branch("secvec_normchi2",&secvec_normchi2);
  mtree->GetBranch("secvec_normchi2")->SetTitle("Vertex normalized chi squared");
  mtree->Branch("secvec_nodf",&secvec_nodf);
  mtree->GetBranch("secvec_nodf")->SetTitle("Track number of degree of freedom");

  mtree->Branch("Esecvec_posx",&Esecvec_posx);
  mtree->GetBranch("Esecvec_posx")->SetTitle("secvec position x (mm)");
  mtree->Branch("Esecvec_posy",&Esecvec_posy);
  mtree->GetBranch("Esecvec_posy")->SetTitle("secvec position y (mm)");
  mtree->Branch("Esecvec_posz",&Esecvec_posz);
  mtree->GetBranch("Esecvec_posz")->SetTitle("secvec position z (mm)");
  mtree->Branch("Esecvec_poserrorx",&Esecvec_poserrorx);
  mtree->GetBranch("Esecvec_poserrorx")->SetTitle("secvec position x error (mm)");
  mtree->Branch("Esecvec_poserrory",&Esecvec_poserrory);
  mtree->GetBranch("Esecvec_poserrory")->SetTitle("secvec position y error (mm)");
  mtree->Branch("Esecvec_poserrorz",&Esecvec_poserrorz);
  mtree->GetBranch("Esecvec_poserrorz")->SetTitle("secvec position z error (mm)");
  mtree->Branch("Esecvec_eleTag",&Esecvec_eleTag);
  mtree->GetBranch("Esecvec_eleTag")->SetTitle("secvec electron Tag");
  mtree->Branch("Esecvec_chi2",&Esecvec_chi2);
  mtree->GetBranch("Esecvec_chi2")->SetTitle("Vertex chi squared");
  mtree->Branch("Esecvec_normchi2",&Esecvec_normchi2);
  mtree->GetBranch("Esecvec_normchi2")->SetTitle("Vertex normalized chi squared");
  mtree->Branch("Esecvec_nodf",&Esecvec_nodf);
  mtree->GetBranch("Esecvec_nodf")->SetTitle("Track number of degree of freedom");
  /*mtree->Branch("electron_Bsecvec",&electron_Bsecvec);
  mtree->GetBranch("electron_Bsecvec")->SetTitle("best match electrons for secondary vertex");*/
  mtree->Branch("electron_BdR",&electron_BdR);
  mtree->GetBranch("electron_BdR")->SetTitle("best match delta R electrons for secondary vertex");
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
  numEsecvec = 0;
  numdisp = 0;
  electron_e.clear();
  electron_pt.clear();
  electron_px.clear();
  electron_py.clear();
  electron_pz.clear();
  electron_vx.clear();
  electron_vy.clear();
  electron_vz.clear();
  //electron_vxerr.clear();
  //electron_vyerr.clear();
  //electron_vzerr.clear();
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
  genelec_DRscore.clear();
  electron_iso.clear();
  electron_isLoose.clear();
  electron_isMedium.clear();
  electron_isTight.clear();
  trk_identity.clear();
  trk_isHQ.clear();
  trk_NChi2.clear();
  trk_d0.clear();
  trk_d0E.clear();
  trk_dZ.clear();
  trk_dZE.clear();
  trk_Nlayer.clear();
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
  secvec_eleTag.clear();
  secvec_chi2.clear();
  secvec_nodf.clear();
  secvec_normchi2.clear();

  Esecvec_posx.clear();
  Esecvec_posy.clear();
  Esecvec_posz.clear();
  Esecvec_poserrorx.clear();
  Esecvec_poserrory.clear();
  Esecvec_poserrorz.clear();
  Esecvec_eleTag.clear();
  Esecvec_chi2.clear();
  Esecvec_nodf.clear();
  Esecvec_normchi2.clear();
  //electron_Bsecvec.clear();
  electron_BdR.clear();
  electron_secN.clear();
  electron_deltaRsim.clear();
  electron_deltaR1sim.clear();
  electron_deltaR2sim.clear();
  electron_deltaRtrue.clear();
  electron_deltaR1true.clear();
  electron_deltaR2true.clear();
  Zjet_pt.clear();

  if(myelectrons.isValid()){
    // get the number of electrons in the event
    //numelectron=myelectrons->size();
    //cout<<"No.Ele: "<<numelectron<<endl;
    for (reco::GsfElectronCollection::const_iterator itElec=myelectrons->begin(); itElec!=myelectrons->end(); ++itElec){
     if(itElec->pt()>-1){
      numelectron++;//Here get actual number of electrons

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
      //cout<<trk->quality(reco::Track::highPurity)<<endl;
      //cout<<"Correspondance:"<<itElec->p4()<<endl;
      electron_e.push_back(itElec->energy());
      electron_pt.push_back(itElec->pt());
      electron_px.push_back(itElec->px());
      electron_py.push_back(itElec->py());
      electron_pz.push_back(itElec->pz());
      electron_px.push_back(itElec->vx());
      electron_py.push_back(itElec->vy());
      electron_pz.push_back(itElec->vz());
      //electron_px.push_back(itElec->vxerror());
      //electron_py.push_back(itElec->vyerror());
      //electron_pz.push_back(itElec->pzvzerror());
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


     }//fin if pt elec
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
    if(itElec1->pt()>-1){
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
	genelec_DRscore.push_back(deltaR(g->p4(),itElec1->p4()));
     }
    }//fin If pt elec
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
   //Handle<PFJetCollection> myjets;//jets collection if neede
   //iEvent.getByLabel("ak5PFJets", myjets);
//displacement variables
  /*float dispx;
  float dispy;
  float dispz;
  float disp;

  float dispR;
  float xerr;
  float yerr;
  float zerr;*/

//iterative identification between track and electron using Delta R, only using tracks pt>20 and electrons from gsf
float saveDR;
int identyTrack[myelectrons->size()];
for(size_t x=0; x!=myelectrons->size();x++){identyTrack[x]=-1;}

int k=0;
//cout<<"size "<<myelectrons->size()<<endl;
for (GsfElectronCollection::const_iterator itElec1=myelectrons->begin(); itElec1!=myelectrons->end(); ++itElec1)
{
 saveDR=100;
 //if(itElec1->pt()>20) {
  int j=0;
  for(TrackCollection::const_iterator itTrack1 = tracks->begin(); itTrack1 != tracks->end(); ++itTrack1)
  {
     if(itTrack1->pt()>1){
        //cout<<"check\n";
      	float DR_ET = deltaR(itElec1->eta(),itElec1->phi(),itTrack1->eta(),itTrack1->phi());
      	if(DR_ET<saveDR){
      	  saveDR=DR_ET;
      	  identyTrack[k]=j;//track idex thta best fit Electrons
        }
        j++;
     }
  }
 //}


//////////////////////////////////////////////////////////////////////

 if(identyTrack[k]!=-1){
  //cout<<"   "<<identyTrack[k]<<' '<<itElec1->pt()<<' ';
  int j=0;
  for(TrackCollection::const_iterator itTrack1 = tracks->begin(); itTrack1 != tracks->end(); ++itTrack1)
  {
     if(itTrack1->pt()>1){
      	if(j==identyTrack[k]){
      	   /*auto trk1=itElec1->gsfTrack();
      	   auto DRtest=deltaR(itElec1->eta(),itElec1->phi(),itTrack1->eta(),itTrack1->phi());
      	   cout<<trk1->pt()<<' '<<j<<' '<< itTrack1->pt()<<' '<<DRtest<<endl;
      	   cout<<"   "<<trk1->quality(reco::Track::highPurity)<<' '<<itTrack1->quality(reco::Track::highPurity);
      	   cout<<' '<<itTrack1->normalizedChi2()<<' ';
      	   cout<<itTrack1->dxy(beamspot.position())<<' '<<itTrack1->dz(beamspot.position())<<' '<<itTrack1->d0Error()<<endl;*/
      	   trk_isHQ.push_back(itTrack1->quality(reco::Track::highPurity));
                 trk_NChi2.push_back(itTrack1->normalizedChi2());
      	   trk_d0.push_back(itTrack1->dxy(beamspot.position()));
      	   trk_d0E.push_back(itTrack1->d0Error());
                 trk_dZ.push_back(itTrack1->dz(beamspot.position()));
                 trk_dZE.push_back(itTrack1->dzError());
      	   trk_Nlayer.push_back(itTrack1->hitPattern().trackerLayersWithMeasurement());
         }

	      }
	      j++;
     //}
  }
 }
 else {
   trk_isHQ.push_back(-1);
         trk_NChi2.push_back(-1);
   trk_d0.push_back(-1);
   trk_d0E.push_back(-1);
         trk_dZ.push_back(-1);
         trk_dZE.push_back(-1);
   trk_Nlayer.push_back(-1);
 }
 k++;
 electron_BdR.push_back(saveDR);

 //electron_Bsecvec.push_back(-1);//initialyzing best secondary vertex vector
}

for(size_t x =0; x<myelectrons->size(); x++){trk_identity.push_back(identyTrack[x]);}

////////////////////////////////////////////////////////////////////



//fin identifier

vector<float> savedisp (myelectrons->size());
vector<TransientVertex> myVertices;
KalmanVertexFitter fitter;
//float savedisp[myelectrons.size()];
int  i=0;
for(TrackCollection::const_iterator itTrack1 = tracks->begin();
       itTrack1 != tracks->end();
       ++itTrack1)
{
  if(itTrack1->pt()>1 && itTrack1->quality(reco::Track::highPurity) ){
      //cout<<itTrack1->quality(reco::Track::highPurity)<<endl;
      int	j2=0;
      for(TrackCollection::const_iterator itTrack2 = itTrack1+1;
         itTrack2 != tracks->end();
         ++itTrack2)
         {

           if( itTrack2->pt()>1 && itTrack2->quality(reco::Track::highPurity) ){

             vector<TransientTrack> trackVec;

             if(t_tks.size()>2 && itTrack1->pt()!=itTrack2->pt()){
          	  //if(itTrack1->pt()!=itTrack3->pt() && itTrack2->pt()!=itTrack3->pt()){

          	   /*cout<<"\npt1: "<<itTrack1->pt()<<" pt2: "<<itTrack2->pt()<<" pt3: "<<itTrack3->pt();
          	   cout<<" deltaR1-2: "<<deltaR(itTrack1->phi(),itTrack1->eta(),itTrack2->phi(),itTrack2->eta());
          	   cout<<" deltaR1-3: "<<deltaR(itTrack3->phi(),itTrack3->eta(),itTrack1->phi(),itTrack1->eta())<<endl;*/

                     //auto trk1 = itElec1->gsfTrack();
                     TransientTrack t_trk1 = (* theB).build(* itTrack1);

                     //auto trk2 = itElec2->gsfTrack();
                     TransientTrack t_trk2 = (* theB).build(* itTrack2);
          	         //TransientTrack t_trk3 = (* theB).build(* itTrack3);

                     trackVec.push_back(t_trk1);
                     trackVec.push_back(t_trk2);
          	         //trackVec.push_back(t_trk3);
                     TransientVertex myVertex = fitter.vertex(trackVec);//reconstruction of secondary vertex Sometimes
          	         trackVec.clear();

                  //int k=0;//este K es para el numero de desplazamientos
                 if(myVertex.isValid()){
                    myVertices.push_back(myVertex);
                  }
          	//}//nonequal Tracks condition 3-2 3-1
           }//Nonequal Tracks condition 1-2
         }
         j2++;
      }
      i++;
  }//final if pt track HQ track
}//fin for Track1
electron_secN.push_back(i);

//cout<<"\n\n\n Number of vertices before merging: "<<myVertices.size()<<"\n\n\n";

///////////////////////Merging vertices////////////////////////////////
vector<TransientVertex> tmpVertices;//temporal vertices storage use for iterations
vector<TransientVertex> finalVertices;//Vertices Fully merged and filtered
vector<TransientTrack> Otrk1;//Original tracks stoage use for iterations
vector<TransientTrack> Otrk2;//Original tracks stoage use for iterations
vector<TransientTrack> Otrktmp;//Merged Original tracks storage use for Vertex fitting
vector<size_t> Vfound;//Storage for vertices that were already compared and similed
vector<size_t> Tfound;//Storage for tracks that were compared and similed
bool rSimil=true;//flag whenever a two similar tracks where found if not found any end merging
bool Vf=false;//refered to Vfound flag tells if the iterand vertex aws already compared
bool Similf=false;// refered to Tfound flag tells if the iterand track is repeated between two vertices. Avoid reapeting tracks in Otrktmp
size_t its=0;// counts number of iterations
size_t vertexcount=0;//counts number of merged vertices
finalVertices.clear();

while(rSimil && its<50){// lasso while stops whenever there are not vertices left to merge or reach the mx number of Iterations
  /////Initialize storage and flags////
  tmpVertices=myVertices;
  myVertices.clear();
  vertexcount=0;
  rSimil=false;
  Vfound.clear();
  ///////////////////////////////////
  its++;//count iterations
  for(size_t x=0; x<tmpVertices.size(); x++){////iteration over the all vertices/////
    //cout<<"1st it flag "<<x<<endl;
    bool Vsimil=false;///// Flag tells if the vertex share tracks with another. guarantees only storing not mergeable vertices

    //bool Vrepeat=false;
    size_t repeatcount=0;//// Counter tells how many times a vertex is repeting in the storage. I.e. whenever two vertices have the same number of tracks and all are equal.
    Vf=false;//initialize
    Otrk1.clear();//initialize

    for(size_t i=0; i<Vfound.size(); i++){if(Vfound.at(i)==x)Vf=true;}//omit comparing a vertex that was already merged

    if(tmpVertices.at(x).normalisedChiSquared()<5 && !Vf){// quality check plus already merged vertex check
      Otrk1=tmpVertices.at(x).originalTracks();// extract original tracks from first compared vertex
      for(size_t y=x+1; y<tmpVertices.size(); y++){/// iteration over all vertices
        //cout<<" 2nd it Flag "<<y<<endl;
        size_t similcount=0;//initialize
        Tfound.clear();//Restart Found tracks

        Vf=false;//Restart alread merged flag
        Otrk2.clear();// restart track container

        for(size_t i=0; i<Vfound.size(); i++){if(Vfound.at(i)==y)Vf=true;}//tell if the current vertex was already compared and merged

        if(tmpVertices.at(y).normalisedChiSquared()<5 && !Vf){ //quality check plus already merged check
          Otrktmp.clear();//restart Temporal track Storage
          Otrk2=tmpVertices.at(y).originalTracks();//Extract original tracks from second compared vertex
          for(size_t z=0; z<Otrk1.size(); z++){// first track set iterations
            for(size_t v=0; v<Otrk2.size(); v++){// second track set iterations
              if( deltaR(Otrk1.at(z).track().eta(),Otrk1.at(z).track().phi(),Otrk2.at(v).track().eta(),Otrk2.at(v).track().phi())==0
                 && Otrk1.at(z).track().pt()==Otrk2.at(v).track().pt() ){ //check whenever a track is shared in both sets
                   //cout<<"X "<<x<<" Y "<<y<<" Z "<<z<<" V "<<v<<endl;
                   //cout<<"  "<<deltaR(Otrk1.at(z).track().eta(),Otrk1.at(z).track().phi(),Otrk2.at(v).track().eta(),Otrk2.at(v).track().phi())
                   //<<"  "<<Otrk1.at(z).track().pt()<<"  "<<Otrk2.at(v).track().pt()<<endl;
                   similcount++;// count number of similar tracks from each set
                   rSimil=true;// FLag there are mergeable vertices
                   //breakSimil=true;
                   Otrktmp=Otrk1;// store all track from first set
                   Vfound.push_back(y);// store index of already compared and merged vertex
                   Tfound.push_back(v);// store index of already compared and similar track
                   Vsimil=true;//Flag the compared vertex is mergeable
                   //break;
                 }
            }
            //if(breakSimil) break;
          }
          if(similcount!=0){//check if two vertices share tracks
            for(size_t w=0; w<Otrk2.size(); w++){//iteration over second set of tracks
              Similf=false;//reset flag track is repeated
              for(size_t s=0; s<Tfound.size(); s++){if(w==Tfound.at(s)) Similf=true;}//checks flag if track is repeated
              if(!Similf) Otrktmp.push_back(Otrk2.at(w));// if track is not repeated store in temporal track storage
            }
          }
          //cout<<" Temporal track size: "<<Otrktmp.size()<<" 1st vertex iterator index: "<<x<<" 2nd vertex iterator index: "<< y <<endl;
          ////poner aqui el Fitter
          //cout<<"Simil Tracks: "<<similcount<<" Otrk1 "<< Otrk1.size()<<" Otrk2 "<<Otrk2.size()<<endl;
          //cout<<" Temporal track size: "<<Otrktmp.size()<<" 1st vertex iterator index: "<<x<<" 2nd vertex iterator index: "<< y <<endl;
          if(similcount==Otrk1.size() && Otrk1.size()==Otrk2.size())repeatcount++;// if two or more vertex are identical count
          /*if(similcount==Otrk1.size() && Otrk1.size()==Otrk2.size()){
            break;
          }*/
          if(Otrktmp.size()>2 && repeatcount<2){//check if temporal track could be fitted check if the identical vertex where already merged
            //cout<<"Simil Tracks: "<<similcount<<" Otrk1 "<< Otrk1.size()<<" Otrk2 "<<Otrk2.size()<<endl;
            //cout<<" Temporal track size: "<<Otrktmp.size()<<" 1st vertex iterator index: "<<x<<" 2nd vertex iterator index: "<< y <<endl;
            TransientVertex myVertex = fitter.vertex(Otrktmp);//fiit vertex
            if(myVertex.isValid()){// sanity check
              vertexcount++;// count number of vertices
              myVertices.push_back(myVertex);//restore merged vertices
            }
          }
        }
      }
      if(!Vsimil && tmpVertices.at(x).normalisedChiSquared()<5) finalVertices.push_back(tmpVertices.at(x));// if track is not mergeable and is good store as final vertex
      //cout<<"Simil vertices: "<<similcount<<endl;
    }
  }
  //cout<<" Valid vertex "<< vertexcount<<endl;
  //cout<<"\n\n\n Number of vertices after merging: "<<myVertices.size()<<" Final vertices: "<<finalVertices.size()<<"\n\n\n";
}
//cout<<"Iterations: "<< its <<" Final vertices: "<<finalVertices.size()<<endl;
/////////////Mergin fin////////////////////////
for(size_t x=0; x<finalVertices.size(); x++){
  TransientVertex myV=finalVertices.at(x);
  Otrk1.clear();
  numsecvec++;
  secvec_posx.push_back(myV.position().x());
  secvec_posy.push_back(myV.position().y());
  secvec_posz.push_back(myV.position().z());
  secvec_poserrorx.push_back(myV.positionError().cxx());
  secvec_poserrory.push_back(myV.positionError().cyy());
  secvec_poserrorz.push_back(myV.positionError().czz());
  secvec_chi2.push_back(myV.totalChiSquared());
  secvec_nodf.push_back(myV.degreesOfFreedom());
  secvec_normchi2.push_back(myV.normalisedChiSquared());
  secvec_eleTag.push_back(-1);

  Otrk1=myV.originalTracks();

  //cout<<"Vertex number of tracks "<<Otrk1.size()<<endl;
  int k=0;
  for(TrackCollection::const_iterator itTrack1 = tracks->begin(); itTrack1 != tracks->end(); ++itTrack1){
    //bool trackId=false;
    if(itTrack1->pt()>1){
      k++;
      for(size_t y=0; y!=myelectrons->size();y++){
        if(identyTrack[y]==k){
          for(size_t z=0; z<Otrk1.size(); z++){
            if(deltaR(Otrk1.at(z).track().eta(),Otrk1.at(z).track().phi(),itTrack1->eta(),itTrack1->phi())==0
              && Otrk1.at(z).track().pt()==itTrack1->pt()){
                secvec_eleTag.back()=y;
              }
          }
        }
      }

    }
  }


  //Delta R entre Track1 Y electron en electron_BdR
  /*secvec_phi.push_back(itTrack1->phi());
  secvec_eta.push_back(itTrack1->eta());
  secvec_deltaR1.push_back(deltaR(itTrack1->eta(),itTrack1->phi(),itTrack2->eta(),itTrack2->phi()));
  secvec_phi1.push_back(itTrack2->phi());
  secvec_eta1.push_back(itTrack2->eta());
  secvec_deltaR2.push_back(deltaR(itTrack1->eta(),itTrack1->phi(),itTrack3->eta(),itTrack3->phi()));
  secvec_phi2.push_back(itTrack3->phi());
  secvec_eta2.push_back(itTrack3->eta());*/
}
//cout<<"\n Nsize: "<<secvec_chi2.size()<<" Tag size: "<<secvec_eleTag.size()<<" electron size "<<myelectrons->size()<<"\n\n";



//////////////////////////////////////////Making secondary vertices with electron tracks///////////////////////////////
///////////////plan iterate all electron track with hq tracks to form vertices and then merge them all
bool Eletrack;
int EleCount=0;
for(GsfElectronCollection::const_iterator itElec1=myelectrons->begin(); itElec1!=myelectrons->end(); ++itElec1){
  Eletrack=false;
  i=0;
  for(TrackCollection::const_iterator itTrack1 = tracks->begin();
         itTrack1 != tracks->end();
         ++itTrack1)
  {
    if( itTrack1->pt()>1){
      if(i==identyTrack[EleCount]) Eletrack=true;
      if( Eletrack && itTrack1->quality(reco::Track::highPurity) ){
        //cout<<itTrack1->quality(reco::Track::highPurity)<<endl;
        int	j2=0;
        for(TrackCollection::const_iterator itTrack2 = itTrack1+1;
           itTrack2 != tracks->end();
           ++itTrack2)
           {

             if( itTrack2->pt()>1 && itTrack2->quality(reco::Track::highPurity) ){

               vector<TransientTrack> trackVec;

               if(t_tks.size()>2 && itTrack1->pt()!=itTrack2->pt()){
                //if(itTrack1->pt()!=itTrack3->pt() && itTrack2->pt()!=itTrack3->pt()){

                 /*cout<<"\npt1: "<<itTrack1->pt()<<" pt2: "<<itTrack2->pt()<<" pt3: "<<itTrack3->pt();
                 cout<<" deltaR1-2: "<<deltaR(itTrack1->phi(),itTrack1->eta(),itTrack2->phi(),itTrack2->eta());
                 cout<<" deltaR1-3: "<<deltaR(itTrack3->phi(),itTrack3->eta(),itTrack1->phi(),itTrack1->eta())<<endl;*/

                       //auto trk1 = itElec1->gsfTrack();
                       TransientTrack t_trk1 = (* theB).build(* itTrack1);

                       //auto trk2 = itElec2->gsfTrack();
                       TransientTrack t_trk2 = (* theB).build(* itTrack2);
                       //TransientTrack t_trk3 = (* theB).build(* itTrack3);

                       trackVec.push_back(t_trk1);
                       trackVec.push_back(t_trk2);
                       //trackVec.push_back(t_trk3);
                       TransientVertex myVertex = fitter.vertex(trackVec);//reconstruction of secondary vertex Sometimes
                       trackVec.clear();

                    //int k=0;//este K es para el numero de desplazamientos
                   if(myVertex.isValid()){
                      myVertices.push_back(myVertex);
                    }
              //}//nonequal Tracks condition 3-2 3-1
             }//Nonequal Tracks condition 1-2
           }
           j2++;
        }
      }
        i++;
    }//final if pt track HQ track
  }//fin for Track1
  electron_secN.push_back(i);
  cout<<"\n\n\n Number of Ele vertices before merging: "<<myVertices.size()<<"\n\n\n";
  ///MErging electrons Vertices
  rSimil=true;//flag whenever a two similar tracks where found if not found any end merging
  Vf=false;//refered to Vfound flag tells if the iterand vertex aws already compared
  Similf=false;// refered to Tfound flag tells if the iterand track is repeated between two vertices. Avoid reapeting tracks in Otrktmp
  its=0;// counts number of iterations
  vertexcount=0;//counts number of merged vertices
  finalVertices.clear();

  while(rSimil && its<50){// lasso while stops whenever there are not vertices left to merge or reach the mx number of Iterations
    /////Initialize storage and flags////
    tmpVertices=myVertices;
    myVertices.clear();
    vertexcount=0;
    rSimil=false;
    Vfound.clear();
    ///////////////////////////////////
    its++;//count iterations
    for(size_t x=0; x<tmpVertices.size(); x++){////iteration over the all vertices/////
      //cout<<"1st it flag "<<x<<endl;
      bool Vsimil=false;///// Flag tells if the vertex share tracks with another. guarantees only storing not mergeable vertices

      //bool Vrepeat=false;
      size_t repeatcount=0;//// Counter tells how many times a vertex is repeting in the storage. I.e. whenever two vertices have the same number of tracks and all are equal.
      Vf=false;//initialize
      Otrk1.clear();//initialize

      for(size_t i=0; i<Vfound.size(); i++){if(Vfound.at(i)==x)Vf=true;}//omit comparing a vertex that was already merged

      if(tmpVertices.at(x).normalisedChiSquared()<5 && !Vf){// quality check plus already merged vertex check
        Otrk1=tmpVertices.at(x).originalTracks();// extract original tracks from first compared vertex
        for(size_t y=x+1; y<tmpVertices.size(); y++){/// iteration over all vertices
          //cout<<" 2nd it Flag "<<y<<endl;
          size_t similcount=0;//initialize
          Tfound.clear();//Restart Found tracks

          Vf=false;//Restart alread merged flag
          Otrk2.clear();// restart track container

          for(size_t i=0; i<Vfound.size(); i++){if(Vfound.at(i)==y)Vf=true;}//tell if the current vertex was already compared and merged

          if(tmpVertices.at(y).normalisedChiSquared()<5 && !Vf){ //quality check plus already merged check
            Otrktmp.clear();//restart Temporal track Storage
            Otrk2=tmpVertices.at(y).originalTracks();//Extract original tracks from second compared vertex
            for(size_t z=0; z<Otrk1.size(); z++){// first track set iterations
              for(size_t v=0; v<Otrk2.size(); v++){// second track set iterations
                if( deltaR(Otrk1.at(z).track().eta(),Otrk1.at(z).track().phi(),Otrk2.at(v).track().eta(),Otrk2.at(v).track().phi())==0
                   && Otrk1.at(z).track().pt()==Otrk2.at(v).track().pt() ){ //check whenever a track is shared in both sets
                     //cout<<"X "<<x<<" Y "<<y<<" Z "<<z<<" V "<<v<<endl;
                     //cout<<"  "<<deltaR(Otrk1.at(z).track().eta(),Otrk1.at(z).track().phi(),Otrk2.at(v).track().eta(),Otrk2.at(v).track().phi())
                     //<<"  "<<Otrk1.at(z).track().pt()<<"  "<<Otrk2.at(v).track().pt()<<endl;
                     similcount++;// count number of similar tracks from each set
                     rSimil=true;// FLag there are mergeable vertices
                     //breakSimil=true;
                     Otrktmp=Otrk1;// store all track from first set
                     Vfound.push_back(y);// store index of already compared and merged vertex
                     Tfound.push_back(v);// store index of already compared and similar track
                     Vsimil=true;//Flag the compared vertex is mergeable
                     //break;
                   }
              }
              //if(breakSimil) break;
            }
            if(similcount!=0){//check if two vertices share tracks
              for(size_t w=0; w<Otrk2.size(); w++){//iteration over second set of tracks
                Similf=false;//reset flag track is repeated
                for(size_t s=0; s<Tfound.size(); s++){if(w==Tfound.at(s)) Similf=true;}//checks flag if track is repeated
                if(!Similf) Otrktmp.push_back(Otrk2.at(w));// if track is not repeated store in temporal track storage
              }
            }
            //cout<<" Temporal track size: "<<Otrktmp.size()<<" 1st vertex iterator index: "<<x<<" 2nd vertex iterator index: "<< y <<endl;
            ////poner aqui el Fitter
            //cout<<"Simil Tracks: "<<similcount<<" Otrk1 "<< Otrk1.size()<<" Otrk2 "<<Otrk2.size()<<endl;
            //cout<<" Temporal track size: "<<Otrktmp.size()<<" 1st vertex iterator index: "<<x<<" 2nd vertex iterator index: "<< y <<endl;
            if(similcount==Otrk1.size() && Otrk1.size()==Otrk2.size())repeatcount++;// if two or more vertex are identical count
            /*if(similcount==Otrk1.size() && Otrk1.size()==Otrk2.size()){
              break;
            }*/
            if(Otrktmp.size()>2 && repeatcount<2){//check if temporal track could be fitted check if the identical vertex where already merged
              //cout<<"Simil Tracks: "<<similcount<<" Otrk1 "<< Otrk1.size()<<" Otrk2 "<<Otrk2.size()<<endl;
              //cout<<" Temporal track size: "<<Otrktmp.size()<<" 1st vertex iterator index: "<<x<<" 2nd vertex iterator index: "<< y <<endl;
              TransientVertex myVertex = fitter.vertex(Otrktmp);//fiit vertex
              if(myVertex.isValid()){// sanity check
                vertexcount++;// count number of vertices
                myVertices.push_back(myVertex);//restore merged vertices
              }
            }
          }
        }
        if(!Vsimil && tmpVertices.at(x).normalisedChiSquared()<5) finalVertices.push_back(tmpVertices.at(x));// if track is not mergeable and is good store as final vertex
        //cout<<"Simil vertices: "<<similcount<<endl;
      }
    }
    //cout<<" Valid vertex "<< vertexcount<<endl;
    //cout<<"\n\n\n Number of vertices after merging: "<<myVertices.size()<<" Final vertices: "<<finalVertices.size()<<"\n\n\n";
  }
  cout<<"Iterations: "<< its <<" Final vertices: "<<finalVertices.size()<<endl;

  for(size_t x=0; x<finalVertices.size(); x++){
    TransientVertex myV=finalVertices.at(x);
    Otrk1.clear();
    numEsecvec++;
    Esecvec_posx.push_back(myV.position().x());
    Esecvec_posy.push_back(myV.position().y());
    Esecvec_posz.push_back(myV.position().z());
    Esecvec_poserrorx.push_back(myV.positionError().cxx());
    Esecvec_poserrory.push_back(myV.positionError().cyy());
    Esecvec_poserrorz.push_back(myV.positionError().czz());
    Esecvec_chi2.push_back(myV.totalChiSquared());
    Esecvec_nodf.push_back(myV.degreesOfFreedom());
    Esecvec_normchi2.push_back(myV.normalisedChiSquared());
    Esecvec_eleTag.push_back(EleCount);
  }

  EleCount++;
}

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
