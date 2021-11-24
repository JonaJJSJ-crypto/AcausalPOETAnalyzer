// -*- C++ -*-
//
// Package:    TriggObjectAnalyzer
// Class:      TriggObjectAnalyzer
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

//classes to extract trigger information
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

//class to trigger result
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

//classes to save data
#include "TTree.h"
#include "TFile.h"
#include<vector>

//
// class declaration
//

class TriggObjectAnalyzer : public edm::EDAnalyzer {
   public:
      explicit TriggObjectAnalyzer(const edm::ParameterSet&);
      ~TriggObjectAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      
      //declare de filter (module) of the trigger
      std::string filterName_;
      HLTConfigProvider hltConfig_;      

      // ----------member data ---------------------------
      
      
      TTree *mtree;
      int numtrigobj; //number of trigger objects in the event
      std::vector<float> trigobj_e;
      std::vector<float> trigobj_pt;
      std::vector<float> trigobj_px;
      std::vector<float> trigobj_py;
      std::vector<float> trigobj_pz;
      std::vector<float> trigobj_eta;
      std::vector<float> trigobj_phi;
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

TriggObjectAnalyzer::TriggObjectAnalyzer(const edm::ParameterSet& iConfig)

{
//now do what ever initialization is needed
	filterName_ = iConfig.getParameter<std::string>("filterName");
	edm::Service<TFileService> fs;
	mtree = fs->make<TTree>("Events", "Events");
	
	
	mtree->Branch("numbertrigobj",&numtrigobj);
	mtree->Branch("trigobj_e",&trigobj_e);
	mtree->Branch("trigobj_pt",&trigobj_pt);
	mtree->Branch("trigobj_px",&trigobj_px);
	mtree->Branch("trigobj_py",&trigobj_py);
	mtree->Branch("trigobj_pz",&trigobj_pz);
	mtree->Branch("trigobj_eta",&trigobj_eta);
	mtree->Branch("trigobj_phi",&trigobj_phi);
	
}


TriggObjectAnalyzer::~TriggObjectAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void
TriggObjectAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   
   InputTag trigEventTag("hltTriggerSummaryAOD","","HLT"); //make sure have correct process on MC
   //data process=HLT, MC depends, Spring11 is REDIGI311X
   Handle<trigger::TriggerEvent> mytrigEvent;
   iEvent.getByLabel(trigEventTag,mytrigEvent);

   Handle<edm::TriggerResults> trigResults; 
   edm::InputTag trigResultsTag("TriggerResults","","HLT");
   iEvent.getByLabel(trigResultsTag,trigResults);

   const edm::TriggerNames& trigNames = iEvent.triggerNames(*trigResults);

   numtrigobj = 0;
   trigobj_e.clear();
   trigobj_pt.clear();
   trigobj_px.clear();
   trigobj_py.clear();
   trigobj_pz.clear();
   trigobj_eta.clear();
   trigobj_phi.clear();

   string toFind[4] = {"HLT_Photon36_R9Id85_Photon22_R9Id85_v","HLT_Photon36_R9Id85_Photon22_CaloId10_Iso50_v", "HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85_v","HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85_v"};

/*** Photon AOD files (used for the electron analysis part) may have any of the above triggers
 * and may also differ from file to file. In the loop below we check which one it is.
 * **/

//int triggerFound=9999;
string filterName = "none";
for(int j = 0; j < 4; j++){
for (unsigned int i = 0; i< trigNames.size(); i++)
{
	
	
	std::string trig = trigNames.triggerName(i);
	const unsigned int triggerIndex(hltConfig_.triggerIndex(trig));
	if ( trig.find(toFind[j]) !=std::string::npos ){
		
		int wr = trigResults->wasrun(triggerIndex);
		int acc = trigResults->accept(triggerIndex);
		int err = trigResults->error(triggerIndex);
		if(wr == 1 && acc == 1 && err == 0 ){
			const unsigned int m(hltConfig_.size(triggerIndex));
  			const vector<string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex));
  			const unsigned int moduleIndex(trigResults->index(triggerIndex));
			assert (moduleIndex<m);
			cout<<"modules for: "<<trig<<endl;
			for (unsigned int j=0; j<=moduleIndex; ++j) {
     				const string& moduleLabel(moduleLabels[j]);
				//cout<<'	'<<moduleLabel<<endl;	
				const unsigned int filterIndex(mytrigEvent->filterIndex(InputTag(moduleLabel,"","HLT")));
				if (filterIndex<mytrigEvent->sizeFilters() && j+1==moduleIndex){
					cout<<moduleLabel<<endl;
					filterName = moduleLabel;
					} 
				}
			//else cout<<"NoTrig "<<wr<<' '<<acc<< ' '<<err<<endl;
			//triggerFound=j;
			//i = trigNames.size();
			//j = 4;
			}
		//else cout<<"NoTrig "<<wr<<' '<<acc<< ' '<<err<<endl;
	}
	//else cout<<"NoTrig "<<wr<<' '<<acc<< ' '<<err<<endl;
		
	}
}
/*** Here we check if the trigger was activated and store the result
 * in the variable passTrig. This is a requirement for the analysis.
 * ***/
//cout<<triggerFound<<endl;
/*string filterName = "none";

if (triggerFound == 0)
{
	filterName = "hltEG22R9Id85DoubleLastFilterUnseeded";
}
else if (triggerFound == 1)
{
	filterName = "hltPhoton36R9Id85Photon22CaloId10Iso50EgammaDoubleLegCombLastFilter";
}
else if (triggerFound ==2) 
{
	filterName = "hltPhoton36CaloId10Iso50Photon22R9Id85EgammaDoubleLegCombLastFilter";
}
else
{
        filterName = "hltEG22CaloId10Iso50TrackIsoDoubleLastFilterUnseeded";
}*/

std::string e_filterName(filterName);

    trigger::size_type filterIndex = mytrigEvent->filterIndex(edm::InputTag(e_filterName,"",trigEventTag.process()));
    if(filterIndex<mytrigEvent->sizeFilters()){
    const trigger::Keys& trigKeys = mytrigEvent->filterKeys(filterIndex);
    const trigger::TriggerObjectCollection & trigObjColl(mytrigEvent->getObjects());

    //now loop of the trigger objects passing filter
    for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt){
      const trigger::TriggerObject trigobj = trigObjColl[*keyIt];

      //do what you want with the trigger objects, you have
      //eta,phi,pt,mass,p,px,py,pz,et,energy accessors
	    trigobj_e.push_back(trigobj.energy());
	    trigobj_pt.push_back(trigobj.pt());
	    trigobj_px.push_back(trigobj.px());
	    trigobj_py.push_back(trigobj.py());
	    trigobj_pz.push_back(trigobj.pz());
	    trigobj_eta.push_back(trigobj.eta());
	    trigobj_phi.push_back(trigobj.phi());

	    numtrigobj=numtrigobj+1;
    }
  }//end filter size check
  
  mtree->Fill();
  return;
  
}

// ------------ method called once each job just before starting event loop  ------------
void
TriggObjectAnalyzer::beginJob()
{}

// ------------ method called once each job just after ending the event loop  ------------
void
TriggObjectAnalyzer::endJob()
{}

// ------------ method called when starting to processes a run  ------------
void
TriggObjectAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
	bool changed(true);
    hltConfig_.init(iRun,iSetup,"HLT",changed);

}

// ------------ method called when ending the processing of a run  ------------
void
TriggObjectAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{}
// ------------ method called when starting to processes a luminosity block  ------------
void
TriggObjectAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
TriggObjectAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TriggObjectAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
 //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggObjectAnalyzer);
