#include "TProfile.h"
#include "ZShape/Base/interface/EfficiencyStatistics.h"
#include "ZShape/EffAcc/interface/ZEfficiencyCalculator.h"


//
// constructors and destructor
//
ZEfficiencyCalculator::ZEfficiencyCalculator(const edm::ParameterSet& iConfig) :
  m_srcTag(iConfig.getUntrackedParameter<edm::InputTag>("zrc",edm::InputTag("source")))
{
	
  //now do what ever initialization is needed

  outFileName_              = iConfig.getUntrackedParameter<std::string>("outHistogramsFile", "");
  writeHistoConservatively_ = iConfig.getUntrackedParameter<bool>("writeHistoBeforeEndJob", false);

  std::vector<edm::ParameterSet> zdefSetups =iConfig.getUntrackedParameter< std::vector<edm::ParameterSet> >("ZDefs");

  for (std::vector<edm::ParameterSet>::iterator i=zdefSetups.begin(); i!=zdefSetups.end(); ++i) {
    std::string name=i->getUntrackedParameter<std::string>("name");
    std::vector<std::string> req1=i->getUntrackedParameter<std::vector<std::string> >("e1");
    std::vector<std::string> req2=i->getUntrackedParameter<std::vector<std::string> >("e2");
    ZShapeZDef* zdef=new ZShapeZDef();
    for (std::vector<std::string>::iterator k=req1.begin(); k!=req1.end(); k++) 
      zdef->addCriterion(ZShapeZDef::crit_E1,*k);
    for (std::vector<std::string>::iterator k=req2.begin(); k!=req2.end(); k++) 
      zdef->addCriterion(ZShapeZDef::crit_E2,*k);
    zdefs_[name]=zdef;
  }

  std::vector<edm::ParameterSet> effs =iConfig.getUntrackedParameter< std::vector<edm::ParameterSet> >("Effs");

  for (std::vector<edm::ParameterSet>::iterator i=effs.begin(); i!=effs.end(); ++i) {
    std::string name=i->getUntrackedParameter<std::string>("name");
    std::string file=i->getUntrackedParameter<std::string>("effFile");
    //std::string var=i->getUntrackedParameter<std::string>("variable");
    loadEfficiency(name,file);
  }

  statsBox_.trials=iConfig.getUntrackedParameter<int>("es_trials",-1);
  if (statsBox_.trials>0) {
    statsBox_.targetEff=iConfig.getUntrackedParameter<std::string>("es_eff");
    statsBox_.targetZDef=iConfig.getUntrackedParameter<std::string>("es_zdef");

    if (efficiencies_.find(statsBox_.targetEff)==efficiencies_.end()) {
      edm::LogError("ZShape") << "Unable to find efficiency '" <<statsBox_.targetEff << "' to vary!";
      statsBox_.trials=0;
    }
    if (zdefs_.find(statsBox_.targetZDef)==zdefs_.end()) {
      edm::LogError("ZShape") << "Unable to find Z-definition '" <<statsBox_.targetZDef << "' to vary!";
      statsBox_.trials=0;
    }

  }

  

}


ZEfficiencyCalculator::~ZEfficiencyCalculator()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void ZEfficiencyCalculator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace cms;
  int pass=0;

  Handle< HepMCProduct > EvtHandle ;
   
  // find initial (unsmeared, unfiltered,...) HepMCProduct
  iEvent.getByLabel(m_srcTag, EvtHandle ) ;
  
  const HepMC::GenEvent* Evt = EvtHandle->GetEvent() ;
  fillEvent(Evt);

  if (evt_.n_elec!=2) return; // need 2 and only 2

  //
  // fill histograms for before any selection
  allCase_.Fill(evt_.elec_[0].p4_, evt_.elec_[1].p4_);
  
  // these stages merely _fill_ bits.  They do not apply cuts!
  acceptanceCuts();
  
  EfficiencyCut* keep=0;

  do {

    if (statsBox_.trials>=1) {
      // swap to alternate universe
      if (pass==0) keep=theCuts_[statsBox_.targetEff];
      theCuts_[statsBox_.targetEff]=statsBox_.alternateCuts[pass];
    }

    for (std::map<std::string,EfficiencyCut*>::iterator cut=theCuts_.begin(); cut!=theCuts_.end(); cut++) {
      for (int ne=0; ne<2; ne++) {
	evt_.elec_[ne].cutResult(cut->first,cut->second->passesCut(evt_.elec_[ne].p4_));
      }
    }

    //------------------------------------------------------------------//
    // all what is before this point goes to a separate producer
    //                           ++ decouple ++
    // all what is after this point goes to a separate analyzer
    //------------------------------------------------------------------//
    
    bool pairing;
    
    if (pass==0) {
      evt_.dump();

      for (std::map<std::string,ZShapeZDef*>::const_iterator q=zdefs_.begin(); q!=zdefs_.end(); ++q) {
	ZPlots* plots=zplots_[q->first];
	// acceptance is always the first cut
	if (!q->second->pass(evt_,1,1,0,&pairing)) continue;
	
	// fill standard histograms after acceptance
	if (!pairing) plots->acceptance_.Fill(evt_.elec_[0].p4_, evt_.elec_[1].p4_);
	else plots->acceptance_.Fill(evt_.elec_[1].p4_, evt_.elec_[0].p4_);
	
	// next n-cuts
	bool ok=true;
	for (int j=1; ok && j<q->second->criteriaCount(ZShapeZDef::crit_E1); j++) {
	  ok=q->second->pass(evt_,j+1,j+1,0,&pairing);
	  if (ok)
	    if (!pairing) 
	      plots->postCut_[j-1].Fill(evt_.elec_[0].p4_, evt_.elec_[1].p4_);
	    else
	      plots->postCut_[j-1].Fill(evt_.elec_[1].p4_, evt_.elec_[0].p4_);
	}
      }
    }
    if (statsBox_.trials>0) {
      ZShapeZDef* zdef=zdefs_[statsBox_.targetZDef];
      if (zdef->pass(evt_,zdef->criteriaCount(ZShapeZDef::crit_E1),zdef->criteriaCount(ZShapeZDef::crit_E2),0,&pairing)) {
	if (!pairing) 
	  statsBox_.hists[pass].Fill(evt_.elec_[0].p4_, evt_.elec_[1].p4_);
	else
	  statsBox_.hists[pass].Fill(evt_.elec_[1].p4_, evt_.elec_[0].p4_);
      }
    }
    pass++;
  } while (statsBox_.trials>0 && pass<statsBox_.trials);

  if (keep!=0) theCuts_[statsBox_.targetEff]=keep;
}


void ZEfficiencyCalculator::fillEvent(const HepMC::GenEvent* Evt) {
  int status;
  int pid;

  ////////////////////////
  // looking for electrons
  int myPid = 11;
  evt_.clear();
  int ne=0;

  //
  // looping on all the particles of the simulated event
  for(HepMC::GenEvent::particle_const_iterator mcpart = Evt->particles_begin();
      mcpart != Evt->particles_end();
      ++ mcpart ) {

    status = (*mcpart)->status();
    pid    = (*mcpart)->pdg_id();
     
    if (status == 3 && abs(pid) == myPid ){

      const HepMC::GenVertex * vertex_=(*mcpart)->production_vertex();

      if (ne==0) {
	evt_.z0_=::math::XYZVector(vertex_->position().x(),vertex_->position().y(),vertex_->position().z() );
      }
      
      math::XYZTLorentzVector  momentum;
      momentum.SetPx(  (*mcpart)->momentum().x() );        momentum.SetPy(  (*mcpart)->momentum().y() );
      momentum.SetPz(  (*mcpart)->momentum().z() );        momentum.SetE(   (*mcpart)->momentum().t() );
      
      if (ne<2) {
	evt_.elec_[ne].p4_=momentum;
      }     
      ne++;

    }

  }
  // end loop on particles
  evt_.n_elec=ne;


  //
  // check than we have 2 and only 2 electrons
 
  if (ne!= 2)
    {
      edm::LogWarning("ZEfficiencyCalculator") << " we need two electrons (pid: " << myPid << " ) while we have: " << ne;
      return;
    }


  //
  // reorder electrons in pT
  float pt1 = evt_.elec_[0].p4_.Pt();
  float pt2 = evt_.elec_[1].p4_.Pt();

  // flipping order only if necessary
  if (pt2 > pt1){
    std::swap( evt_.elec_[0], evt_.elec_[1]);
  }
}  
  





// ------------ method called once each job just before starting event loop  ------------
void 
ZEfficiencyCalculator::beginJob(const edm::EventSetup&)
{

  if (statsBox_.trials>0) createAlternateEfficiencies();

  // file to write out the histograms produced by the ZEfficiencyCalculator
  histoFile_   = new TFile(outFileName_.c_str(),"RECREATE");

  TDirectory* td=histoFile_->mkdir("All");
  td->cd();
  allCase_.Book();

  for (std::map<std::string,ZShapeZDef*>::const_iterator q=zdefs_.begin(); q!=zdefs_.end(); ++q) {
    ZPlots* zplots=new ZPlots(q->second->criteriaCount(ZShapeZDef::crit_E1)-1); // -1 for acceptance
    
    TDirectory* pd = histoFile_->mkdir(q->first.c_str());

    // organize in directories in order to have histograms at each step

    td = pd->mkdir("Acceptance");
    td->cd();
    zplots->acceptance_.Book();

    for (int i=1; i<q->second->criteriaCount(ZShapeZDef::crit_E1); ++i) {
      char dirname[1024];
      std::string c1=q->second->criteria(ZShapeZDef::crit_E1,i);
      //Change this, it assumes there must always be the same number of criteria for both electrons
      //  one could add in a dummy criteria
      std::string c2=q->second->criteria(ZShapeZDef::crit_E2,i);
      if (c1==c2) 
	sprintf(dirname,"C%02d-%s",i,c1.c_str());
      else 
	sprintf(dirname,"C%02d-%s-%s",i,c1.c_str(),c2.c_str());
      td = pd->mkdir(dirname);
      td->cd();
      zplots->postCut_[i-1].Book();
    }

    zplots_[q->first]=zplots;
  }
 
  // book trials histograms
  if (statsBox_.trials>0) {
    char dirname[1024];
    sprintf(dirname,"EffStats_%s_%s",statsBox_.targetEff.c_str(),statsBox_.targetZDef.c_str());     
    TDirectory* pd = histoFile_->mkdir(dirname);
    statsBox_.hists=std::vector<EffHistos>(statsBox_.trials);
    for (int j=0; j<statsBox_.trials;j++) {
     sprintf(dirname,"Trial%d",j+1);
     td = pd->mkdir(dirname);
     td->cd();
     statsBox_.hists[j].Book();
     statsBox_.rawHists[j]->SetDirectory(td);
    }
    statsBox_.cutsProfile->SetDirectory(pd);
  }

  // if one wanted to be very safe, start writing to file histograms immediately
  if (writeHistoConservatively_) histoFile_->Write();

}


void ZEfficiencyCalculator::loadEfficiency(const std::string& name, const std::string& fname) {

  if (fname.find(".root")!=std::string::npos) {
    TFile effFile(fname.c_str(),"READ");
    
    edm::LogInfo("ZShape") << "Reading ROOT " << fname << " for efficiency " << name;
    EfficiencyStore* effs       = new EfficiencyStore(&effFile, name);
    TH1F* histo = effs->getValuesHisto1D();
    EfficiencyCut* theCut  = new   EfficiencyCut(histo);
    efficiencies_[name]=effs;
    theCuts_[name]=theCut;
  } else {
    edm::LogInfo("ZShape") << "Reading TEXT " << fname << " for efficiency " << name;
    EfficiencyStore* effs       = new EfficiencyStore(fname);
    TH1F* histo = effs->getValuesHisto1D();
    EfficiencyCut* theCut  = new   EfficiencyCut(histo);
    efficiencies_[name]=effs;
    theCuts_[name]=theCut;
  }
}


// ------------ method called once each job just after ending the event loop  ------------
void 
ZEfficiencyCalculator::endJob() {

  histoFile_->Write();  
  delete histoFile_;

}

void ZEfficiencyCalculator::createAlternateEfficiencies() {
  if (statsBox_.trials<1) return;
  EfficiencyStore* es=efficiencies_[statsBox_.targetEff];

  TH1* val=es->getValuesHisto1D();
  TH1* denom=es->getDenominatorHisto1D();

  for (int i=0; i<statsBox_.trials; ++i) {
    char name[1024];
    snprintf(name,1024,"StatEff_%s_%s_%d",statsBox_.targetEff.c_str(),statsBox_.targetZDef.c_str(),i+1);
    TH1F* randy=new TH1F(name,name,val->GetXaxis()->GetNbins(),val->GetXaxis()->GetXmin(),val->GetXaxis()->GetXmax());
    if (i==0) {
      snprintf(name,1024,"Variations_%s_%s",statsBox_.targetEff.c_str(),statsBox_.targetZDef.c_str());
      statsBox_.cutsProfile=new TProfile(name,name,val->GetXaxis()->GetNbins(),val->GetXaxis()->GetXmin(),val->GetXaxis()->GetXmax());
    }
    for (int j=1; j<=val->GetXaxis()->GetNbins(); j++) {
      double ave=val->GetBinContent(j);
      int den=int(denom->GetBinContent(j));
      double rv=0;
      if (den>0 && ave>0) { // HACK! (ave should be able to be zero) 
	EfficiencyStatistics esd(ave,den);
	esd.setRandom(&randomNum);
	rv=esd.rand();
      }
      randy->Fill(val->GetBinCenter(j),rv);
      statsBox_.cutsProfile->Fill(val->GetBinCenter(j),rv);
    }
    statsBox_.rawHists.push_back(randy);
    statsBox_.alternateCuts.push_back(new EfficiencyCut(randy));
  }
}

void ZEfficiencyCalculator::acceptanceCuts() {
  for (int i=0; i<2; i++) {
    evt_.elec_[i].cutResult("ACC(EB)",theEtaSelector_.isInAcceptance(evt_.elec_[i].p4_,EtaAcceptance::zone_EB));
    evt_.elec_[i].cutResult("ACC(EE)",theEtaSelector_.isInAcceptance(evt_.elec_[i].p4_,EtaAcceptance::zone_EE));
    evt_.elec_[i].cutResult("ACC(ECAL)",theEtaSelector_.isInAcceptance(evt_.elec_[i].p4_,EtaAcceptance::zone_ECAL));
    evt_.elec_[i].cutResult("ACC(HF)",theEtaSelector_.isInAcceptance(evt_.elec_[i].p4_,EtaAcceptance::zone_HF));
    evt_.elec_[i].cutResult("ACC(ANY)",theEtaSelector_.isInAcceptance(evt_.elec_[i].p4_,EtaAcceptance::zone_ANY));
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZEfficiencyCalculator);
