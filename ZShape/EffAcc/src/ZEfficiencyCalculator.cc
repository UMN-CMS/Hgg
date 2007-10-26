#include "ZShape/EffAcc/interface/ZEfficiencyCalculator.h"


//
// constructors and destructor
//
ZEfficiencyCalculator::ZEfficiencyCalculator(const edm::ParameterSet& iConfig)
{
	
  //now do what ever initialization is needed

  outFileName_              = iConfig.getUntrackedParameter<std::string>("outHistogramsFile", "");
  writeHistoConservatively_ = iConfig.getUntrackedParameter<bool>("writeHistoBeforeEndJob", false);

  inputEffFileName_         = iConfig.getUntrackedParameter<std::string>("inputEffHistoFile", "");

  inputEffFileNameCut1_    = iConfig.getUntrackedParameter<std::string>("inputRootFile1", "");
  inputEffFileNameCut2_    = iConfig.getUntrackedParameter<std::string>("inputRootFile2", "");
  inputEffFileNameCut3_    = iConfig.getUntrackedParameter<std::string>("inputRootFile3", "");

  cut1HistoName_            = iConfig.getUntrackedParameter<std::string>("cut1Histo", "seeding alone");
  cut2HistoName_            = iConfig.getUntrackedParameter<std::string>("cut2Histo", "tracking alone");
  cut3HistoName_            = iConfig.getUntrackedParameter<std::string>("cut3Histo", "eReco alone");

  myElectrons_.clear();

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

  Handle< HepMCProduct > EvtHandle ;
   
  // find initial (unsmeared, unfiltered,...) HepMCProduct
  iEvent.getByLabel( "source", EvtHandle ) ;
  
  const HepMC::GenEvent* Evt = EvtHandle->GetEvent() ;
  
  int status;
  int pid;

  ////////////////////////
  // looking for electrons
  int myPid = 11;
  myElectrons_.clear();    myElectrons_.reserve(2);

  //
  // looping on all the particles of the simulated event
  for(HepMC::GenEvent::particle_const_iterator mcpart = Evt->particles_begin();
      mcpart != Evt->particles_end();
      ++ mcpart ) {

    status = (*mcpart)->status();
    pid    = (*mcpart)->pdg_id();
     
    if (status == 3 && abs(pid) == myPid ){

      const HepMC::GenVertex * vertex_=(*mcpart)->production_vertex();

      XYZTLorentzVector  vertex;
      vertex.SetPx( vertex_->position().x() );        vertex.SetPy( vertex_->position().y() );
      vertex.SetPz( vertex_->position().z() );        vertex.SetE(  0 );
       
      XYZTLorentzVector  momentum;
      momentum.SetPx(  (*mcpart)->momentum().x() );        momentum.SetPy(  (*mcpart)->momentum().y() );
      momentum.SetPz(  (*mcpart)->momentum().z() );        momentum.SetE(   (*mcpart)->momentum().t() );

     
      myElectrons_.push_back( RawParticle(momentum, vertex) ) ;

    }

  }
  // end loop on particles


  //
  // check than we have 2 and only 2 electrons
 
  if (myElectrons_.size() != 2)
    {
      edm::LogWarning("ZEfficiencyCalculator") << " we need two electrons (pid: " << myPid << " ) while we have: " << myElectrons_.size();
      return;
    }
  else
    {
      edm::LogWarning("ZEfficiencyCalculator") << " there are two electrons (pid: " << myPid << " )";
    }



  //
  // reorder electrons in pT
  float pt1 = myElectrons_[0].et();
  float pt2 = myElectrons_[1].et();

  // flipping order only if necessary
  if (pt2 > pt1){
    RawParticle A = myElectrons_[0];
    RawParticle B = myElectrons_[1];
    myElectrons_.clear();
    myElectrons_.push_back(B);
    myElectrons_.push_back(A);
  }
  
  



  //------------------------------------------------------------------//
  // all what is before this point goes to a separate producer
  //                           ++ decouple ++
  // all what is after this point goes to a separate analyzer
  //------------------------------------------------------------------//




  //
  // fill histograms for before any selection
  allCase_.Fill(myElectrons_[0], myElectrons_[1]);

  // apply acceptance cuts
  if ( !theEtaSelector_.isInAcceptance( myElectrons_[0])  ) return;  
  if ( !theEtaSelector_.isInAcceptance( myElectrons_[1])  ) return;

  // fill standard histograms after acceptance
  acceptance_.Fill(myElectrons_[0], myElectrons_[1]);

  //
  // apply efficiency step by step, histograms after each step!

  // seeding
  if(
     seedingCut_->passesCut (myElectrons_[0].eta() )    &&
     seedingCut_->passesCut (myElectrons_[1].eta() )
     )
    {//      edm::LogWarning("ZEfficiencyCalculator") << "stage 102: seeding cut passed ";
    } else
    {
      //      edm::LogWarning("ZEfficiencyCalculator") << "stage 103: seeding cut not passed ";
      return;
    }


  // fill histograms after cut on seeding
  aftercut1_.Fill(myElectrons_[0], myElectrons_[1]);
   
  //
  // tracking
  if(
     trackingCut_->passesCut (myElectrons_[0].eta() )    &&
     trackingCut_->passesCut (myElectrons_[1].eta() )
     )
    {//      edm::LogWarning("ZEfficiencyCalculator") << "stage 112: tracking cut passed ";
    } else
    {
      //      edm::LogWarning("ZEfficiencyCalculator") << "stage 113: tracking cut not passed ";
      return;
    }

  // fill histograms after cut on tracking
  aftercut2_.Fill(myElectrons_[0], myElectrons_[1]);

  //
  // eReco
  if(
     eRecoCut_->passesCut (myElectrons_[0].eta() )    &&
     eRecoCut_->passesCut (myElectrons_[1].eta() )
     )
    {//      edm::LogWarning("ZEfficiencyCalculator") << "stage 122: eReco cut passed ";
    } else
    {
      //      edm::LogWarning("ZEfficiencyCalculator") << "stage 123: eReco cut not passed ";
      return;
    }

  // fill histograms after cut on tracking
  aftercut3_.Fill(myElectrons_[0], myElectrons_[1]);

}





// ------------ method called once each job just before starting event loop  ------------
void 
ZEfficiencyCalculator::beginJob(const edm::EventSetup&)
{

  // file to write out the histograms produced by the ZEfficiencyCalculator
  histoFile_   = new TFile(outFileName_.c_str(),"RECREATE");


  // organize in directories in order to have histograms at each step
  TDirectory* td = histoFile_->mkdir("All");
  td->cd();
  allCase_.Book();

  td = histoFile_->mkdir("Accepted");
  td->cd();
  acceptance_.Book();

  td = histoFile_->mkdir("AfterCut1");
  td->cd();
  aftercut1_.Book();

  td = histoFile_->mkdir("AfterCut2");
  td->cd();
  aftercut2_.Book();

  td = histoFile_->mkdir("AfterCut3");
  td->cd();
  aftercut3_.Book();

  // if one wanted to be very safe, start writing to file histograms immediately
  if (writeHistoConservatively_) histoFile_->Write();



  // opening file where one gets efficiencies
  effHistoFile_ = new TFile(inputEffFileName_.c_str(),"READ");

  eff1File_ = new TFile(inputEffFileNameCut1_.c_str(),"READ");
  std::cout << "eff1File_ is: " << eff1File_ << std::endl;
  eff1_       = new EfficiencyStore(eff1File_);
  histoSeeding_ = eff1_->getValuesHisto1D();
  std::cout << "histoSeeding_  is: " << histoSeeding_ << std::endl;
  //   std::cout << "in ZEfficiencyCalculator found histo named: " << histoSeeding_->GetTitle()
  // 	    << " with " << histoSeeding_->GetNbinsX() << " bins" << std::endl;
  seedingCut_  = new   EfficiencyCut( histoSeeding_ );
  eff1File_->Close();


  eff2File_ = new TFile(inputEffFileNameCut2_.c_str(),"READ");
  eff2_       = new EfficiencyStore(eff2File_);
  histoTracking_ = eff2_->getValuesHisto1D();
  //  std::cout << "in ZEfficiencyCalculator found histo named: " << histoTracking_->GetTitle()
  //	    << " with " << histoTracking_->GetNbinsX() << " bins" << std::endl;
  trackingCut_  = new   EfficiencyCut( histoTracking_ );
  eff2File_->Close();


  eff3File_ = new TFile(inputEffFileNameCut3_.c_str(),"READ");
  eff3_       = new EfficiencyStore(eff3File_);
  histoEReco_ = eff3_->getValuesHisto1D();
  std::cout << "in ZEfficiencyCalculator found histo named: " <<histoEReco_->GetTitle()
	    << " with " << histoEReco_->GetNbinsX() << " bins" << std::endl;
  eRecoCut_  = new   EfficiencyCut( histoEReco_ );
  eff3File_->Close();

//    histoTracking_        = (TH1F *)  effHistoFile_ -> Get(cut2HistoName_.c_str());
//    trackingCut_ = new   EfficiencyCut( histoTracking_ );


//    histoEReco_   = (TH1F *)  effHistoFile_ -> Get(cut3HistoName_.c_str());
//    eRecoCut_     = new   EfficiencyCut(histoEReco_  );
}


// ------------ method called once each job just after ending the event loop  ------------
void 
ZEfficiencyCalculator::endJob() {

  effHistoFile_ -> Close();

  histoFile_->Write();  
  delete histoFile_;

}

//define this as a plug-in
DEFINE_FWK_MODULE(ZEfficiencyCalculator);
