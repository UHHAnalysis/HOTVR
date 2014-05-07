// $Id: MyCycle.cxx,v 1.10 2012/12/07 14:21:51 peiffer Exp $

#include <iostream>


using namespace std;

// Local include(s):

#include "include/MistagCycle.h"
#include "include/TopFitCalc.h"
#include "../SFrameTools/include/TMVA_tagger.h"
#include <TMVA/Reader.h>



ClassImp( MistagCycle );

MistagCycle::MistagCycle()
   : AnalysisCycle() {

  // constructor, declare additional variables that should be 
  // obtained from the steering-xml file
  
  // set the integrated luminosity per bin for the lumi-yield control plots
  
  
  /*    tmva_tagger= new TMVA_tagger();
  tmva_tagger->Set_Reader("Qref_weight");
  std::cout<<"BLAAAAAAA"<<std::endl;*/
  // DeclareProperty( "ApplyMttbarGenCut", m_mttgencut );
  SetIntLumiPerBin(500.);

}

MistagCycle::~MistagCycle() 
{
  // destructor
}

void MistagCycle::BeginCycle() throw( SError ) 
{
  // Start of the job, general set-up and definition of 
  // objects are done here

  // Important: first call BeginCycle of base class
  AnalysisCycle::BeginCycle();

  return;

}

void MistagCycle::EndCycle() throw( SError ) 
{
  // clean-up, info messages and final calculations after the analysis

  
  // call the base cycle class for all standard methods
  AnalysisCycle::EndCycle();

  return;

}

void MistagCycle::BeginInputData( const SInputData& id ) throw( SError ) 
{
  
  // declaration of histograms and selections

  // Important: first call BeginInputData of base class
  AnalysisCycle::BeginInputData( id );

  
   Selection* postselection = new Selection("postselection");
   /*  postselection->addSelectionModule(new TriggerSelection("HLT_HT750_v"));
 postselection->addSelectionModule(new HThadCut(0,int_infinity(),1000,int_infinity()));
 //  postselection->addSelectionModule(new NTopJetSelection(1,int_infinity(),300,int_infinity()));
  postselection->addSelectionModule(new NTopJetSelection(2,int_infinity(),400,int_infinity()));
  //cuts from paper 
  // postselection->addSelectionModule(new NJetSelection(2,int_infinity(),400));
  //   postselection->addSelectionModule(new NJetdeltaySelection(1.0));
  postselection->addSelectionModule(new NTopJetdeltaphiSelection(2.1));*/
 RegisterSelection(postselection);

  // -------------------- set up the selections ---------------------------

  /* Selection* BSel = new Selection( "BSelection");
  BSel->addSelectionModule(new NBTagSelection(1)); //at least one b tag

  Selection* NoBSel = new Selection( "NoBSelection");
  NoBSel->addSelectionModule(new NBTagSelection(0,0)); //no b tags

  Selection* chi2_selection= new Selection("chi2_selection");
  static Chi2Discriminator* m_chi2discr = new Chi2Discriminator();
  chi2_selection->addSelectionModule(new HypothesisDiscriminatorCut( m_chi2discr, -1*double_infinity(), 10));
  chi2_selection->addSelectionModule(new MttbarGenCut(0,700));

  Selection* TopSel = new Selection("TopSelection");
  //DO NOT use trigger selection in PROOF mode at the moment
  //TopSel->addSelectionModule(new TriggerSelection("HLT_PFJet320_v"));
  TopSel->addSelectionModule(new NTopJetSelection(1,int_infinity(),350,2.5));
  TopSel->addSelectionModule(new NTopTagSelection(1,int_infinity()));




  RegisterSelection(BSel);
  RegisterSelection(NoBSel);
  RegisterSelection(TopSel);
  RegisterSelection(chi2_selection);*/

  // ---------------- set up the histogram collections --------------------

  // histograms without any cuts
  // RegisterHistCollection( new MyHists("NoCuts") );

  //  static Chi2Discriminator* m_chi2discr = new Chi2Discriminator();
   /* RegisterHistCollection( new TopTagHists("test3"));
   RegisterHistCollection( new TopTagcontrol("preSelection"));
  RegisterHistCollection( new TopTagcontrol("postSelection"));
  RegisterHistCollection( new TopTagcontrol("postbSelection"));
  RegisterHistCollection( new TopTagcontrol("postNOBSelection"));
  RegisterHistCollection( new HypothesisHists("Chi2_BTag", m_chi2discr ) );*/


  //RegisterHistCollection( new HypothesisHists("Chi2_Presel", m_chi2discr ) );
    static Chi2Discriminator* m_chi2discr = new Chi2Discriminator();
    RegisterHistCollection( new TopTagHists("test3"));
    RegisterHistCollection( new SubstructureHists("substructure_alljets"));
    RegisterHistCollection( new SubstructureHists("substructure_cms"));
     RegisterHistCollection( new SubstructureHists("substructure_hep"));
     RegisterHistCollection( new SubstructureHists("substructure_mva"));
     RegisterHistCollection( new SubstructureHists("substructure_tobias"));
 RegisterHistCollection( new SubstructureHists("substructure_pt700"));
    RegisterHistCollection( new Mistagcontrol("preSelection"));
    RegisterHistCollection( new Mistagcontrol("postSelection"));
    RegisterHistCollection( new HypothesisHists("Chi2_BTag", m_chi2discr ) );
    RegisterHistCollection( new HypothesisHists("Chi2_Presel", m_chi2discr ) );
    RegisterHistCollection( new EventHists("Event_Presel") );
    RegisterHistCollection( new JetHists("Jets_Presel") );
    RegisterHistCollection( new ElectronHists("Electron_Presel") );
    RegisterHistCollection( new MuonHists("Muon_Presel") );
    RegisterHistCollection( new TauHists("Tau_Presel") );
    RegisterHistCollection( new TopJetHists("TopJets_Presel") );
    RegisterHistCollection( new EventHists("Event_Postsel") );
    RegisterHistCollection( new JetHists("Jets_Postsel") );
    RegisterHistCollection( new ElectronHists("Electron_Postsel") );
    RegisterHistCollection( new MuonHists("Muon_Postsel") );
    RegisterHistCollection( new TauHists("Tau_Postsel") );
    RegisterHistCollection( new TopJetHists("TopJets_Postsel") );
  /*  RegisterHistCollection( new HypothesisHists("Chi2_NoCuts", m_chi2discr ) );

  //histograms with and without b tagging
  RegisterHistCollection( new MyHists("BTag") );
  RegisterHistCollection( new MyHists("NoBTag") );
  RegisterHistCollection( new HypothesisHists("Chi2_BTag", m_chi2discr ) );
  RegisterHistCollection( new HypothesisHists("Chi2_NoBTag", m_chi2discr ) );

  // histograms after the top selection
  RegisterHistCollection( new MyHists("TopSel") );
  RegisterHistCollection( new HypothesisHists("Chi2_TopSel", m_chi2discr ) );*/

    

  InitHistos();
  tmva_tagger=new TMVA_tagger();
    
  return;

}

void MistagCycle::EndInputData( const SInputData& id ) throw( SError ) 
{
  AnalysisCycle::EndInputData( id );
  return;

}

void MistagCycle::BeginInputFile( const SInputData& id ) throw( SError ) 
{
  // Connect all variables from the Ntuple file with the ones needed for the analysis
  // The variables are commonly stored in the BaseCycleContaincer

  // important: call to base function to connect all variables to Ntuples from the input tree
  AnalysisCycle::BeginInputFile( id );
  
  return;

}

void MistagCycle::ExecuteEvent( const SInputData& id, Double_t weight) throw( SError ) 
{
  AnalysisCycle::ExecuteEvent( id, weight );
  //  tmva_tagger=TMVA_tagger::Instance();
   // this is the most important part: here the full analysis happens
  // user should implement selections, filling of histograms and results
  /*   BaseHists* HistsTopTag = GetHistCollection("test3");
   BaseHists* HistspreSelection = GetHistCollection("preSelection");
 BaseHists* HistspostSelection = GetHistCollection("postSelection");
 BaseHists* HistspostbSelection = GetHistCollection("postbSelection");
 BaseHists* HistspostNOBSelection = GetHistCollection("postNOBSelection");*/

//standard control plots
   BaseHists* HistsTopTag = GetHistCollection("test3");
   BaseHists* HistsSubstructure_alljets= GetHistCollection("substructure_alljets");
   BaseHists* HistsSubstructure_cms= GetHistCollection("substructure_cms");
   BaseHists* HistsSubstructure_hep= GetHistCollection("substructure_hep");
   BaseHists* HistsSubstructure_mva= GetHistCollection("substructure_mva");
    BaseHists* HistsSubstructure_tobias= GetHistCollection("substructure_tobias");
     BaseHists* HistsSubstructure_pt700= GetHistCollection("substructure_pt700");
   BaseHists* HistspreSelection = GetHistCollection("preSelection");
   BaseHists* HistspostSelection = GetHistCollection("postSelection");
   BaseHists* eventhists = GetHistCollection((std::string)("Event_Presel"));
    BaseHists* jethists = GetHistCollection((std::string)("Jets_Presel"));
    BaseHists* elehists = GetHistCollection((std::string)("Electron_Presel"));
    BaseHists* muonhists = GetHistCollection((std::string)("Muon_Presel"));
    BaseHists* tauhists = GetHistCollection((std::string)("Tau_Presel"));
    BaseHists* topjethists = GetHistCollection((std::string)("TopJets_Presel"));
    BaseHists* eventhists_post = GetHistCollection((std::string)("Event_Postsel"));
    BaseHists* jethists_post = GetHistCollection((std::string)("Jets_Postsel"));
    BaseHists* elehists_post = GetHistCollection((std::string)("Electron_Postsel"));
    BaseHists* muonhists_post = GetHistCollection((std::string)("Muon_Postsel"));
    BaseHists* tauhists_post = GetHistCollection((std::string)("Tau_Postsel"));
    BaseHists* topjethists_post = GetHistCollection((std::string)("TopJets_Postsel"));
     eventhists->Fill();
    jethists->Fill();
    elehists->Fill();
    muonhists->Fill();
    tauhists->Fill();
    topjethists->Fill();
    HistspreSelection->Fill();


     static Selection* Spostselection = GetSelection("postselection");
     //  if(!Spostselection->passSelection())  throw SError( SError::SkipEvent );

    Cleaner cleaner;
    EventCalc* calc = EventCalc::Instance();
    BaseCycleContainer* bcc = calc->GetBaseCycleContainer();

    std::vector<Jet> uncleaned_jets;
    for(unsigned int i=0; i<bcc->jets->size(); ++i) {
        uncleaned_jets.push_back(bcc->jets->at(i));
    }
    if(bcc->jets) cleaner.JetCleaner(100,2.1);

    //selection

TopJet antitaggedJet;
  bool antitag=false;
  for(unsigned int j=0; j<bcc->topjets->size();++j){
    TopJet topjet =bcc->topjets->at(j);
 
    double mmin3=0;
    double mjet3=0;
    int nsubjets3=0;
    TopTag(topjet,mjet3,nsubjets3,mmin3);
    if(mmin3<30 && mjet3>140 &&mjet3<250 && !antitag) {
      antitag=true;
      antitaggedJet=topjet;
    }
  } 
 
 if(antitag) for(unsigned int j=0; j<bcc->topjets->size();++j){
      TopJet topjet =bcc->topjets->at(j);
      if(topjet.v4()!=antitaggedJet.v4()){
     eventhists_post->Fill();
    jethists_post->Fill();
    elehists_post->Fill();
    muonhists_post->Fill();
    tauhists_post->Fill();
    topjethists_post->Fill();
    HistspostSelection->Fill();
     double mmin=0;
    double mjet=0;
    int nsubjets=0;
    double mva_value=0;
    if(topjet.pt()>700) HistsSubstructure_pt700->Fill2(topjet,mva_value);
    if(tmva_tagger->IsTagged("NPVweight",topjet,0.840,mva_value)) HistsSubstructure_mva->Fill2(topjet,mva_value);
     if(TopTag(topjet,mjet,nsubjets,mmin)) HistsSubstructure_cms->Fill2(topjet,mva_value);
     if(HepTopTagFull(topjet,calc->GetPFParticles())) HistsSubstructure_hep->Fill2(topjet,mva_value);
     if(tmva_tagger->IsTobiasTagged(topjet)) HistsSubstructure_tobias->Fill2(topjet,mva_value);
    HistsSubstructure_alljets->Fill2(topjet,mva_value);
    antitag=false;
     }
 }
  //Chi2Discriminator* tagchi2discr;
  // first step: call Execute event of base class to perform basic consistency checks
  // also, the good-run selection is performed there and the calculator is reset
   
    
      
	   	




      //     HistsTopTag->Fill();




    // if(bcc->taus) cleaner.TauCleaner(20,2.1);
    //   WriteOutputTree();
  return;
  
}


