// $Id: RocCycle.cxx,v 1.10 2012/12/07 14:21:51 peiffer Exp $

#include <iostream>

using namespace std;

#include "include/JetDisplayCycle.h"
#include "include/JetDisplayHists.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/contrib/HOTVR.hh"
#include "SFrameAnalysis/include/SelectionModules.h"
//#include "SFrameAnalysis/include/HypothesisHists.h"
using namespace fastjet;
using namespace contrib;
ClassImp( JetDisplayCycle );

JetDisplayCycle::JetDisplayCycle()
   : AnalysisCycle() {

  // constructor, declare additional variables that should be 
  // obtained from the steering-xml file
  
  // set the integrated luminosity per bin for the lumi-yield control plots
  SetIntLumiPerBin(500.);

}

JetDisplayCycle::~JetDisplayCycle() 
{
  // destructor
}

void JetDisplayCycle::BeginCycle() throw( SError ) 
{
  // Start of the job, general set-up and definition of 
  // objects are done here

  // Important: first call BeginCycle of base class
  AnalysisCycle::BeginCycle();

  return;

}

void JetDisplayCycle::EndCycle() throw( SError ) 
{
  // clean-up, info messages and final calculations after the analysis

  
  // call the base cycle class for all standard methods
  AnalysisCycle::EndCycle();

  return;

}

void JetDisplayCycle::BeginInputData( const SInputData& id ) throw( SError ) 
{
  // declaration of histograms and selections.
  // AnalysisCyle expects Selections and HistCollections to be registered here.
  // Their memory will be released in AnalysisCycle::EndInputData.

  // Important: first call BeginInputData of base class
  AnalysisCycle::BeginInputData( id );

 

  // ---------------- set up the histogram collections --------------------
    for(unsigned int i=0;i<_event_max;i++){
    TString hname = TString::Format("JetDisplay_event%d",i);
    RegisterHistCollection( new JetDisplayHists(hname));
    hname = TString::Format("JetDisplay_akt_event%d",i);
    RegisterHistCollection( new JetDisplayHists(hname));
    hname = TString::Format("JetDisplay_ca_event%d",i);
    RegisterHistCollection( new JetDisplayHists(hname));
    hname = TString::Format("JetDisplay_highpT_event%d",i);
    RegisterHistCollection( new JetDisplayHists(hname));
    hname = TString::Format("JetDisplay_highpT_akt_event%d",i);
    RegisterHistCollection( new JetDisplayHists(hname));
    hname = TString::Format("JetDisplay_highpT_ca_event%d",i);
    RegisterHistCollection( new JetDisplayHists(hname));
    cout<<"Initializing "<<(double)i*100./(double)_event_max<<"%"<<std::endl;
    }
  m_counter=0;
  m_counter2=0;
  toptagger=new TopTagger;
  clustering = new Clustering(m_clustering);
  clustering_ca= new Clustering("/afs/desy.de/user/t/tlapsien/Analysis/Sframex/HOTVR/config/ca.config");
  clustering_akt= new Clustering("/afs/desy.de/user/t/tlapsien/Analysis/Sframex/HOTVR/config/akt.config");
  IR_Saftey= new Infrared_Saftey();
  // important: initialise histogram collections after their definition
  InitHistos();

}

void JetDisplayCycle::EndInputData( const SInputData& id ) throw( SError ) 
{
  AnalysisCycle::EndInputData( id );
}

void JetDisplayCycle::BeginInputFile( const SInputData& id ) throw( SError ) 
{
  // Connect all variables from the Ntuple file with the ones needed for the analysis
  // The variables are commonly stored in the BaseCycleContaincer
  // important: call to base function to connect all variables to Ntuples from the input tree
  AnalysisCycle::BeginInputFile( id );
}

void JetDisplayCycle::ExecuteEvent( const SInputData& id, Double_t weight) throw( SError ) 
{
  // this is the most important part: here the full analysis happens
  // user should implement selections, filling of histograms and results

  // first step: call Execute event of base class to perform basic consistency checks
  // also, the good-run selection is performed there and the calculator is reset
  AnalysisCycle::ExecuteEvent( id, weight );

  // get the histogram collections. NOTE: this could be done more performant by making
  // all thse BaseHists* vairables private member variables of RocCycle and
  // setting them in BeginInputData. Then, there is no need here to call GetHistColletion ...

  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();
  bool IsRealData = calc->IsRealData();

  // BaseHists* Jetdisplay_hists_event[100];
  // for(unsigned int i=0;i<5;i++) {
  if(m_counter>=_event_max && m_counter2>=_event_max) throw SError( SError::SkipEvent );
  //if(calc->GetEventNum()!=5215 && calc->GetEventNum()!=8633)  throw SError( SError::SkipEvent );
     std::cout<<calc->GetEventNum()<<std::endl;
  //std::cout<<"\r"<<"Finding events.. \n";



 
 
 
 
  matching= new Matching();
  std::vector<GenParticle>* genparticles = calc->GetGenParticles();
  matching->Run_matching(genparticles);
  std::vector<fastjet::PseudoJet> parts= matching->get_hadrons();
  std::vector<fastjet::PseudoJet> denominator_jets=matching->get_denominator_jets(id.GetVersion());

  if(denominator_jets.size()<2) throw SError( SError::SkipEvent );
  
  std::vector<fastjet::PseudoJet> hotvr_jets, soft_jets, rejected_jets,rejected_subjets, ca_jets, akt_jets;
  hotvr_jets=clustering->get_clustered_jets(parts);

  

 

  if(hotvr_jets.size()<2)  throw SError( SError::SkipEvent );
 
  if(hotvr_jets.at(0).user_info<HOTVRinfo>().nsubjets()!=3)   throw SError( SError::SkipEvent );

 if(hotvr_jets.at(1).user_info<HOTVRinfo>().nsubjets()!=3) throw SError( SError::SkipEvent );

  if(hotvr_jets.at(0).m()<140 || hotvr_jets.at(0).m()>220 ) throw SError( SError::SkipEvent );
 
  if(hotvr_jets.at(1).m()<140 || hotvr_jets.at(1).m()>220 ) throw SError( SError::SkipEvent );

   if(abs(hotvr_jets.at(0).phi_std())>2. ||  abs(hotvr_jets.at(1).phi_std())>2.)  throw SError( SError::SkipEvent );

  if(!(hotvr_jets.at(0).pt()>200 && hotvr_jets.at(0).pt()<450) && !(hotvr_jets.at(0).pt()>800 && hotvr_jets.at(0).pt()<1000))  throw SError( SError::SkipEvent );

 
  std::cout<<"----------------------------"<<std::endl;

   IR_Saftey->add_grid(parts,75,75);
  clustering->Reset();
  hotvr_jets=clustering->get_clustered_jets(parts);
  soft_jets=clustering->get_soft_jets();
  rejected_jets=clustering->get_rejected_jets();
  rejected_subjets=clustering->get_rejected_subjets();
  ca_jets=clustering_ca->get_clustered_jets(parts);
  akt_jets=clustering_akt->get_clustered_jets(parts);
  
 
 
   
    
    //matching_radius=1.5;
 
  if(hotvr_jets.at(0).pt()>200 && hotvr_jets.at(0).pt()<450)	  
	  if(m_counter<_event_max){  
	     string hname = "JetDisplay_event"+ std::to_string(m_counter);
	    string hname2 = "JetDisplay_akt_event"+ std::to_string(m_counter);
	    string hname3 = "JetDisplay_ca_event"+ std::to_string(m_counter);
	    m_Jetdisplay_hists_event[m_counter]= GetHistCollection(hname);
	    m_Jetdisplay_akt_hists_event[m_counter]= GetHistCollection(hname2);
	    m_Jetdisplay_ca_hists_event[m_counter]= GetHistCollection(hname3);
	    
	    ((JetDisplayHists*) m_Jetdisplay_hists_event[m_counter])->FillEvent(hotvr_jets,parts,soft_jets,rejected_jets, rejected_subjets);
	    ((JetDisplayHists*) m_Jetdisplay_akt_hists_event[m_counter])->FillEvent(akt_jets,parts,soft_jets,rejected_jets, rejected_subjets);
	    ((JetDisplayHists*) m_Jetdisplay_ca_hists_event[m_counter])->FillEvent(ca_jets,parts,soft_jets,rejected_jets, rejected_subjets);
	      m_counter++;



	  }

 if(hotvr_jets.at(0).pt()>800 && hotvr_jets.at(0).pt()<1000)
  if(m_counter2<_event_max){  
     string hname = "JetDisplay_highpT_event"+ std::to_string(m_counter2);
    string hname2 = "JetDisplay_highpT_akt_event"+ std::to_string(m_counter2);
    string hname3 = "JetDisplay_highpT_ca_event"+ std::to_string(m_counter2);
    m_Jetdisplay_hists_highpT_event[m_counter2]= GetHistCollection(hname);
    m_Jetdisplay_akt_hists_highpT_event[m_counter2]= GetHistCollection(hname2);
    m_Jetdisplay_ca_hists_highpT_event[m_counter2]= GetHistCollection(hname3);
    
    ((JetDisplayHists*) m_Jetdisplay_hists_highpT_event[m_counter2])->FillEvent(hotvr_jets,parts,soft_jets,rejected_jets, rejected_subjets);
    ((JetDisplayHists*) m_Jetdisplay_akt_hists_highpT_event[m_counter2])->FillEvent(akt_jets,parts,soft_jets,rejected_jets, rejected_subjets);
    ((JetDisplayHists*) m_Jetdisplay_ca_hists_highpT_event[m_counter2])->FillEvent(ca_jets,parts,soft_jets,rejected_jets, rejected_subjets);
       m_counter2++;
  }




  /* if(m_counter<50){  
    string hname = "JetDisplay_event"+ std::to_string(m_counter);
    string hname2 = "JetDisplay_akt_event"+ std::to_string(m_counter);
    string hname3 = "JetDisplay_ca_event"+ std::to_string(m_counter);
    m_Jetdisplay_hists_event[m_counter]= GetHistCollection(hname);
    m_Jetdisplay_akt_hists_event[m_counter]= GetHistCollection(hname2);
    m_Jetdisplay_ca_hists_event[m_counter]= GetHistCollection(hname3);
   
    if(((JetDisplayHists*) m_Jetdisplay_hists_event[m_counter])->FillEvent(hotvr_jets,parts,soft_jets,rejected_jets)){
      // ((JetDisplayHists*) m_Jetdisplay_akt_hists_event[m_counter])->FillEvent(200,450,"akt",0.8);
      //((JetDisplayHists*) m_Jetdisplay_ca_hists_event[m_counter])->FillEvent(200,450,"ca",0.8);
       m_counter++;
    }
  }

  if(m_counter2<50){  
    string hname = "JetDisplay_highpT_event"+ std::to_string(m_counter2);
    string hname2 = "JetDisplay_highpT_akt_event"+ std::to_string(m_counter2);
    string hname3 = "JetDisplay_highpT_ca_event"+ std::to_string(m_counter2);
    m_Jetdisplay_hists_highpT_event[m_counter2]= GetHistCollection(hname);
    m_Jetdisplay_akt_hists_highpT_event[m_counter2]= GetHistCollection(hname2);
    m_Jetdisplay_ca_hists_highpT_event[m_counter2]= GetHistCollection(hname3);
  
    if(((JetDisplayHists*) m_Jetdisplay_hists_highpT_event[m_counter2])->FillEvent(hotvr_jets,parts,soft_jets,rejected_jets)){
      //((JetDisplayHists*) m_Jetdisplay_akt_hists_highpT_event[m_counter2])->FillEvent(800,10000,"akt",0.8);
      //((JetDisplayHists*) m_Jetdisplay_ca_hists_highpT_event[m_counter2])->FillEvent(800,10000,"ca",0.8);
       m_counter2++;
    }
  }
  */

 if(m_counter2==22) std::cout<<"22 "<<calc->GetEventNum()<<std::endl;
 if(m_counter==37) std::cout<<"38 "<<calc->GetEventNum()<<std::endl;

 std::cout	  <<"\r"<<"Progress "<<(m_counter2+m_counter)/2./_event_max*100.<<"%"<<std::flush;
 
  clustering->Reset();
  clustering_ca->Reset();
  clustering_akt->Reset();
  
  return;
  
}


