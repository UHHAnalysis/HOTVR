// $Id: RocCycle.cxx,v 1.10 2012/12/07 14:21:51 peiffer Exp $

#include <iostream>

using namespace std;

#include "include/JetDisplayAnimationCycle.h"
#include "include/JetDisplayHists.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/contrib/HOTVR.hh"
#include "SFrameAnalysis/include/SelectionModules.h"
//#include "SFrameAnalysis/include/HypothesisHists.h"
using namespace fastjet;
using namespace contrib;
ClassImp( JetDisplayAnimationCycle );

JetDisplayAnimationCycle::JetDisplayAnimationCycle()
   : AnalysisCycle() {

  // constructor, declare additional variables that should be 
  // obtained from the steering-xml file
  
  // set the integrated luminosity per bin for the lumi-yield control plots
  SetIntLumiPerBin(500.);

}

JetDisplayAnimationCycle::~JetDisplayAnimationCycle() 
{
  // destructor
}

void JetDisplayAnimationCycle::BeginCycle() throw( SError ) 
{
  // Start of the job, general set-up and definition of 
  // objects are done here

  // Important: first call BeginCycle of base class
  AnalysisCycle::BeginCycle();

  return;

}

void JetDisplayAnimationCycle::EndCycle() throw( SError ) 
{
  // clean-up, info messages and final calculations after the analysis

  
  // call the base cycle class for all standard methods
  AnalysisCycle::EndCycle();

  return;

}

void JetDisplayAnimationCycle::BeginInputData( const SInputData& id ) throw( SError ) 
{
  // declaration of histograms and selections.
  // AnalysisCyle expects Selections and HistCollections to be registered here.
  // Their memory will be released in AnalysisCycle::EndInputData.

  // Important: first call BeginInputData of base class
  AnalysisCycle::BeginInputData( id );

 

  // ---------------- set up the histogram collections --------------------
  Book(TH2F("JetDisplay_decay","Jet event display",50,-PI,PI,50,-PI,PI));
  
  m_counter=0;
  m_counter2=0;
  toptagger=new TopTagger;
  // clustering = new Clustering(m_clustering);
  clustering_ca= new Clustering("/afs/desy.de/user/t/tlapsien/Analysis/Sframex/HOTVR/config/ca.config");
  clustering_akt= new Clustering("/afs/desy.de/user/t/tlapsien/Analysis/Sframex/HOTVR/config/akt.config");
   clustering = new Clustering(m_clustering);
  IR_Saftey= new Infrared_Saftey();
  // important: initialise histogram collections after their definition
  InitHistos();

}

void JetDisplayAnimationCycle::EndInputData( const SInputData& id ) throw( SError ) 
{
  AnalysisCycle::EndInputData( id );
}

void JetDisplayAnimationCycle::BeginInputFile( const SInputData& id ) throw( SError ) 
{
  // Connect all variables from the Ntuple file with the ones needed for the analysis
  // The variables are commonly stored in the BaseCycleContaincer
  // important: call to base function to connect all variables to Ntuples from the input tree
  AnalysisCycle::BeginInputFile( id );
}

bool JetDisplayAnimationCycle::Jet_contains_hard_particles(fastjet::PseudoJet jet,double ptcut){
  int counter=0;
  for(int i=0;i<jet.constituents().size();i++) if(jet.constituents().at(i).pt()>ptcut) counter++; 
  if(counter>1) return true;
  else return false;
}

void JetDisplayAnimationCycle::ExecuteEvent( const SInputData& id, Double_t weight) throw( SError ) 
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
 
  TTbarGen2* Decay = calc->GetTTbarGen2();
    
  std::vector<GenParticle> Had_Tops;
  std::vector<GenParticle> decay_products;
   bool isHadronic=false;
   if (Decay->IsTopHadronicDecay()) {
     Had_Tops.push_back(Decay->Top());
     decay_products.push_back(Decay->Wdecay1());
     decay_products.push_back(Decay->Wdecay2());
     decay_products.push_back(Decay->bTop());
   }
   if (Decay->IsAntiTopHadronicDecay()){ 
      Had_Tops.push_back(Decay->Antitop());
      decay_products.push_back(Decay->WMinusdecay1());
      decay_products.push_back(Decay->WMinusdecay2());    
      decay_products.push_back(Decay->bAntitop());
    }
  // BaseHists* Jetdisplay_hists_event[100];
  // for(unsigned int i=0;i<5;i++) {
  if(m_counter>=_event_max && m_counter2>=_event_max) throw SError( SError::SkipEvent );
  
  
  //std::cout<<"\r"<<"Finding events.. \n";



 
 
 
 
  matching= new Matching();
  std::vector<GenParticle>* genparticles = calc->GetGenParticles();
  matching->Run_matching(genparticles);
  std::vector<fastjet::PseudoJet> parts= matching->get_hadrons();
  std::vector<fastjet::PseudoJet> denominator_jets=matching->get_denominator_jets(id.GetVersion());

  if(denominator_jets.size()<2) throw SError( SError::SkipEvent );
  
  std::vector<fastjet::PseudoJet> hotvr_jets, soft_jets, rejected_jets, ca_jets, akt_jets;
  
    IR_Saftey->add_grid(parts,75,75);
   hotvr_jets=clustering->get_clustered_jets(parts);
   /// VariableRPlugin pluginAKT(rho, min_r, max_r, VariableRPlugin::AKTLIKE);
  //fastjet::JetDefinition jet_defCA(&pluginAKT);


  // VariableRPlugin plugin(600, 0.1, 1.5,VariableRPlugin::CALIKE );
  //fastjet::JetDefinition jet_def(&plugin);
  //fastjet::ClusterSequence cs2(parts, jet_def);
  // hotvr_jets = sorted_by_pt(cs2.inclusive_jets(100));

 fastjet::ClusterSequence* cs=clustering->get_clustsq();


  std::vector<int> history=cs->unique_history_order();
 std::cout<<"!!!!"<<std::endl;
  std::vector<fastjet::PseudoJet> jets=cs->jets();
  std::vector<fastjet::PseudoJet> clustered_jets;
  std::vector<fastjet::PseudoJet> childless=cs->childless_pseudojets();

  int nhistos=0;
  
  std::cout<<"childless "<<childless.size()<<std::endl;

  std::vector<fastjet::PseudoJet> jets1=cs->inclusive_jets(100);

   std::cout<<history.size()<<std::endl;
  std::cout<<cs->jets().size()<<std::endl;
  for(int i=0;i<jets.size();i++) if(cs->object_in_jet(jets.at(i),jets1.at(0))) std::cout<<"ja"<<std::endl;

  /*fastjet::PseudoJet parent1, parent2;
  for(int i=0;i<jets.size();i++) {
    //std::cout<<i<<" "<<cs->has_parents(cs->jets().at(i),parent1, parent2)<<std::endl;
    if(!cs->has_parents(cs->jets().at(i),parent1, parent2)) nhistos++;
    else break;
 }
  


  for(int i=0;i<cs->jets().size();i++){
    TString hname2 = TString::Format("JetDisplay_jet%i", i);
    Book(TH2F(hname2,"Jet event display",50,-PI,PI,50,-PI,PI));
    for(int j=nhistos;j<cs->jets().size();j++){
      // std::cout<<history.at(i)<<std::endl;
      if(j==history.at(i)) clustered_jets.push_back(cs->jets().at(i));
    }
    for(int j=0;j<clustered_jets.size();j++) ((TH2D*)Hist(hname2))->Fill(clustered_jets.at(j).phi_std(),clustered_jets.at(j).eta(),clustered_jets.at(j).pt());
    }
  */
  /* for(int i=0;i<jets.size();i++){
    std::cout<<i<<" "<<jets.at(i).constituents().size()<<std::endl;
    }*/


  fastjet::PseudoJet parent1, parent2;

  nhistos=cs->n_particles();
  std::cout<<cs->jets().size()-nhistos<<std::endl;
  Book(TH2F("JetDisplay_pf_all","Jet event display",50,-PI,PI,50,-PI,PI));
  for(int p=0;p<cs->jets().size()/*-nhistos+1*/;p++) {
    //TString hname2 = TString::Format("JetDisplay_jet%i", p);
      //  Book(TH2F(hname2,"Jet event display",50,-PI,PI,50,-PI,PI));
    }
  int nhistos2=cs->jets().size()-nhistos+1;
  for(int i=0;i<nhistos;i++){
    ((TH2D*)Hist("JetDisplay_pf_all"))->Fill(jets.at(i).phi_std(),jets.at(i).eta(),jets.at(i).pt());
  }
  int counter=1;
  fastjet::PseudoJet oldjet;
  std::cout<<"pt "<<oldjet.pt()<<std::endl;
  for(int j=nhistos;j<jets.size();j++)

     {
       if(!Jet_contains_hard_particles(jets.at(j),0.01)) continue;
       for(int b=0;b<nhistos;b++){
	 if(cs->object_in_jet(jets.at(b),jets.at(j)) && jets.at(j).m()>5/*&& jets.at(j).constituents().size()==2*/){
	   oldjet=jets.at(j);
	   jets.at(b).set_user_index(j);
	   std::cout<<"J "<<oldjet.pt()<<std::endl;
	   clustered_jets.push_back(jets.at(b));
	   //nhistos--;
	   //jets.erase(jets.begin()+b);
	   // std::cout<<"deleted"<<std::endl;
	    
	   
	 }
	 }
       if(clustered_jets.size()!=0) {
	 
	 for(int i=0;i<clustered_jets.size();i++){
	   
	   TString hname2 = TString::Format("JetDisplay_jet%i", counter);
	   
	    Book(TH2F(hname2,"Jet event display",50,-PI,PI,50,-PI,PI));
	   ((TH2D*)Hist(hname2))->Fill(clustered_jets.at(i).phi_std(),clustered_jets.at(i).eta(),clustered_jets.at(i).pt());
	   
	 }
	 counter++;
	 clustered_jets.clear();
       }
     }
 for(int o=0;o<decay_products.size();o++) ((TH2D*)Hist("JetDisplay_decay"))->Fill(decay_products[o].phi(),decay_products[o].eta(),decay_products[o].pt());


  /*if(hotvr_jets.size()<2)  throw SError( SError::SkipEvent );
 
  if(hotvr_jets.at(0).user_info<HOTVRinfo>().nsubjets()!=3)   throw SError( SError::SkipEvent );

 if(hotvr_jets.at(1).user_info<HOTVRinfo>().nsubjets()!=3) throw SError( SError::SkipEvent );

  if(hotvr_jets.at(0).m()<140 || hotvr_jets.at(0).m()>220 ) throw SError( SError::SkipEvent );
 
  if(hotvr_jets.at(1).m()<140 || hotvr_jets.at(1).m()>220 ) throw SError( SError::SkipEvent );

   if(abs(hotvr_jets.at(0).phi_std())>2. ||  abs(hotvr_jets.at(1).phi_std())>2.)  throw SError( SError::SkipEvent );

  if(!(hotvr_jets.at(0).pt()>200 && hotvr_jets.at(0).pt()<450) && !(hotvr_jets.at(0).pt()>800 && hotvr_jets.at(0).pt()<1000))  throw SError( SError::SkipEvent );

 


  IR_Saftey->add_grid(parts,75,75);
  clustering->Reset();
  hotvr_jets=clustering->get_clustered_jets(parts);
  soft_jets=clustering->get_soft_jets();
  rejected_jets=clustering->get_rejected_jets();
  ca_jets=clustering_ca->get_clustered_jets(parts);
  akt_jets=clustering_akt->get_clustered_jets(parts);
  */
 
 
   
    
    //matching_radius=1.5;
 
 
    


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



 std::cout	  <<"\r"<<"Progress "<<(m_counter2+m_counter)/2./_event_max*100.<<"%"<<std::flush;
 
  clustering->Reset();
  clustering_ca->Reset();
  clustering_akt->Reset();
  
  return;
  
}


