// $Id: RocCycle.cxx,v 1.10 2012/12/07 14:21:51 peiffer Exp $

#include <iostream>


#include "include/HOTVRCycle.h"
#include "include/SubstructureHists.h"
#include "include/HOTVRHists.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/contrib/HOTVR.hh"
#include "SFrameAnalysis/include/SelectionModules.h"
//#include "SFrameAnalysis/include/HypothesisHists.h"
using namespace std;
using namespace fastjet;
using namespace contrib;
ClassImp( HOTVRCycle );

HOTVRCycle::HOTVRCycle()
   : AnalysisCycle() {

  // constructor, declare additional variables that should be 
  // obtained from the steering-xml file
  
  // set the integrated luminosity per bin for the lumi-yield control plots
  SetIntLumiPerBin(500.);

}

HOTVRCycle::~HOTVRCycle() 
{
  // destructor
}

void HOTVRCycle::BeginCycle() throw( SError ) 
{
  // Start of the job, general set-up and definition of 
  // objects are done here

  // Important: first call BeginCycle of base class
  AnalysisCycle::BeginCycle();

  return;

}

void HOTVRCycle::EndCycle() throw( SError ) 
{
  // clean-up, info messages and final calculations after the analysis

  
  // call the base cycle class for all standard methods
  AnalysisCycle::EndCycle();

  return;

}

void HOTVRCycle::BeginInputData( const SInputData& id ) throw( SError ) 
{
  // declaration of histograms and selections.
  // AnalysisCyle expects Selections and HistCollections to be registered here.
  // Their memory will be released in AnalysisCycle::EndInputData.

  // Important: first call BeginInputData of base class
  AnalysisCycle::BeginInputData( id );

 

  // ---------------- set up the histogram collections --------------------
  RegisterHistCollection( new SubstructureHists("HOTVR_hists200-400"));
  RegisterHistCollection( new SubstructureHists("HOTVR_hists400-600"));
  RegisterHistCollection( new SubstructureHists("HOTVR_hists600-800"));
  RegisterHistCollection( new HOTVRHists("HOTVR_hists_eff"));
   

  // important: initialise histogram collections after their definition
  InitHistos();

}

void HOTVRCycle::EndInputData( const SInputData& id ) throw( SError ) 
{
  AnalysisCycle::EndInputData( id );
}

void HOTVRCycle::BeginInputFile( const SInputData& id ) throw( SError ) 
{
  // Connect all variables from the Ntuple file with the ones needed for the analysis
  // The variables are commonly stored in the BaseCycleContaincer
  // important: call to base function to connect all variables to Ntuples from the input tree
  AnalysisCycle::BeginInputFile( id );
}

void HOTVRCycle::ExecuteEvent( const SInputData& id, Double_t weight) throw( SError ) 
{
  // this is the most important part: here the full analysis happens
  // user should implement selections, filling of histograms and results

  // first step: call Execute event of base class to perform basic consistency checks
  // also, the good-run selection is performed there and the calculator is reset
  AnalysisCycle::ExecuteEvent( id, weight );

  // get the histogram collections. NOTE: this could be done more performant by making
  // all thse BaseHists* vairables private member variables of RocCycle and
  // setting them in BeginInputData. Then, there is no need here to call GetHistColletion ...

  BaseHists* HOTVR_hists200 = GetHistCollection("HOTVR_hists200-400");
   BaseHists* HOTVR_hists400 = GetHistCollection("HOTVR_hists400-600");
   BaseHists* HOTVR_hists600 = GetHistCollection("HOTVR_hists600-800");
 BaseHists* HOTVR_hists_eff = GetHistCollection("HOTVR_hists_eff");


  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();
  bool IsRealData = calc->IsRealData();
 

  double rho(600.);
  double mu(30.), theta(0.7), pt_cut(20.);
  double min_r(0.1), max_r(1.5);
 
 
  // important: get the event weight
 

  //implement tagger choice

  //new class tagger mit constructor Tagger("cms")
  /* double mjet,mmin;
	  int nsubjets;
	  std::vector<fastjet::PseudoJet> subjets;
	  fastjet::PseudoJet CMSjet;
	 
	  // jet=varjets.at(bestjetindex2);
	  double chi=0;
	  double Mmicrojet=0;
	  int Nmicrojets=0;
	  std::vector<fastjet::PseudoJet> microjets;
	  
	  double microconesize=0.3;
	  if(Had_Tops[j].pt()>500) microconesize=0.2;
	  //if(Had_Tops[j].pt()>700) microconesize=0.1;
	    bool CMStag=false;
	    CMStag=CMSTopTagFull_pseudo_CA(jet,3.0,0.05,0.4,0.0004,mjet, nsubjets,mmin,subjets,CMSjet);
	  
	   
	   //chi = Showerdeconstruction_taggerV2->ChiMicro_pseudo(jet,Nmicrojets,Mmicrojet,microconesize,microjets);
	    double Rmin,mass_Rmin,pt_Rmin,mass_diff,pt_for_exp,Rmin_exp,fw_out;
	  double m12,m13, m23, m123, m_pruned, m_unfiltered;
	  double fw=0.15;
	  // if(HepTopTagFull_pseudo(jet, m12, m13, m23, fw, m123, m_pruned, m_unfiltered))   CMStag=true;
	  // if(MultiRTopTag_pseudo(jet, calc->GetPFParticles(),fw, Rmin,mass_Rmin,pt_Rmin,mass_diff,pt_for_exp,Rmin_exp,fw_out) &&Rmin-Rmin_exp<0.07 && mass_Rmin>140 && mass_Rmin<250) CMStag=true;
	  //  if(MultiRTopTag_pseudo(jet, calc->GetPFParticles(), Rmin,mass_Rmin,pt_Rmin,mass_diff,pt_for_exp,Rmin_exp) &&Rmin-Rmin_exp<0.07 && mass_Rmin>140 && mass_Rmin<250) CMStag=true;
	  // if(log(chi)>2) CMStag=true;
	  // else CMStag=false;
	  //std::cout<<"tag "<<CMStag<<std::endl;*/



  int bestjetindex2=-1;
  matching= new Matching();
  std::vector<GenParticle>* genparticles = calc->GetGenParticles();
  matching->Run_matching(genparticles);
  std::vector<fastjet::PseudoJet> parts= matching->get_hadrons();
 
  double jet_radius=0.8;
  double ptmin=100;
  std::vector<fastjet::PseudoJet> hotvr_jets;
  // hotvr_jets=clustering->get_clustered_jets(parts,Clustering::E_algorithm::e_akt,jet_radius,ptmin);
  hotvr_jets=clustering->get_clustered_hotvr_jets(parts,Clustering::E_algorithm::e_akt,  ptmin, rho,min_r ,max_r, mu, theta, pt_cut);
 
  std::vector<fastjet::PseudoJet> denominator_jets=matching->get_denominator_jets(id.GetVersion());
((SubstructureHists*)HOTVR_hists200)->SetIdVersion(id.GetVersion());

  ((SubstructureHists*)HOTVR_hists400)->SetIdVersion(id.GetVersion());

((SubstructureHists*)HOTVR_hists600)->SetIdVersion(id.GetVersion());

 for(int j=0;j<denominator_jets.size();j++){
   ((HOTVRHists*) HOTVR_hists_eff)->Fill_denominator(denominator_jets[j],weight);
   fastjet::PseudoJet matched_jet=matching->get_closest_jet(hotvr_jets,denominator_jets[j]);
     
     //if(matching->IsMatched(idVersion,hotvr_jets[i],hotvr_jets[i].user_info<HOTVRinfo>().radius(),matched_jet,denominator_jets)) std::cout<<"MATCHING"<<std::endl;;
     
    
     //std::cout<<" number of denominator jets "<<denominator_jets.size()<<std::endl;
     double matching_radius;
     if(matched_jet.has_user_info<HOTVRinfo>()) matching_radius=matched_jet.user_info<HOTVRinfo>().radius();
     else matching_radius=jet_radius;
     if (matching->IsMatched(matched_jet,matching_radius,denominator_jets[j])) 
       {
	 if(denominator_jets[j].pt()>200 && denominator_jets[j].pt()<400)  ((SubstructureHists*)HOTVR_hists200)->Fill(matched_jet,matching_radius,denominator_jets[j],weight);
	 if(denominator_jets[j].pt()>400 && denominator_jets[j].pt()<600)  ((SubstructureHists*)HOTVR_hists400)->Fill(matched_jet,matching_radius,denominator_jets[j],weight);
	  if(denominator_jets[j].pt()>600 && denominator_jets[j].pt()<800)  ((SubstructureHists*)HOTVR_hists600)->Fill(matched_jet,matching_radius,denominator_jets[j],weight);
	  ((HOTVRHists*) HOTVR_hists_eff)->Fill_nominator(matched_jet, jet_radius,denominator_jets[j],weight);
       }
     
     


     
     //((HOTVRHists*)HOTVR_hists)->SetIdVersion(id.GetVersion());
     
     //HOTVR_hists->Fill();
     

   
 }
 //for(int i=0;i<denominator_jets.size();i++) ((HOTVRHists*) HOTVR_hists_eff)->Fill_denominator(denominator_jets[i],weight);
  
  
  return;
  
}


