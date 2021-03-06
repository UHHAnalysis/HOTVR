// $Id: RocCycle.cxx,v 1.10 2012/12/07 14:21:51 peiffer Exp $

#include <iostream>


#include "include/HOTVRCycle.h"
#include "include/SubstructureHists.h"
#include "include/HOTVRHists.h"
#include "include/ClusteringHists.h"
#include "include/Clusteringinfo.hh"
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
  
  RegisterHistCollection( new ClusteringHists("Clustering"));
  RegisterHistCollection( new SubstructureHists("HOTVR_hists200-400"));
  RegisterHistCollection( new SubstructureHists("HOTVR_hists400-600"));
  RegisterHistCollection( new SubstructureHists("HOTVR_hists600-800"));
  RegisterHistCollection( new SubstructureHists("HOTVR_hists800-1000"));
  RegisterHistCollection( new SubstructureHists("HOTVR_hists1000-1200"));
  RegisterHistCollection( new SubstructureHists("HOTVR_hists_all"));

  RegisterHistCollection( new SubstructureHists("HOTVR_hists200-400_beforetag"));
  RegisterHistCollection( new SubstructureHists("HOTVR_hists400-600_beforetag"));
  RegisterHistCollection( new SubstructureHists("HOTVR_hists600-800_beforetag"));
  RegisterHistCollection( new SubstructureHists("HOTVR_hists800-1000_beforetag"));
  RegisterHistCollection( new SubstructureHists("HOTVR_hists1000-1200_beforetag"));
  RegisterHistCollection( new SubstructureHists("HOTVR_hists_all_beforetag"));

  RegisterHistCollection( new SubstructureHists("HOTVR_hists200-400_nsub3"));
  RegisterHistCollection( new SubstructureHists("HOTVR_hists400-600_nsub3"));
  RegisterHistCollection( new SubstructureHists("HOTVR_hists600-800_nsub3"));
  RegisterHistCollection( new SubstructureHists("HOTVR_hists800-1000_nsub3"));
  RegisterHistCollection( new SubstructureHists("HOTVR_hists1000-1200_nsub3"));

  RegisterHistCollection( new SubstructureHists("HOTVR_hists200-400_ptfraction"));
  RegisterHistCollection( new SubstructureHists("HOTVR_hists400-600_ptfraction"));
  RegisterHistCollection( new SubstructureHists("HOTVR_hists600-800_ptfraction"));
  RegisterHistCollection( new SubstructureHists("HOTVR_hists800-1000_ptfraction"));
  RegisterHistCollection( new SubstructureHists("HOTVR_hists1000-1200_ptfraction"));
  

  RegisterHistCollection( new HOTVRHists("HOTVR_hists_eff"));
   

  // important: initialise histogram collections after their definition
  InitHistos();
  toptagger=new TopTagger;
  clustering = new Clustering(m_clustering);
  // clustering->jets(m_clustering);
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
  BaseHists* Clustering = GetHistCollection("Clustering");
  BaseHists* HOTVR_hists200 = GetHistCollection("HOTVR_hists200-400");
   BaseHists* HOTVR_hists400 = GetHistCollection("HOTVR_hists400-600");
   BaseHists* HOTVR_hists600 = GetHistCollection("HOTVR_hists600-800");
   BaseHists* HOTVR_hists800 = GetHistCollection("HOTVR_hists800-1000");
   BaseHists* HOTVR_hists1000 = GetHistCollection("HOTVR_hists1000-1200");
   BaseHists* HOTVR_hists_all = GetHistCollection("HOTVR_hists_all");
   BaseHists* HOTVR_hists200_beforetag = GetHistCollection("HOTVR_hists200-400_beforetag");
   BaseHists* HOTVR_hists400_beforetag = GetHistCollection("HOTVR_hists400-600_beforetag");
   BaseHists* HOTVR_hists600_beforetag = GetHistCollection("HOTVR_hists600-800_beforetag");
   BaseHists* HOTVR_hists800_beforetag = GetHistCollection("HOTVR_hists800-1000_beforetag");
   BaseHists* HOTVR_hists1000_beforetag = GetHistCollection("HOTVR_hists1000-1200_beforetag");
   BaseHists* HOTVR_hists_all_beforetag = GetHistCollection("HOTVR_hists_all_beforetag");

   BaseHists* HOTVR_hists200_nsub3 = GetHistCollection("HOTVR_hists200-400_nsub3");
   BaseHists* HOTVR_hists400_nsub3 = GetHistCollection("HOTVR_hists400-600_nsub3");
   BaseHists* HOTVR_hists600_nsub3 = GetHistCollection("HOTVR_hists600-800_nsub3");
   BaseHists* HOTVR_hists800_nsub3 = GetHistCollection("HOTVR_hists800-1000_nsub3");
   BaseHists* HOTVR_hists1000_nsub3 = GetHistCollection("HOTVR_hists1000-1200_nsub3");

   BaseHists* HOTVR_hists200_ptfraction = GetHistCollection("HOTVR_hists200-400_ptfraction");
   BaseHists* HOTVR_hists400_ptfraction = GetHistCollection("HOTVR_hists400-600_ptfraction");
   BaseHists* HOTVR_hists600_ptfraction = GetHistCollection("HOTVR_hists600-800_ptfraction");
   BaseHists* HOTVR_hists800_ptfraction = GetHistCollection("HOTVR_hists800-1000_ptfraction");
   BaseHists* HOTVR_hists1000_ptfraction = GetHistCollection("HOTVR_hists1000-1200_ptfraction");
 BaseHists* HOTVR_hists_eff = GetHistCollection("HOTVR_hists_eff");


  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();
  bool IsRealData = calc->IsRealData();
 

  int bestjetindex2=-1;
  matching= new Matching();
  std::vector<GenParticle>* genparticles = calc->GetGenParticles();
  matching->Run_matching(genparticles);
  std::vector<fastjet::PseudoJet> parts= matching->get_hadrons();
  
  
  std::vector<fastjet::PseudoJet> hotvr_jets;
  hotvr_jets=clustering->get_clustered_jets(parts);
  
  
   
  std::vector<fastjet::PseudoJet> denominator_jets=matching->get_denominator_jets(id.GetVersion());
  
  
  for(int j=0;j<denominator_jets.size();j++){
    
    ((HOTVRHists*) HOTVR_hists_eff)->Fill_denominator(denominator_jets[j],weight);
    if(hotvr_jets.size()==0) continue;
    
    fastjet::PseudoJet matched_jet=matching->get_closest_jet(hotvr_jets,denominator_jets[j]);
    if(!matched_jet.has_user_info<HOTVRinfo>()) matched_jet.set_user_info(new HOTVRinfo(matched_jet,matched_jet.constituents()));
    double matching_radius;
    if(matched_jet.has_user_info<HOTVRinfo>()) matching_radius=matched_jet.user_info<HOTVRinfo>().radius();
    
   
    
    //matching_radius=1.5;
    if (matching->IsMatched(matched_jet,matching_radius,denominator_jets[j])) 
      {
	((ClusteringHists*)Clustering)->Fill(matched_jet,weight);
	
	
	if(denominator_jets[j].pt()>200 && denominator_jets[j].pt()<400)  ((SubstructureHists*)HOTVR_hists200_beforetag)->Fill(matched_jet,matching_radius,denominator_jets[j],weight);
	if(denominator_jets[j].pt()>200 && denominator_jets[j].pt()<400 &&matched_jet.user_info<HOTVRinfo>().nsubjets()>2 )  ((SubstructureHists*)HOTVR_hists200_nsub3)->Fill(matched_jet,matching_radius,denominator_jets[j],weight);
	if(denominator_jets[j].pt()>200 && denominator_jets[j].pt()<400 &&matched_jet.user_info<HOTVRinfo>().ptfraction(0)<0.8 )  ((SubstructureHists*)HOTVR_hists200_ptfraction)->Fill(matched_jet,matching_radius,denominator_jets[j],weight);
	if(denominator_jets[j].pt()>400 && denominator_jets[j].pt()<600)  ((SubstructureHists*)HOTVR_hists400_beforetag)->Fill(matched_jet,matching_radius,denominator_jets[j],weight);
	if(denominator_jets[j].pt()>400 && denominator_jets[j].pt()<600 &&matched_jet.user_info<HOTVRinfo>().nsubjets()>2 )  ((SubstructureHists*)HOTVR_hists400_nsub3)->Fill(matched_jet,matching_radius,denominator_jets[j],weight);
	if(denominator_jets[j].pt()>400 && denominator_jets[j].pt()<600 &&matched_jet.user_info<HOTVRinfo>().ptfraction(0)<0.8 )  ((SubstructureHists*)HOTVR_hists400_ptfraction)->Fill(matched_jet,matching_radius,denominator_jets[j],weight);
	 if(denominator_jets[j].pt()>600 && denominator_jets[j].pt()<800)  ((SubstructureHists*)HOTVR_hists600_beforetag)->Fill(matched_jet,matching_radius,denominator_jets[j],weight);
	 if(denominator_jets[j].pt()>600 && denominator_jets[j].pt()<800 &&matched_jet.user_info<HOTVRinfo>().nsubjets()>2 )  ((SubstructureHists*)HOTVR_hists600_nsub3)->Fill(matched_jet,matching_radius,denominator_jets[j],weight);
	 if(denominator_jets[j].pt()>600 && denominator_jets[j].pt()<800 &&matched_jet.user_info<HOTVRinfo>().ptfraction(0)<0.8 )  ((SubstructureHists*)HOTVR_hists600_ptfraction)->Fill(matched_jet,matching_radius,denominator_jets[j],weight);
	 if(denominator_jets[j].pt()>800 && denominator_jets[j].pt()<1000)  ((SubstructureHists*)HOTVR_hists800_beforetag)->Fill(matched_jet,matching_radius,denominator_jets[j],weight);
	 if(denominator_jets[j].pt()>800 && denominator_jets[j].pt()<1000 &&matched_jet.user_info<HOTVRinfo>().nsubjets()>2 )  ((SubstructureHists*)HOTVR_hists800_nsub3)->Fill(matched_jet,matching_radius,denominator_jets[j],weight);
	 if(denominator_jets[j].pt()>800 && denominator_jets[j].pt()<1000 &&matched_jet.user_info<HOTVRinfo>().ptfraction(0)<0.8 )  ((SubstructureHists*)HOTVR_hists800_ptfraction)->Fill(matched_jet,matching_radius,denominator_jets[j],weight);
	 if(denominator_jets[j].pt()>1000 && denominator_jets[j].pt()<1200)  ((SubstructureHists*)HOTVR_hists1000_beforetag)->Fill(matched_jet,matching_radius,denominator_jets[j],weight);
	 if(denominator_jets[j].pt()>1000 && denominator_jets[j].pt()<1200 &&matched_jet.user_info<HOTVRinfo>().nsubjets()>2 )  ((SubstructureHists*)HOTVR_hists1000_nsub3)->Fill(matched_jet,matching_radius,denominator_jets[j],weight);
	 if(denominator_jets[j].pt()>1000 && denominator_jets[j].pt()<1200 &&matched_jet.user_info<HOTVRinfo>().ptfraction(0)<0.8 )  ((SubstructureHists*)HOTVR_hists1000_ptfraction)->Fill(matched_jet,matching_radius,denominator_jets[j],weight);


	 ((SubstructureHists*)HOTVR_hists_all_beforetag)->Fill(matched_jet,matching_radius,denominator_jets[j],weight);
	 if(toptagger->Is_tagged(m_tagger, matched_jet)){
	   if(denominator_jets[j].pt()>200 && denominator_jets[j].pt()<400)  ((SubstructureHists*)HOTVR_hists200)->Fill(matched_jet,matching_radius,denominator_jets[j],weight);
	   if(denominator_jets[j].pt()>400 && denominator_jets[j].pt()<600)  ((SubstructureHists*)HOTVR_hists400)->Fill(matched_jet,matching_radius,denominator_jets[j],weight);
	   if(denominator_jets[j].pt()>600 && denominator_jets[j].pt()<800)  ((SubstructureHists*)HOTVR_hists600)->Fill(matched_jet,matching_radius,denominator_jets[j],weight);
	   if(denominator_jets[j].pt()>800 && denominator_jets[j].pt()<1000)  ((SubstructureHists*)HOTVR_hists800)->Fill(matched_jet,matching_radius,denominator_jets[j],weight);
	   if(denominator_jets[j].pt()>1000 && denominator_jets[j].pt()<1200)  ((SubstructureHists*)HOTVR_hists1000)->Fill(matched_jet,matching_radius,denominator_jets[j],weight);
	   ((SubstructureHists*)HOTVR_hists_all)->Fill(matched_jet,matching_radius,denominator_jets[j],weight);
	   ((HOTVRHists*) HOTVR_hists_eff)->Fill_nominator(matched_jet, matching_radius,denominator_jets[j],weight);
	 }
      }
     
     


     
     
     

   
 }

 //if(hotvr_jets.size()!=0) delete hotvr_jets.at(0).associated_cluster_sequence();
 delete matching;
 //delete clustering;
 clustering->Reset();
 // delete toptagger;
  return;
  
}


