#include <time.h>
#include "include/HOTVRHists.h"
#include "SFrameTools/include/EventCalc.h"
//#include "include/ZprimeFullHadTools.h"
#include "NtupleWriter/include/JetProps.h"
//#include "include/TopTagfunctions.h"

#include "fastjet/ClusterSequenceArea.hh"
 #include "fastjet/ClusterSequencePassiveArea.hh"
#include "NtupleWriter/include/JetPropsPseudo.h"
#include "fastjet/contrib/ClusteringVetoPlugin.hh"
#include "fastjet/contrib/HHTopTagger.hh"
#include "fastjet/contrib/HOTVR.hh"
#include <iostream>
#include "TH3.h"
using namespace std;
using namespace fastjet;
using namespace contrib;
HOTVRHists::HOTVRHists(const char* name) : BaseHists(name)
{
  //  std::cout<<name<<std::endl;
;
  // named default constructor
   
}

HOTVRHists::~HOTVRHists()
{
  // default destructor, does nothing
}

void HOTVRHists::Init()
{
  // book all histograms here

 
 
 
  Book( TH1F( "pT_denominator"," p_{T} topjets (selection)",75,0,3000));
  Book( TH1F( "eta_denominator"," eta topjets (selection)",20,-3,3));
 
  Book( TH1F( "pT_nominator"," p_{T} topjets w mass (hotvr tagged)",75,0,3000));
  Book( TH1F( "eta_nominator"," eta topjets w mass (hotvr tagged)",20,-3,3));
    Book( TH2D("nsub_nominator","eff nsub",75,0,3000,100,-0.005,0.995));
  Book( TH2D("nsub_denominator","eff nsub",75,0,3000,100,-0.005,0.995));
  
  
  

  
  
 

}

void HOTVRHists::SetIdVersion(TString s)
{
  idVersion=s;
}

void HOTVRHists::Fill(fastjet::PseudoJet jet,double jet_radius, fastjet::PseudoJet matched_jet, double weight)
{

}
	  
	   
  
      
void HOTVRHists::Fill_nominator(fastjet::PseudoJet jet,double jet_radius, fastjet::PseudoJet matched_jet,double weight){
  Hist("pT_nominator")->Fill(matched_jet.pt(),weight);
    if(jet.constituents().size()>0){
	   double tau1,tau2,tau3,tau4;
	   JetPropsPseudo jp(&jet);
	  
	   tau2 = jp.GetNsubjettiness(2, Njettiness::onepass_kt_axes, 1., jet_radius);
	   tau3 = jp.GetNsubjettiness(3, Njettiness::onepass_kt_axes, 1., jet_radius );
	   
    
    for(int roc=0;roc<100;roc++){
      if(tau3/tau2<roc/100.) ((TH2D*)Hist("nsub_nominator"))->Fill(matched_jet.pt(),roc/100.,weight);
    }	   

    

    }
}


    void HOTVRHists::Fill_denominator(fastjet::PseudoJet denominator_jet,double weight){
      
      Hist("pT_denominator")->Fill(denominator_jet.pt(),weight); //fill denominator hists for efficiecies
      for(int roc=0;roc<100;roc++) {((TH2D*)Hist("nsub_denominator"))->Fill(denominator_jet.pt(),roc/100.,weight);}//Fill denominator hist for nsubjettiness scan
      
    }

   

  

 




void HOTVRHists::Finish()
{
  // final calculations, like division and addition of certain histograms
  EventCalc* calc = EventCalc::Instance();
  bool IsRealData = calc->IsRealData();
  if (IsRealData){
    Hist("N_pv_perLumiBin")->Divide( Hist("N_pv_perLumiBin"), Hist("N_events_perLumiBin"));
    Hist( "N_pv_perLumiBin")->GetYaxis()->SetTitle("Events/Lumi");
  }

}

