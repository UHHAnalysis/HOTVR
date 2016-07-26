#include <time.h>
#include "include/ClusteringHists.h"
#include "SFrameTools/include/EventCalc.h"
//#include "include/ZprimeFullHadTools.h"
#include "NtupleWriter/include/JetProps.h"
//#include "include/TopTagfunctions.h"
#include "include/Clusteringinfo.hh"
#include "fastjet/ClusterSequenceArea.hh"
 #include "fastjet/ClusterSequencePassiveArea.hh"
#include "NtupleWriter/include/JetPropsPseudo.h"
#include "fastjet/contrib/ClusteringVetoPlugin.hh"
#include "fastjet/contrib/HHTopTagger.hh"
#include "fastjet/contrib/HOTVR.hh"
#include "fastjet/contrib/Jetinfo.hh"
#include <iostream>
#include "TH3.h"
using namespace std;
using namespace fastjet;
using namespace contrib;
ClusteringHists::ClusteringHists(const char* name) : BaseHists(name)
{
  //  std::cout<<name<<std::endl;
;
  // named default constructor
   
}

ClusteringHists::~ClusteringHists()
{
  // default destructor, does nothing
}

void ClusteringHists::Init()
{
  // book all histograms here

 
 
 
  
  
  
  

  
  
  Book( TH1D( "weight", ";weight;Events", 400, 0., 1000.));
  
  Book(TH1D("Clustering","Clustering algorithmus",10,0,10));
  Book(TH1D("Clustering_radius","Clustering radius",30,0,3));
  Book(TH1D("minimum_jet_pT","minimum jet p_{T}",100,0,500));
  Book(TH1D("Rho","Rho",100,0,1000));
  Book(TH1D("Mu","Mu",100,0,200));
  Book(TH1D("Theta","Theta",10,0,1));
  Book(TH1D("Rmin","R_{min}",30,0,3));
  Book(TH1D("Rmax","R_{max}",30,0,3));
  Book(TH1D("Modus","Modus",10,0,10));
  Book(TH1D("pT_cut","p_{T} cut subjets",100,0,100));
  
  
 
 
  

}



void ClusteringHists::Fill(fastjet::PseudoJet jet, double weight)
{
  if(jet.has_user_info<Jetinfo>()){
    Hist("Clustering")->Fill(jet.user_info<Jetinfo>().clustering().c_str(),weight);
    Hist("Clustering_radius")->Fill(jet.user_info<Jetinfo>().radius(),weight);
    Hist("minimum_jet_pT")->Fill(jet.user_info<Jetinfo>().ptmin(),weight);
    Hist("Rho")->Fill(jet.user_info<Jetinfo>().rho(),weight);
    Hist("Mu")->Fill(jet.user_info<Jetinfo>().mu(),weight);
    Hist("Theta")->Fill(jet.user_info<Jetinfo>().theta(),weight);
    Hist("Rmin")->Fill(jet.user_info<Jetinfo>().rmin(),weight);
    Hist("Rmax")->Fill(jet.user_info<Jetinfo>().rmax(),weight);
     Hist("Modus")->Fill(jet.user_info<Jetinfo>().modus().c_str(),weight);
    Hist("pT_cut")->Fill(jet.user_info<Jetinfo>().ptcut(),weight);

  }
  else std::cout<<"no clustering info"<<std::endl;
 
} 
   
	  
	   
  
     




   

  

 




void ClusteringHists::Finish()
{
  // final calculations, like division and addition of certain histograms
  EventCalc* calc = EventCalc::Instance();
  bool IsRealData = calc->IsRealData();
  if (IsRealData){
    Hist("N_pv_perLumiBin")->Divide( Hist("N_pv_perLumiBin"), Hist("N_events_perLumiBin"));
    Hist( "N_pv_perLumiBin")->GetYaxis()->SetTitle("Events/Lumi");
  }

}

