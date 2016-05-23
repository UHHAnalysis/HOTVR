#include <time.h>
#include "include/SubstructureHists.h"
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
SubstructureHists::SubstructureHists(const char* name) : BaseHists(name)
{
  //  std::cout<<name<<std::endl;
;
  // named default constructor
   
}

SubstructureHists::~SubstructureHists()
{
  // default destructor, does nothing
}

void SubstructureHists::Init()
{
  // book all histograms here

 
 
 
  
  
  
  

  
  
  Book( TH1D( "weight", ";weight;Events", 400, 0., 1000.));
  
  
  
  //Variable R

  Book(TH2F("VariableR","VariableR",1000,0,2000,20,0,40));
 
  Book(TH1D("jet_pT","jet p_{T}",100,0,1000));
 
  Book(TH1D("mmin","massjump mmin",100,0,200));
 

  Book(TH1D("mjet","mjet",100,0,300));
 

  

 
   
 
 
  Book(TH1D("nsubjets","nsubjets massjump",20,-0.5,19.5));
  
 
    
  Book(TH1D("tau1","tau1",100,0,1.0));
  Book(TH1D("tau2","tau2",100,0,1.0));
   Book(TH1D("tau3","tau3",100,0,1.0));
  Book(TH1D("tau4","tau4",100,0,1.0));
  Book(TH1D("tau2tau1","tau2tau1",100,0,1.0));
  Book(TH1D("tau3tau2","tau2tau1",100,0,1.0));
  Book(TH1D("tau3tau1","tau2tau1",100,0,1.0));

 
  
  //subjetmass
  Book(TH1D("subjet01","subjet01 mass",100,0,200));
  Book(TH1D("subjet02","subjet02 mass",100,0,200));
  Book(TH1D("subjet12","subjet12 mass",100,0,200));
 
 
  
   //pt
  

    
  //pt subjets
  Book( TH1F( "pT_subjet1"," p_{T} leading subjet",100,0,2000));
  Book( TH1F( "pT_subjetfrac1"," p_{T} leading subjet",100,0,1));
  Book( TH1F( "pT_subjet1_ly"," p_{T} leading subjet",100,0,2000));
  Book( TH1F( "pT_subjet2","p_{T} 2nd subjet",100,0,2000));
  Book( TH1F( "pT_subjetfrac2"," p_{T} leading subjet",100,0,1));
  Book( TH1F( "pT_subjet2_ly","p_{T} 2nd subjet",100,0,2000));
  Book( TH1F( "pT_subjet3","p_{T} 3rd subjet",100,0,1000));
  Book( TH1F( "pT_subjetfrac3"," p_{T} leading subjet",100,0,3));
  Book( TH1F( "pT_subjet3_ly","p_{T} 3rd subjet",100,0,1000));
  Book( TH1F( "pT_subjet4","p_{T} 4th subjet",100,0,800));
  Book( TH1F( "pT_subjetfrac4","p_{T} 4th subjet",100,0,800));
  Book( TH1F( "pT_subjet4_ly","p_{T} 4th subjet",100,0,800));
  Book( TH1F( "eta_subjet1","#eta leading subjet",100,-3,3));
  Book( TH1F( "eta_subjet1_ly","#eta leading subjet",100,-3,3));
  Book( TH1F( "eta_subjet2","#eta 2nd subjet",100,-3,3));
  Book( TH1F( "eta_subjet2_ly","#eta 2nd subjet",100,-3,3));
  Book( TH1F( "eta_subjet3","#eta 3rd subjet",100,-3,3));
  Book( TH1F( "eta_subjet3_ly","#eta 3rd subjet",100,-3,3));
  Book( TH1F( "eta_subjet4","#eta 4th subjet",100,-3,3));
  Book( TH1F( "eta_subjet4_ly","#eta 4th subjet",100,-3,3));
  Book( TH1F( "phi_subjet1","#phi leading subjet",100,-M_PI,M_PI));
  Book( TH1F( "phi_subjet1_ly","#phi leading subjet",100,-M_PI,M_PI));
  Book( TH1F( "phi_subjet2","#phi 2nd subjet",100,-M_PI,M_PI));
  Book( TH1F( "phi_subjet2_ly","#phi 2nd subjet",100,-M_PI,M_PI));
  Book( TH1F( "phi_subjet3","#phi 3rd subjet",100,-M_PI,M_PI));
  Book( TH1F( "phi_subjet3_ly","#phi 3rd subjet",100,-M_PI,M_PI));
  Book( TH1F( "phi_subjet4","#phi 4th subjet",100,-M_PI,M_PI));
  Book( TH1F( "phi_subjet4_ly","#phi 4th subjet",100,-M_PI,M_PI));
  Book( TH1F("mass_subjet1","mass subjet1",100,0,250));
  Book( TH1F("mass_subjet2","mass subjet2",100,0,250));
  Book( TH1F("mass_subjet3","mass subjet3",100,0,250));
  Book( TH1F("mass_subjet4","mass subjet4",100,0,250));
 
  Book(TH1D("ptfraction1","p_{T,sub1}/p_{T,jet} [GeV]",20,0,1));
 
  
  clustering=new Clustering();

}

void SubstructureHists::SetIdVersion(TString s)
{
  idVersion=s;
}

void SubstructureHists::Fill(fastjet::PseudoJet jet,double jet_radius, fastjet::PseudoJet matched_jet, double weight)
{

 	{
       
	  std::vector<fastjet::PseudoJet> SortedSubJets;
	  if(jet.has_user_info<HOTVRinfo>()) SortedSubJets=sorted_by_pt(jet.user_info<HOTVRinfo>().subjets());//Get subjets of hotvr jets
	 
	  
	  
	  
	  //calculate mmin
	  double mmin=0;
	  double m12=0;
	  double m01=0;
	  double m02=0;
	  if(SortedSubJets.size()>2){
	    m01 = 0;
	    m01=(SortedSubJets[0]+SortedSubJets[1]).m();
	    m02= 0;
	    m02=(SortedSubJets[0]+SortedSubJets[2]).m();
	    m12 = 0;
	    m12 = (SortedSubJets[1]+SortedSubJets[2]).m();
	    Hist("subjet01")->Fill(m01,weight);
	    Hist("subjet02")->Fill(m02,weight);
	    Hist("subjet12")->Fill(m12,weight);
	  }
	  mmin = std::min(m01,std::min(m02,m12));

	 	  
	  double jet_mass=jet.m();
	  double jet_mmin=mmin;
	  int jet_nsubjets=SortedSubJets.size();
	
	  double jet_subjet1pt=0;
	  double jet_subjet2pt=0;
	  if(SortedSubJets.size()>0) jet_subjet1pt=SortedSubJets.at(0).pt();//Fill subjet pT
	  if(SortedSubJets.size()>1) jet_subjet2pt=SortedSubJets.at(1).pt();//Fill subjet pT
      
     
	  Hist("jet_pT")->Fill(jet.pt(),weight); //Fill jet pT
      
	 //Fill number of subjets
	 Hist("nsubjets")->Fill(SortedSubJets.size(),weight);
	 //	 std::cout<<"subjet const "<<SortedSubJets.at(0).constituents().size()<<std::endl;

	
	 //Fill mmin
	 Hist("mmin")->Fill(mmin,weight);
	
	 //Fill subjet pT
	 for (unsigned int bk =0; bk<=3; ++bk)
	   {
	     if (SortedSubJets.size()> bk)
	       {
	  	 TString hname = TString::Format("pT_subjet%d", bk+1);
		 Hist(hname)->Fill(SortedSubJets.at(bk).pt(),weight);
		 hname = TString::Format("pT_subjetfrac%d", bk+1);
		 Hist(hname)->Fill(SortedSubJets.at(bk).pt()/jet.pt(),weight);
		 TString hname_ly = TString::Format("pT_subjet%d_ly", bk+1);
		 Hist(hname_ly)->Fill(SortedSubJets.at(bk).pt(),weight);
		 TString hname_eta = TString::Format("eta_subjet%d", bk+1);
		 Hist(hname_eta)->Fill(SortedSubJets.at(bk).eta(),weight);
		 TString hname_eta_ly = TString::Format("eta_subjet%d_ly", bk+1);
		 Hist(hname_eta_ly)->Fill(SortedSubJets.at(bk).eta(),weight);
		 TString hname_phi = TString::Format("phi_subjet%d", bk+1);
		 Hist(hname_phi)->Fill(SortedSubJets.at(bk).phi(),weight);
		 TString hname_phi_ly = TString::Format("phi_subjet%d_ly", bk+1);
		 Hist(hname_phi_ly)->Fill(SortedSubJets.at(bk).phi(),weight);
		 TString hname_mass= TString::Format("mass_subjet%d",bk+1);
		 Hist(hname_mass)->Fill(SortedSubJets.at(bk).m(),weight);

	       }
	   }
	 
	 //Fill jet mass
	 double mjet;
	 mjet=jet.m();
	 Hist("mjet")->Fill(mjet,weight);
	


	  std::cout<<"here3"<<std::endl;
	 //Fill pT fraction
	  if(SortedSubJets.size()>0) Hist("ptfraction1")->Fill( SortedSubJets.at(0).pt()/jet.pt(),weight);

	 
	 //Fill jet radius as function of pT
	  ((TH2D*)Hist("VariableR"))->Fill(jet.perp(),jet_radius*10,weight); 
       
	 //calculate nsubjettiness
	  std::cout<<"here2"<<std::endl;
	 if(jet.constituents().size()>0){
	   double tau1,tau2,tau3,tau4;
	   JetPropsPseudo jp(&jet);
	   jet_radius=1.5;
	   tau1 = jp.GetNsubjettiness(1, Njettiness::onepass_kt_axes, 1., jet_radius);
	   tau2 = jp.GetNsubjettiness(2, Njettiness::onepass_kt_axes, 1., jet_radius);
	   tau3 = jp.GetNsubjettiness(3, Njettiness::onepass_kt_axes, 1., jet_radius );
	   tau4 = jp.GetNsubjettiness(4, Njettiness::onepass_kt_axes, 1., jet_radius);
	   
	   //Fill nsubjettiness
	  
	   Hist("tau1")->Fill(tau1,weight);
	 
	   Hist("tau2")->Fill(tau2,weight);
	   	  
	   Hist("tau3")->Fill(tau3,weight);
	  
	   Hist("tau4")->Fill(tau4,weight);
	  
	   Hist("tau2tau1")->Fill(tau2/tau1,weight);
	   
	   Hist("tau3tau2")->Fill(tau3/tau2,weight);
	  
	 
	 //Fill more nominator hists
	
	 
	 }
	 
       
     }
 
} 
   
	  
	   
  
      
    




   

  

 




void SubstructureHists::Finish()
{
  // final calculations, like division and addition of certain histograms
  EventCalc* calc = EventCalc::Instance();
  bool IsRealData = calc->IsRealData();
  if (IsRealData){
    Hist("N_pv_perLumiBin")->Divide( Hist("N_pv_perLumiBin"), Hist("N_events_perLumiBin"));
    Hist( "N_pv_perLumiBin")->GetYaxis()->SetTitle("Events/Lumi");
  }

}

