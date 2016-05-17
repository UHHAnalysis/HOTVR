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

 
 
  Book( TH1F( "pT_s1"," p_{T} topjets (selection)",75,0,3000));
  Book( TH1F( "pT_s1_all"," p_{T} topjets (selection)",75,0,3000));
  Book( TH1F( "eta_s1_all"," eta topjets (selection)",20,-3,3));
  Book( TH1F( "pT_s1_mis"," p_{T} topjets (selection)",75,0,3000));
  Book( TH1F( "pT_s1400"," p_{T} topjets (selection)",75,0,3000));
  Book( TH1F( "pT_s1_cms_tagged400"," p_{T} topjets (cms tagged)",75,0,3000));
  Book( TH1F( "pT_s1_cms_nsub_tagged"," p_{T}  nsubjettiness (cms tagged)",75,0,3000));
  Book( TH1F( "pT_s1_hotvr_tagged"," p_{T} topjets (hotvr tagged)",75,0,3000));
  Book( TH1F( "pT_s1_hotvr2_tagged"," p_{T} topjets (hotvr tagged)",75,0,3000));
  Book( TH1F( "pT_s1_hotvr3_tagged"," p_{T} topjets (hotvr tagged)",75,0,3000));
  Book( TH1F( "pT_s1_hotvr4_tagged"," p_{T} topjets sd mass (hotvr tagged)",75,0,3000));
  Book( TH1F( "pT_s1_hotvr5_tagged"," p_{T} topjets nsubjettiness (hotvr tagged)",75,0,3000));
  Book( TH1F( "pT_s1_hotvr6_tagged"," p_{T} topjets w mass (hotvr tagged)",75,0,3000));
  Book( TH1F( "eta_s1_hotvr6_tagged"," eta topjets w mass (hotvr tagged)",20,-3,3));
  
  
  
  

  
  
  Book( TH1D( "weight", ";weight;Events", 400, 0., 1000.));
  
  
  
  //Variable R

  Book(TH2F("VariableR","VariableR",1000,0,2000,20,0,40));
 
 
 
  Book(TH1D("mmin","massjump mmin",50,0,100));
 

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

  Book( TH2D("nsub_eff","eff nsub",75,0,3000,100,-0.005,0.995));
  Book( TH2D("nsub_eff_norm","eff nsub",75,0,3000,100,-0.005,0.995));
 
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

	  std::cout<<"mmin1 "<<mmin<<std::endl;
	  std::cout<<"mmin2 "<<jet.user_info<HOTVRinfo>().mmin()<<std::endl;

	  
	  double jet_mass=jet.m();
	  double jet_mmin=mmin;
	  int jet_nsubjets=SortedSubJets.size();
	  double jet_ptfraction=SortedSubJets.at(0).pt()/jet.pt();
	  double jet_subjet1pt=0;
	  double jet_subjet2pt=0;
	  if(SortedSubJets.size()>0) jet_subjet1pt=SortedSubJets.at(0).pt();//Fill subjet pT
	  if(SortedSubJets.size()>1) jet_subjet2pt=SortedSubJets.at(1).pt();//Fill subjet pT
      
     

      
	 //Fill number of subjets
	 Hist("nsubjets")->Fill(SortedSubJets.size(),weight);
	 std::cout<<"nsubjets1 "<<SortedSubJets.size()<<std::endl;
	  std::cout<<"nsubjets2 "<<jet.user_info<HOTVRinfo>().nsubjets()<<std::endl;

     
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

	       }
	   }
	 
	 //Fill jet mass
	 double mjet;
	 mjet=jet.m();
	 Hist("mjet")->Fill(mjet,weight);
	


	 
	 //Fill pT fraction
	 Hist("ptfraction1")->Fill( SortedSubJets.at(0).pt()/jet.pt(),weight);
	 std::cout<<"ptfraction1 "<< SortedSubJets.at(0).pt()/jet.pt()<<" "<< SortedSubJets.at(0).pt()/jet.phi()<<std::endl;
	 std::cout<<"ptfraction2 "<<jet.user_info<HOTVRinfo>().ptfraction(0)<<std::endl;
	 
	 //Fill nominator hists
	 // if(SortedSubJets.size()>2)  Hist("pT_s1_hotvr_tagged")->Fill(matched_jet.pt(),weight);
	 /* if(hotvr_jets[i].m()>mtopLow &&  hotvr_jets[i].m()<mtopHigh)  Hist("pT_s1_hotvr2_tagged")->Fill(matched_jet.pt(),weight);
	 if(SortedSubJets.size()>2 && hotvr_jets[i].m()>mtopLow &&  hotvr_jets[i].m()<mtopHigh)  Hist("pT_s1_hotvr3_tagged")->Fill(matched_jet.pt(),weight);
	 if(SortedSubJets.size()>2 && hotvr_jets[i].m()>mtopLow &&  hotvr_jets[i].m()<mtopHigh    &&SortedSubJets.at(0).pt()>30 &&SortedSubJets.at(1).pt()>30&&  SortedSubJets.at(0).pt()/hotvr_jets[i].pt()<0.8 )  Hist("pT_s1_hotvr6_tagged")->Fill(matched_jet.pt(),weight);
	 if(SortedSubJets.size()>2 && hotvr_jets[i].m()>mtopLow &&  hotvr_jets[i].m()<mtopHigh && mmin>50 && hotvr_jets.size()>2)  Hist("pT_s1_hotvr4_tagged")->Fill(matched_jet.pt(),weight);
	 */
	 //Fill jet radius as function of pT
	  ((TH2D*)Hist("VariableR"))->Fill(jet.perp(),jet_radius*10,weight); 
       
	 //calculate nsubjettiness
	
	 if(jet.constituents().size()>0){
	   double tau1,tau2,tau3,tau4;
	   JetPropsPseudo jp(&jet);
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
	
	   /*if(SortedSubJets.size()>2 &&  hotvr_jets[i].m()>mtopLow &&   hotvr_jets[i].m()<mtopHigh &&tau3/tau2<0.7  &&SortedSubJets.at(0).pt()>30 &&SortedSubJets.at(1).pt()>30&&  SortedSubJets.at(0).pt()/hotvr_jets[i].pt()<0.8 &&mmin>50)  Hist("pT_s1_hotvr5_tagged")->Fill(matched_jet.pt(),weight);
	 if(SortedSubJets.size()>2 &&  hotvr_jets[i].m()>mtopLow &&   hotvr_jets[i].m()<mtopHigh   &&SortedSubJets.at(0).pt()>30 &&SortedSubJets.at(1).pt()>30&&  SortedSubJets.at(0).pt()/hotvr_jets[i].pt()<0.8 &&mmin>50) for(int roc=0;roc<100;roc++){
	     if(tau3/tau2<roc/100.){
	       ((TH2D*)Hist("nsub_eff"))->Fill(matched_jet.pt(),roc/100.,weight);
	    
	     }}*/
	 }
	 
       
     }
 
} // std::cout<<denominator_jets.size()<<std::endl;
   /*for(int y=0;y<denominator_jets.size();y++){
	Hist("pT_s1_all")->Fill(denominator_jets.at(y).pt(),weight); //fill denominator hists for efficiecies
	for(int roc=0;roc<100;roc++) {((TH2D*)Hist("nsub_eff_norm"))->Fill(denominator_jets.at(y).pt(),roc/100.,weight);}//Fill denominator hist for nsubjettiness scan

	}*/
   
	  
	   
  
      
    




   

  

 




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

