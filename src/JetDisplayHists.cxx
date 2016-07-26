#include <time.h>
#include "include/JetDisplayHists.h"
#include "SFrameTools/include/EventCalc.h"
//#include "include/ZprimeFullHadTools.h"
#include "NtupleWriter/include/JetProps.h"
//#include "include/TopTagfunctions.h"
 #include "fastjet/PseudoJet.hh"
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
JetDisplayHists::JetDisplayHists(const char* name) : BaseHists(name)
{
  //  std::cout<<name<<std::endl;
  // m_eventnumber=eventnumber;
  // named default constructor
   
}

JetDisplayHists::~JetDisplayHists()
{
  // default destructor, does nothing
}

void JetDisplayHists::Init()
{
  // book all histograms here
 //jet display
  int bins=50;
  for(int o=0;o<10;o++) {
    TString hname3 = TString::Format("JetDisplay_jet%i_subjetmass",o);
    Book(TH2F(hname3,hname3,10,0,9,300,0,300));
    hname3 = TString::Format("JetDisplay_jet%i_subjetpT",o);
    Book(TH2F(hname3,hname3,10,0,9,300,0,3000));
    for(int p=0;p<20;p++) {
      TString hname2 = TString::Format("JetDisplay_jet%i_subjet%i", o,p);
      Book(TH2F(hname2,"Jet event display",bins,-PI,PI,bins,-PI,PI));
    }
  }
  
  Book(TH2F("JetDisplay_beam","Jet event display",bins,-PI,PI,bins,-PI,PI));
  Book(TH2F("JetDisplay_radiation","Jet event display",bins,-PI,PI,bins,-PI,PI));
  Book(TH2F("JetDisplay_rejected_subjets","Jet event display",bins,-PI,PI,bins,-PI,PI));
  Book(TH2F("JetDisplay","Jet event display",bins,-PI,PI,bins,-PI,PI));
  Book(TH2F("JetDisplay1","Jet event display",bins,-PI,PI,bins,-PI,PI));
  Book(TH2F("JetDisplay_fatjet","Jet event display",bins,-PI,PI,bins,-PI,PI));
  Book(TH2F("JetDisplay_subjets","Jet event display",bins,-PI,PI,bins,-PI,PI));
  Book(TH2F("JetDisplay_pf0","Jet event display",bins,-PI,PI,bins,-PI,PI));
  Book(TH2F("JetDisplay_pf1","Jet event display",bins,-PI,PI,bins,-PI,PI));
  Book(TH2F("JetDisplay_pf2","Jet event display",bins,-PI,PI,bins,-PI,PI));
  Book(TH2F("JetDisplay_pf_all","Jet event display",bins,-PI,PI,bins,-PI,PI));
  Book(TH2F("JetDisplay_top","Jet event display",bins,-PI,PI,bins,-PI,PI));
  Book(TH2F("JetDisplay_gluon","Jet event display",bins,-PI,PI,bins,-PI,PI));
  Book(TH2F("JetDisplay_light","Jet event display",bins,-PI,PI,bins,-PI,PI));
  Book(TH2F("JetDisplay_decay","Jet event display",bins,-PI,PI,bins,-PI,PI));
 

  Book(TH2F("Jetmass","jetmass",10,0,9,300,0,300));
  Book(TH2F("JetpT","jet p_{T}",10,0,9,300,0,3000));
  

  IR_Saftey= new Infrared_Saftey();
  clustering= new Clustering();
}

void JetDisplayHists::SetIdVersion(TString s)
{
  idVersion=s;
}

void  JetDisplayHists::Fill(){
}

bool JetDisplayHists::FillEvent(std::vector<fastjet::PseudoJet> jets, std::vector<fastjet::PseudoJet> parts,std::vector<fastjet::PseudoJet> rejected_jets, std::vector<fastjet::PseudoJet> soft_jets, std::vector<fastjet::PseudoJet> rejected_subjets)
{


  EventCalc* calc = EventCalc::Instance();
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



	//cluster with HOTVR
    
   for(int y=0;y<rejected_jets.size();y++) for(int k=0; k<rejected_jets[y].constituents().size();k++)  ((TH2D*)Hist("JetDisplay_beam"))->Fill(rejected_jets[y].constituents().at(k).phi_std(),rejected_jets[y].constituents().at(k).eta(), rejected_jets[y].constituents().at(k).perp()     );
   for(int y=0;y<soft_jets.size();y++) for(int k=0;k<soft_jets[y].constituents().size();k++)   ((TH2D*)Hist("JetDisplay_radiation"))->Fill(soft_jets[y].constituents().at(k).phi_std(),soft_jets[y].constituents().at(k).eta(),soft_jets[y].constituents().at(k).perp());
    for(int y=0;y<rejected_subjets.size();y++) for(int k=0;k<rejected_subjets[y].constituents().size();k++)   ((TH2D*)Hist("JetDisplay_rejected_subjets"))->Fill(rejected_subjets[y].constituents().at(k).phi_std(),rejected_subjets[y].constituents().at(k).eta(),rejected_subjets[y].constituents().at(k).perp());

    
   for(int o=0;o<jets.size();o++) {
     ((TH2D*)Hist("JetDisplay_fatjet"))->Fill(jets[o].phi_std(),jets[o].eta(),jets[o].perp());
     double mass=sqrt(jets[o].m2());
    
    
     ((TH2D*)Hist("Jetmass"))->Fill(o,mass);
     ((TH2D*)Hist("JetpT"))->Fill(o,jets[o].pt());
     if(jets[o].has_user_info<HOTVRinfo>()) { 
	 std::vector<fastjet::PseudoJet> SortedSubJets=jets[o].user_info<HOTVRinfo>().subjets();
	
	 for(int p=0;p<SortedSubJets.size();p++) {
	   TString hname3 = TString::Format("JetDisplay_jet%i_subjetmass",o);
	   Hist(hname3)->Fill(p,SortedSubJets.at(p).m());
	   hname3 = TString::Format("JetDisplay_jet%i_subjetpT",o);
	   Hist(hname3)->Fill(p,SortedSubJets.at(p).pt());
	   TString hname2 = TString::Format("JetDisplay_jet%i_subjet%i", o,p);
	   
	   for(int h=0;h<SortedSubJets.at(p).constituents().size();h++) ((TH2D*)Hist(hname2))->Fill(SortedSubJets.at(p).constituents().at(h).phi_std(),SortedSubJets.at(p).constituents().at(h).eta(),SortedSubJets.at(p).constituents().at(h).perp() );
	   
	 }
       } 
     for(int h=0;h<jets[o].constituents().size();h++) {
       if(o==0)  ((TH2D*)Hist("JetDisplay_pf0"))->Fill(jets[o].constituents().at(h).phi_std(),jets[o].constituents().at(h).eta(),jets[o].constituents().at(h).perp());
       if(o==1) ((TH2D*)Hist("JetDisplay_pf1"))->Fill(jets[o].constituents().at(h).phi_std(),jets[o].constituents().at(h).eta(),jets[o].constituents().at(h).perp());
       if(o==2) ((TH2D*)Hist("JetDisplay_pf2"))->Fill(jets[o].constituents().at(h).phi_std(),jets[o].constituents().at(h).eta(),jets[o].constituents().at(h).perp());
       
     }
   }
   for(int j=0; j< Had_Tops.size();j++) ((TH2D*)Hist("JetDisplay_top"))->Fill(Had_Tops[j].phi(),Had_Tops[j].eta(),Had_Tops[j].pt());
   for(int o=0;o<decay_products.size();o++) ((TH2D*)Hist("JetDisplay_decay"))->Fill(decay_products[o].phi(),decay_products[o].eta(),decay_products[o].pt());
   for(int o=0;o<parts.size();o++) ((TH2D*)Hist("JetDisplay_pf_all"))->Fill(parts[o].phi_std(),parts[o].eta(),parts[o].pt());
   



   return true;

}
    
      
         
	  
    
    
     
     
  




void JetDisplayHists::Finish()
{
  // final calculations, like division and addition of certain histograms
 

}

