#include "SFrameTools/include/EventCalc.h"
#include "include/SubstructureHists.h"
#include "include/TopFitCalc.h"
#include "NtupleWriter/include/JetProps.h"
#include <iostream>

using namespace std;

SubstructureHists::SubstructureHists(const char* name) : BaseHists(name)
{
  // named default constructor

}

SubstructureHists::~SubstructureHists()
{
  // default destructor, does nothing
}

void SubstructureHists::Init()
{
  // book all histograms here
  Book(TH1D("test","test",100,0,100));
  Book(TH1D("npv","number of primary vertices",51,0,50));
  Book( TH1D( "pT"," p_{T} topjets",100,0,2000));
  Book( TH1D( "pT_ly"," p_{T} topjets",100,0,2000));
  Book( TH1D( "eta"," #eta topjets",100,-3,3));
  Book( TH1F( "phi"," #phi topjets",100,-M_PI,M_PI));
  Book(TH1D("nsubjets","number subjets",11,0,10));
  Book(TH1D("subjet1pt","p_{T} leading subjet",70,0,1000));
  Book(TH1D("subjet1pt_ly","p_{T} leading subjet",70,0,1000));
  Book(TH1D("subjet2pt","p_{T} 2nd subjet",70,0,500));
  Book(TH1D("subjet2pt_ly","p_{T} 2nd subjet",70,0,500));
  Book(TH1D("subjet3pt","p_{T} 3rd subjet",70,0,500));
  Book(TH1D("subjet3pt_ly","p_{T} 3rd subjet",70,0,500));
  Book(TH1D("Jetmass","mass of CA-Jet",50,0,500));
  Book(TH1D("Subjetmass12","subjet mass_{12}",100,0,300));
  Book(TH1D("Subjetmass13","subjet mass_{13}",100,0,300));
  Book(TH1D("Subjetmass23","subjet_mass_{23}",100,0,300));
  Book(TH1D("subjet_mmin","subjet m_{min}",100,0,160));
  Book(TH1D("tau1","N-subjettiness #tau_{1}",50,0,1));
  Book(TH1D("tau2","N-subjettiness #tau_{2}",50,0,1));  
  Book(TH1D("tau3","N-subjettiness #tau_{3}",50,0,1));
  Book(TH1D("tau3tau2","N-subjettiness #tau_{3}/#tau_{2}",50,0,1));
  Book(TH1D("tau2tau1","N-subjettiness #tau_{2}/#tau_{1}",50,0,1));
  Book(TH1D("Jetmasspruned","mass of CA-Jet pruned",50,0,500));
  Book(TH1D("Subjetmass12pruned","subjet mass_{12p} pruned",100,0,300));
  Book(TH1D("Subjetmass13pruned","subjet mass_{13p} pruned",100,0,300));
  Book(TH1D("Subjetmass23pruned","subjet mass_{23p} pruned",100,0,300));
  Book(TH1D("tau1pruned","N-subjettiness #tau_{1p} pruned",50,0,1));
  Book(TH1D("tau2pruned","N-subjettiness #tau_{2p} pruned",50,0,1));  
  Book(TH1D("tau3pruned","N-subjettiness #tau_{3p} pruned",50,0,1));
  Book(TH1D("tau3tau2pruned","N-subjettiness #tau_{3p}/#tau_{2p} pruned",50,0,1));
  Book(TH1D("tau2tau1pruned","N-subjettiness #tau_{2p}/#tau_{1p} pruned",50,0,1));
  Book(TH1D("Qvolatility","Q-jets volatilty",30,0,0.6));
  Book(TH1D("HelAng12","Helicity angle subjet 1,2",35,0,3.5));
  Book(TH1D("HelAng13","Helicity angle subjet 1,3",35,0,3.5));
  Book(TH1D("HelAng23","Helicity angle subjet 2,3",35,0,3.5));
  Book( TH1F( "psi02_ly","Jetshapes #Psi (0.2)",70,0,1.4));
  Book( TH1F( "psi04_ly","Jetshapes #Psi (0.4)",70,0,1.4));
  Book( TH1F( "psi06_ly","Jetshapes #Psi (0.6)",70,0,1.4));
  Book( TH1F( "psi08_ly","Jetshapes #Psi (0.8)",70,0,1.4));
  Book( TH1F( "psi10_ly","Jetshapes #Psi (1.0)",70,0,1.4));
  Book( TH1F( "number_of_constituents","Number of constituents",100,0,400));
  Book( TH1F( "number_of_charged_constituents","Number of charged constituents",50,0,200));
  Book( TH1F( "jet_charge","charge of topjets",100,-50,50));
  Book( TH1F( "weighted_jet_charge","topjet weighted charge",55,-1,1));
  Book( TH1F( "weighted_jet_charge_02","topjet weighted charge k=0.2",100,-10,10));
  Book( TH1F( "weighted_jet_charge_04","topjet weighted charge k=0.4",105,-5,5));
  Book( TH1F( "weighted_jet_charge_06","topjet weighted charge k=0.6",105,-2,2));
  Book( TH1F( "weighted_jet_charge_08","topjet weighted charge k=0.8",105,-2,2));
  Book( TH1F( "first_jet_moment", "topjet first moment",100,0,1));
  Book( TH1F( "second_jet_moment", "topjet second moment",100,0,1));
  Book(TH1F("MVA","MVA values after selection",100,-2,2));
  Book(TH1F("MVA_ly","MVA values after selection",100,-2,2));
  /*  Book(TH1F("MVA","MVA values after selection",400,-2,2));
  Book(TH1F("MVA_ly","MVA values after selection",400,-2,2));
  Book(TH1F("MVA350","MVA values after selection",400,-2,2));
  Book(TH1F("MVA350_ly","MVA values after selection",400,-2,2));
  Book(TH1F("MVA2","MVA values before selection",400,-2,2));
  Book(TH1F("MVA2_ly","MVA values before selection",400,-2,2));
  Book( TH1F( "NbJets", "number of b-jets", 8, -0.5, 7.5 ) );
  Book( TH1F( "NbJets_ly", "number of b-jets", 8, -0.5, 7.5 ) );
  Book( TH1F( "psi02_ly","Jetshapes #Psi (0.2)",140,0,1.4));
  Book( TH1F( "psi04_ly","Jetshapes #Psi (0.4)",140,0,1.4));
  Book( TH1F( "psi06_ly","Jetshapes #Psi (0.6)",140,0,1.4));
  Book( TH1F( "psi08_ly","Jetshapes #Psi (0.8)",140,0,1.4));
  Book( TH1F( "psi10_ly","Jetshapes #Psi (1.0)",140,0,1.4));
  Book( TH1F( "number_of_constituents","Number of constituents",400,0,400));
  Book( TH1F( "number_of_charged_constituents","Number of charged constituents",200,0,200));
  Book( TH1F( "jet_charge","charge of topjets",100,-50,50));
  Book( TH1F( "weighted_jet_charge","topjet weighted charge",220,-2,2));
  Book( TH1F( "weighted_jet_charge_02","topjet weighted charge k=0.2",820,-10,10));
  Book( TH1F( "weighted_jet_charge_04","topjet weighted charge k=0.4",420,-5,5));
  Book( TH1F( "weighted_jet_charge_06","topjet weighted charge k=0.6",210,-2,2));
  Book( TH1F( "weighted_jet_charge_08","topjet weighted charge k=0.8",210,-2,2));
  Book( TH1F( "first_jet_moment", "topjet first moment",200,0,1));
  Book( TH1F( "second_jet_moment", "topjet second moment",200,0,1));	
  Book( TH1F( "NTopJets", "number of topjets", 8, -0.5, 7.5 ) );
  Book( TH1F( "NTopJets_ly", "number of topjets", 8, -0.5, 7.5 ) );
  Book( TH1F( "pT"," p_{T} topjets",100,0,2000));
  Book( TH1F( "pT_ly"," p_{T} topjets",100,0,2000));
  Book( TH1F( "pT_s1"," p_{T} topjets (selection)",50,0,2000));
  Book( TH1F( "pT_s1_ly"," p_{T} topjets (selection)",50,0,2000));
  Book( TH1F( "pT_s1_cms_tagged"," p_{T} topjets (cms tagged)",50,0,2000));
  Book( TH1F( "pT_s1_cms_tagged_ly"," p_{T} topjets (cms tagged)",50,0,2000));
   Book( TH1F( "pT_s1_arne_tagged"," p_{T} topjets (arne tagged)",50,0,2000));
  Book( TH1F( "pT_s1_arne_tagged_ly"," p_{T} topjets (arne tagged)",50,0,2000));
  Book( TH1F( "eta"," #eta topjets",100,-3,3));
  Book( TH1F( "eta_ly"," #eta topjets",100,-3,3));
  Book( TH1F( "phi"," #phi topjets",100,-M_PI,M_PI));
  Book( TH1F( "phi_ly"," #phi topjets",100,-M_PI,M_PI));
  Book( TH1F( "MJet", "m_{jet}", 100,0,400 ) );
  Book( TH1F( "MJet_ly", "m_{jet}", 100,0,400 ) );
  Book( TH1F( "Mmin", "m_{min}", 100,0,160 ) );
  Book( TH1F( "Mmin_ly", "m_{min}", 100,0,160 ) );
  Book( TH1F( "NSubjets", "number of subjets", 6,-0.5,5.5) );
  Book( TH1F( "NSubjets_ly", "number of subjets", 6,-0.5,5.5 ) );
  Book( TH1F( "pT_1"," p_{T} leading topjet",100,0,2000));
  Book( TH1F( "pT_1_ly"," p_{T} leading topjet",100,0,2000));
  Book( TH1F( "pT_2","p_{T} 2nd topjet",100,0,2000));
  Book( TH1F( "pT_2_ly","p_{T} 2nd topjet",100,0,2000));
  Book( TH1F( "pT_3","p_{T} 3rd topjet",100,0,1000));
  Book( TH1F( "pT_3_ly","p_{T} 3rd topjet",100,0,1000));
  Book( TH1F( "pT_4","p_{T} 4th topjet",100,0,800));
  Book( TH1F( "pT_4_ly","p_{T} 4th topjet",100,0,800));
  Book( TH1F( "eta_1","#eta leading topjet",100,-3,3));
  Book( TH1F( "eta_1_ly","#eta leading topjet",100,-3,3));
  Book( TH1F( "eta_2","#eta 2nd topjet",100,-3,3));
  Book( TH1F( "eta_2_ly","#eta 2nd topjet",100,-3,3));
  Book( TH1F( "eta_3","#eta 3rd topjet",100,-3,3));
  Book( TH1F( "eta_3_ly","#eta 3rd topjet",100,-3,3));
  Book( TH1F( "eta_4","#eta 4th topjet",100,-3,3));
  Book( TH1F( "eta_4_ly","#eta 4th topjet",100,-3,3));
  Book( TH1F( "phi_1","#phi leading topjet",100,-M_PI,M_PI));
  Book( TH1F( "phi_1_ly","#phi leading topjet",100,-M_PI,M_PI));
  Book( TH1F( "phi_2","#phi 2nd topjet",100,-M_PI,M_PI));
  Book( TH1F( "phi_2_ly","#phi 2nd topjet",100,-M_PI,M_PI));
  Book( TH1F( "phi_3","#phi 3rd topjet",100,-M_PI,M_PI));
  Book( TH1F( "phi_3_ly","#phi 3rd topjet",100,-M_PI,M_PI));
  Book( TH1F( "phi_4","#phi 4th topjet",100,-M_PI,M_PI));
  Book( TH1F( "phi_4_ly","#phi 4th topjet",100,-M_PI,M_PI));
  Book( TH1F( "MJet_1", "m_{jet} leading topjet", 100,0,400 ) );
  Book( TH1F( "MJet_1_ly", "m_{jet} leading topjet", 100,0,400 ) );
  Book( TH1F( "MJet_2", "m_{jet} 2nd topjet", 100,0,400 ) );
  Book( TH1F( "MJet_2_ly", "m_{jet} 2nd topjet", 100,0,400 ) );
  Book( TH1F( "MJet_3", "m_{jet} 3rd topjet", 100,0,400 ) );
  Book( TH1F( "MJet_3_ly", "m_{jet} 3rd topjet", 100,0,400 ) );
  Book( TH1F( "MJet_4", "m_{jet} 4th topjet", 100,0,400 ) );
  Book( TH1F( "MJet_4_ly", "m_{jet} 4th topjet", 100,0,400 ) );
  Book( TH1F( "Mmin_1", "m_{min} leading topjet", 100,0,160 ) );
  Book( TH1F( "Mmin_1_ly", "m_{min}, leading topjet", 100,0,160 ) );
  Book( TH1F( "Mmin_2", "m_{min} 2nd topjet", 100,0,160 ) );
  Book( TH1F( "Mmin_2_ly", "m_{min} 2nd topjet", 100,0,160 ) );
  Book( TH1F( "Mmin_3", "m_{min} 3rd topjet", 100,0,160 ) );
  Book( TH1F( "Mmin_3_ly", "m_{min} 3rd topjet", 100,0,160 ) );
  Book( TH1F( "Mmin_4", "m_{min} 4th topjet", 100,0,160 ) );
  Book( TH1F( "Mmin_4_ly", "m_{min} 4th topjet", 100,0,160 ) );
  Book( TH1F( "NSubjets_1", "number of subjets leading topjet", 6,-0.5,5.5) );
  Book( TH1F( "NSubjets_1_ly", "number of subjets leading topjet", 6,-0.5,5.5 ) );
  Book( TH1F( "NSubjets_2", "number of subjets 2nd topjet", 6,-0.5,5.5 ) );
  Book( TH1F( "NSubjets_2_ly", "number of subjets 2nd topjet", 6,-0.5,5.5 ) );
  Book( TH1F( "NSubjets_3", "number of subjets 3rd topjet", 6,-0.5,5.5 ) );
  Book( TH1F( "NSubjets_3_ly", "number of subjets 3rd topjet", 6,-0.5,5.5 ) );
  Book( TH1F( "NSubjets_4", "number of subjets 4th topjet", 6,-0.5,5.5 ) ); 
  Book( TH1F( "NSubjets_4_ly", "number of subjets 4th topjet", 6,-0.5,5.5 ) ); 

  Book( TH1F( "QjetsVol_1", "Qjets volatility leading topjet", 50, 0.0, 2.0) );
  Book( TH1F( "QjetsVol_1_ly", "Qjets volatility leading topjet", 50, 0.0, 2.0 ) );

  Book( TH1F( "Nsubjettiness1_1", "#tau_{1} leading topjet", 50, 0.0, 1.0) );
  Book( TH1F( "Nsubjettiness1_1_ly", "#tau_{1} leading topjet", 50, 0.0, 1.0 ) );

  Book( TH1F( "Nsubjettiness2_1", "#tau_{2} leading topjet", 50, 0.0, 1.0) );
  Book( TH1F( "Nsubjettiness2_1_ly", "#tau_{2} leading topjet", 50, 0.0, 1.0 ) );
  
  Book( TH1F( "Nsubjettiness3_1", "#tau_{3} leading topjet", 50, 0.0, 1.0) );
  Book( TH1F( "Nsubjettiness3_1_ly", "#tau_{3} leading topjet", 50, 0.0, 1.0 ) );

  Book( TH1F( "Nsubjettiness3_2_1", "#tau_{2}/#tau_{3} leading topjet", 50, 0.0, 1.0) );
  Book( TH1F( "Nsubjettiness3_2_1_ly", "#tau_{2}/#tau_{3} leading topjet", 50, 0.0, 1.0 ) );
  Book( TH2F( "MVA_pt","MVA vs pt",200,-1,1,50,0,2000));
  Book( TH1F( "MVA_eff", "efficiency of MVA",200,-1,1));*/
 

}


void SubstructureHists::Fill2(TopJet topjet, double mva_value)
{
  EventCalc* calc = EventCalc::Instance();
  double weight = calc->GetWeight();
  double npv = calc->GetPrimaryVertices()->size();
   Hist("MVA")->Fill(mva_value,weight);
   Hist("MVA_ly")->Fill(mva_value,weight);
  Hist("npv")->Fill(npv,weight);
  
  int Nsubjets=topjet.numberOfDaughters();
  Hist("nsubjets")->Fill(Nsubjets,weight);
  LorentzVector allsubjets(0,0,0,0);
  
  for(int j=0; j<topjet.numberOfDaughters(); ++j){
    allsubjets += topjet.subjets()[j].v4();
  }
  double Topjetmass;
  if(!allsubjets.isTimelike()){
    Topjetmass=0;
  } else {
    Topjetmass = allsubjets.M();
  }
  Hist("Jetmass")->Fill(Topjetmass,weight);
  
  
  //subjet mass
  double Subjet12mass;
  double Subjet13mass;
  double Subjet23mass;
  std::vector<Particle> subjet = topjet.subjets();
  sort(subjet.begin(),subjet.end(),HigherPt());
 if(Nsubjets>=2){
    Subjet12mass=(subjet[0].v4()+subjet[1].v4()).mass();
     }
  if(Nsubjets>=3){
    Subjet12mass=(subjet[0].v4()+subjet[1].v4()).mass();
    Subjet13mass=(subjet[0].v4()+subjet[2].v4()).mass();
    Subjet23mass=(subjet[1].v4()+subjet[2].v4()).mass();
  }
  if(Subjet12mass!=0)  Hist("Subjetmass12")->Fill(Subjet12mass,weight);
  if(Subjet13mass!=0) Hist("Subjetmass13")->Fill(Subjet13mass,weight);
  if(Subjet23mass!=0) Hist("Subjetmass23")->Fill(Subjet23mass,weight);
  
   //nsubjetiness
   JetProps jp(&topjet, calc->GetPFParticles() );
  double tau1 = jp.GetNsubjettiness(1, Njettiness::onepass_kt_axes, 1., 0.8);
  double tau2 = jp.GetNsubjettiness(2, Njettiness::onepass_kt_axes, 1., 0.8);
  double tau3 = jp.GetNsubjettiness(3, Njettiness::onepass_kt_axes, 1., 0.8);
  Hist("tau1")->Fill(tau1,weight);
  Hist("tau2")->Fill(tau2,weight);
  Hist("tau3")->Fill(tau3,weight);
  //pruned subjetiness
  double tau1pruned = jp.GetPrunedNsubjettiness(1, Njettiness::onepass_kt_axes, 1., 0.8);
  double tau2pruned = jp.GetPrunedNsubjettiness(2, Njettiness::onepass_kt_axes, 1., 0.8);
  double tau3pruned = jp.GetPrunedNsubjettiness(3, Njettiness::onepass_kt_axes, 1., 0.8);
  Hist("tau1pruned")->Fill(tau1pruned,weight);
  Hist("tau2pruned")->Fill(tau2pruned,weight);
  Hist("tau3pruned")->Fill(tau3pruned,weight);
  
  //ratio subjetiness
  Hist("tau3tau2")->Fill(tau3/tau2,weight);
  Hist("tau2tau1")->Fill(tau2/tau1,weight);
  Hist("tau3tau2pruned")->Fill(tau3pruned/tau2pruned,weight);
  Hist("tau2tau1pruned")->Fill(tau2pruned/tau1pruned,weight);
  
  double Q_volatility = topjet.qjets_volatility();
  Hist("Qvolatility")->Fill(Q_volatility,weight);
 
  
 // calculate pruned masses
  
  double Topjetmass_pruned;
    double Subjet12mass_pruned;
    double Subjet13mass_pruned;
    double Subjet23mass_pruned;
  std::vector<fastjet::PseudoJet> jets = jp.GetFastJet(2.0);   // something large to make sure jet is inside radius
  if(jets.empty()){
      m_logger << WARNING << "TMVATreeFiller::FillTopJetProperties: no jet found!" << SLogger::endmsg; 
  }
  else{
    fastjet::PseudoJet pjet = jp.GetPrunedJet(jets[0]);
    std::vector<fastjet::PseudoJet> prunedsubjets;
    if (pjet.constituents().size()>=2){
        prunedsubjets = pjet.exclusive_subjets(2);
    }
    if (pjet.constituents().size()>=3){
        prunedsubjets = pjet.exclusive_subjets(3);
    }

    unsigned int pnsubs = prunedsubjets.size();

    fastjet::PseudoJet psubjets(0,0,0,0);
    
    for(unsigned int j=0; j<pnsubs; ++j){
        psubjets += pjet.pieces()[j];
    }
    Topjetmass_pruned = psubjets.m();
    if (pnsubs>=2) {
        Subjet12mass_pruned = (prunedsubjets[0]+prunedsubjets[1]).m();
    }

    if (pnsubs>=3) {
        Subjet13mass_pruned = (prunedsubjets[0]+prunedsubjets[2]).m();
        Subjet23mass_pruned = (prunedsubjets[1]+prunedsubjets[2]).m();
    }
  }
  Hist("Jetmasspruned")->Fill(Topjetmass_pruned,weight);
  if(Subjet12mass_pruned!=0) Hist("Subjetmass12pruned")->Fill(Subjet12mass_pruned,weight);
   if(Subjet13mass_pruned!=0) Hist("Subjetmass13pruned")->Fill(Subjet13mass_pruned,weight);
   if(Subjet23mass_pruned!=0) Hist("Subjetmass23pruned")->Fill(Subjet23mass_pruned,weight);


  
   
  //HelicityAngles && subjet pt
  double subjet1pt;
  double subjet2pt;
  double subjet3pt;
  double HelAng12;
  double HelAng13;
  double HelAng23;
  if(Nsubjets>0) subjet1pt=sqrt(pow(subjet[0].v4().Px(),2)+pow(subjet[0].v4().Py(),2));
  if(Nsubjets>1){
    TLorentzVector subjet1(subjet[0].v4().Px(),subjet[0].v4().Py(),subjet[0].v4().Pz(),subjet[0].v4().E());
    TLorentzVector subjet2(subjet[1].v4().Px(),subjet[1].v4().Py(),subjet[1].v4().Pz(),subjet[1].v4().E());
    subjet1pt=sqrt(pow(subjet[0].v4().Px(),2)+pow(subjet[0].v4().Py(),2));
    subjet2pt=sqrt(pow(subjet[1].v4().Px(),2)+pow(subjet[1].v4().Py(),2));
    HelAng12= HelicityAngle(subjet1,subjet2);
  }
  if(Nsubjets>2){
    TLorentzVector subjet1(subjet[0].v4().Px(),subjet[0].v4().Py(),subjet[0].v4().Pz(),subjet[0].v4().E());
    TLorentzVector subjet2(subjet[1].v4().Px(),subjet[1].v4().Py(),subjet[1].v4().Pz(),subjet[1].v4().E());
    TLorentzVector subjet3(subjet[2].v4().Px(),subjet[2].v4().Py(),subjet[2].v4().Pz(),subjet[2].v4().E());
    subjet1pt=sqrt(pow(subjet[0].v4().Px(),2)+pow(subjet[0].v4().Py(),2));
    subjet2pt=sqrt(pow(subjet[1].v4().Px(),2)+pow(subjet[1].v4().Py(),2));
    subjet3pt=sqrt(pow(subjet[2].v4().Px(),2)+pow(subjet[2].v4().Py(),2));
    HelAng12= HelicityAngle(subjet1,subjet2);
    HelAng13= HelicityAngle(subjet1,subjet3);
    HelAng23= HelicityAngle(subjet2,subjet3);
    }
  if(subjet1pt!=0){
  Hist("subjet1pt")->Fill(subjet1pt,weight);
  Hist("subjet1pt_ly")->Fill(subjet1pt,weight);}
   if(subjet2pt!=0){
  Hist("subjet2pt")->Fill(subjet2pt,weight);
  Hist("subjet2pt_ly")->Fill(subjet2pt,weight);}
    if(subjet3pt!=0){
  Hist("subjet3pt")->Fill(subjet3pt,weight);
  Hist("subjet3pt_ly")->Fill(subjet3pt,weight);}

    if(HelAng12!=0) Hist("HelAng12")->Fill(HelAng12,weight);
  if(HelAng13!=0) Hist("HelAng13")->Fill(HelAng13,weight);
 if(HelAng23!=0) Hist("HelAng23")->Fill(HelAng23,weight);

 
  //calculate Jetshapes
  double m_psi_02 = calc->IntegratedJetShape( &topjet, 0.2, 0.0 , e_CA8);
  double m_psi_04 = calc->IntegratedJetShape( &topjet, 0.4, 0.0 , e_CA8);
  double m_psi_06 = calc->IntegratedJetShape( &topjet, 0.6, 0.0 , e_CA8);
  double m_psi_08 = calc->IntegratedJetShape( &topjet, 0.8, 0.0 , e_CA8);
  double m_psi_10 = calc->IntegratedJetShape( &topjet, 1.0, 0.0 , e_CA8);
  Hist("psi02_ly")->Fill(m_psi_02,weight);
  Hist("psi04_ly")->Fill(m_psi_04,weight);
  Hist("psi06_ly")->Fill(m_psi_06,weight);
  Hist("psi08_ly")->Fill(m_psi_08,weight);
  Hist("psi10_ly")->Fill(m_psi_10,weight);
       
  //calculate jet properties
  std::vector<PFParticle> jetconsts = calc->GetJetPFParticles(&topjet);
  int m_number_of_constituents = jetconsts.size();
  int m_number_of_charged_constituents = 0;
  double m_jet_charge = calc->JetCharge(&topjet);
  Hist("number_of_constituents")->Fill(m_number_of_constituents,weight);
  Hist("jet_charge")->Fill(m_jet_charge,weight);
  for(int j=0; j< m_number_of_constituents; ++j){
    if( fabs(jetconsts[j].charge())>0.01) {
      m_number_of_charged_constituents ++;
    }
  }
  Hist("number_of_charged_constituents")->Fill(m_number_of_charged_constituents,weight);
  double m_weighted_jet_charge = calc->EnergyWeightedJetCharge(&topjet);
  double m_weighted_jet_charge_02 = calc->EnergyWeightedJetCharge(&topjet, 0.2);
  double m_weighted_jet_charge_04 = calc->EnergyWeightedJetCharge(&topjet, 0.4);
  double m_weighted_jet_charge_06 = calc->EnergyWeightedJetCharge(&topjet, 0.6);
  double m_weighted_jet_charge_08 = calc->EnergyWeightedJetCharge(&topjet, 0.8);
  Hist("weighted_jet_charge")->Fill(m_weighted_jet_charge,weight);
  Hist("weighted_jet_charge_02")->Fill(m_weighted_jet_charge_02,weight);
  Hist("weighted_jet_charge_04")->Fill(m_weighted_jet_charge_04,weight);
  Hist("weighted_jet_charge_06")->Fill(m_weighted_jet_charge_06,weight);
  Hist("weighted_jet_charge_08")->Fill(m_weighted_jet_charge_08,weight);
  double m_first_jet_moment = calc->JetMoment(&topjet,1);
  double m_second_jet_moment = calc->JetMoment(&topjet,2);
  Hist("first_jet_moment")->Fill(m_first_jet_moment,weight);
  Hist("second_jet_moment")->Fill(m_second_jet_moment,weight);
   
  //pt topjet    
  Hist("pT") -> Fill(topjet.pt(), weight);
  Hist("pT_ly") -> Fill(topjet.pt(), weight);
  Hist("eta") -> Fill(topjet.eta(),weight);
  Hist("phi") -> Fill(topjet.phi(),weight);
  double mmin=0;
  double mjet=0;
  int nsubjets=0;
  TopTag(topjet,mjet,nsubjets,mmin);
  if(mmin!=0) Hist("subjet_mmin")->Fill(mmin,weight);
 
 
}

void SubstructureHists::Fill()
{
}
 
