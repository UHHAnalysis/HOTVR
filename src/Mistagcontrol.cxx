#include "SFrameTools/include/EventCalc.h"
#include "include/Mistagcontrol.h"
#include "include/TopFitCalc.h"
#include <iostream>

using namespace std;

Mistagcontrol::Mistagcontrol(const char* name) : BaseHists(name)
{
  // named default constructor

}

Mistagcontrol::~Mistagcontrol()
{
  // default destructor, does nothing
}

void Mistagcontrol::Init()
{
  // book all histograms here
  Book(TH1F("MVA","MVA values after selection",400,-2,2));
  Book(TH1F("MVA_ly","MVA values after selection",400,-2,2));
  Book(TH1F("MVA350","MVA values after selection",400,-2,2));
  Book(TH1F("MVA350_ly","MVA values after selection",400,-2,2));
  Book(TH1F("MVA2","MVA values before selection",400,-2,2));
  Book(TH1F("MVA2_ly","MVA values before selection",400,-2,2));
  Book(TH1F("SD","chi (microjets) shower deconstruction",21,-10,10));
  Book(TH1F("SD2","chi(subjets) shower deconstruction",21,-10,10));
  Book( TH1F( "NMicrojets", "number of microjets", 11,-0.5,10.5) ); 
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
  Book( TH1F( "pT_s2"," p_{T} CA15 (selection)",50,0,2000));
  Book( TH1F( "pT_s2_ly"," p_{T} CA15 (selection)",50,0,2000));
  Book( TH1F( "pT_s1_cms_tagged"," p_{T} topjets (cms tagged)",50,0,2000));
  Book( TH1F( "pT_s1_cms_tagged_ly"," p_{T} topjets (cms tagged)",50,0,2000));
   Book( TH1F( "pT_s1_hep_tagged"," p_{T} topjets (hep tagged)",50,0,2000));
  Book( TH1F( "pT_s1_hep_tagged_ly"," p_{T} topjets (hep tagged)",50,0,2000));
   Book( TH1F( "pT_s1_arne_tagged"," p_{T} topjets (arne tagged)",50,0,2000));
  Book( TH1F( "pT_s1_arne_tagged_ly"," p_{T} topjets (arne tagged)",50,0,2000));
   Book( TH1F( "pT_s1_tobias_tagged"," p_{T} topjets (tobias tagged)",50,0,2000));
  Book( TH1F( "pT_s1_tobias_tagged_ly"," p_{T} topjets (tobias tagged)",50,0,2000));
  Book( TH1F( "pT_s2_hep2_tagged_ly"," p_{T} CA15 (hep tagged)",50,0,2000));
  Book( TH1F( "pT_s2_hep2_tagged"," p_{T} CA15 (hep tagged)",50,0,2000));
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
  Book( TH1F( "MVA_eff", "efficiency of MVA",200,-1,1));
  // tmva_tagger=TMVA_tagger::Instance();
  tmva_tagger=new TMVA_tagger();
   Showerdeconstruction_tagger= new Showerdeconstruction();
  //tmva_tagger->Set_Reader("Uncorr+3");
  // tmva_tagger->Set_Reader("NPVweight");
  //tmva_tagger->Set_Reader("bestown70");
  //tmva_tagger->Set_Reader("pT");
  /*TMVA::Reader *readerx=new TMVA::Reader();
  float Tmass;
  float S1mass;
  float S2mass;
  float S3mass;
  float N;
  float Q;
  float t3t2;
  reader->AddVariable( "TopJet_mass", &Tmass );
   reader->AddVariable( "Subjets12_mass", &S1mass);
    reader->AddVariable( "Subjets13_mass", &S2mass);
    reader->AddVariable( "Subjets23_mass", &S3mass);
    reader->AddVariable( "TopJet_Nsubjets", &N);
    reader->AddVariable( "TopJet_Qjets_volatility", &Q);
    reader->AddVariable("t3/t2 := TopJet_tau3/TopJet_tau2", &t3t2);
    reader->BookMVA("BDTG", "/nfs/dust/cms/user/tlapsien/TMVA/weights/Qref_weight_BDTG.weights.xml");
    Double_t mvaValue = reader->EvaluateMVA("BDTG");*/
    

}


void Mistagcontrol::Fill()
{
   // important: get the event weight
  EventCalc* calc = EventCalc::Instance();
  double weight = calc -> GetWeight();
  // std::cout<<weight<<std::endl;
  TopFitCalc* topfit = TopFitCalc::Instance();
  //TMVA_tagger* tmva_tagger=TMVA_tagger::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();
  int NTopJets = bcc-> topjets -> size();
  Hist("NTopJets")->Fill(NTopJets, weight);
  Hist("NTopJets_ly")->Fill(NTopJets, weight);
  int nbtags=0;
  for(unsigned int i =0; i<bcc->jets->size();++i){
    Jet jet = bcc->jets->at(i);
    if(IsTagged(jet,e_CSVT)) nbtags++;
  }
  Hist("NbJets")->Fill(nbtags,weight);
  Hist("NbJets_ly")->Fill(nbtags,weight); 

  //TMVA
  Double_t mva_value;
  
  //TMVA::Reader *r=new TMVA::Reader();
  //tmva_tagger->Set_Reader("Qref_weight");
  TopJet antitaggedJet;
  bool antitag=false;
  for(unsigned int j=0; j<bcc->topjets->size();++j){
    TopJet topjet =bcc->topjets->at(j);
 
    double mmin3=0;
    double mjet3=0;
    int nsubjets3=0;
    TopTag(topjet,mjet3,nsubjets3,mmin3);
    if(mmin3<30 && mjet3>140 &&mjet3<250 && !antitag) {
      antitag=true;
      antitaggedJet=topjet;
    }
  }  

 
  
  if(antitag) for (unsigned int i =0; i<bcc->topjets->size(); ++i)
    {
      TopJet topjet =  bcc->topjets->at(i);
      /*  tmva_tagger->push_variables(topjet);
	  mva_value=tmva_tagger->GetMVA_value();*/
     
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
  
      
      Hist("pT") -> Fill(topjet.pt(), weight);
      Hist("pT_ly") -> Fill(topjet.pt(), weight);

      //selection to test efficiency
      bool jet_distance=true;
      bool selection_thad=true;
      if(topjet.v4()!=antitaggedJet.v4()){
      Hist("pT_s1")->Fill(topjet.pt(),weight);
      Hist("pT_s1_ly")->Fill(topjet.pt(),weight);
	     
       //std::cout<<mva_value<<std::endl;
      Hist("MVA")->Fill(mva_value,weight);
      Hist("MVA_ly")->Fill(mva_value,weight);
      if(topjet.pt()>350){
      Hist("MVA350")->Fill(mva_value,weight);
      Hist("MVA350_ly")->Fill(mva_value,weight);}
      Hist("MVA_pt")->Fill(mva_value,topjet.pt());
      double mmin2=0;
      double mjet2=0;
      int nsubjets2=0;
      if(TopTag(topjet,mjet2,nsubjets2,mmin2)) {
	Hist("pT_s1_cms_tagged")->Fill(topjet.pt(),weight);
	Hist("pT_s1_cms_tagged_ly")->Fill(topjet.pt(),weight);
      }
      /* if(HepTopTagFull(topjet,calc->GetPFParticles())) {
	Hist("pT_s1_hep_tagged")->Fill(topjet.pt(),weight);
	Hist("pT_s1_hep_tagged_ly")->Fill(topjet.pt(),weight);
	}*/
        for (unsigned int k =0; k<bcc->higgstagjets->size(); ++k){
	  TopJet CA15jet =  bcc->higgstagjets->at(k);
	  if((sqrt(pow(topjet.phi()-CA15jet.phi(),2)+pow(topjet.eta()-CA15jet.eta(),2))<1.5) && HepTopTag(CA15jet)){
	    Hist("pT_s1_hep_tagged")->Fill(topjet.pt(),weight);
	    Hist("pT_s1_hep_tagged_ly")->Fill(topjet.pt(),weight);
	 }
	
	}

      if(tmva_tagger->IsTagged("NPVweight",topjet,0.840,mva_value)){
	Hist("pT_s1_arne_tagged")->Fill(topjet.pt(),weight);
	Hist("pT_s1_arne_tagged_ly")->Fill(topjet.pt(),weight);  
      }
       if(tmva_tagger->IsTobiasTagged(topjet)){
	Hist("pT_s1_tobias_tagged")->Fill(topjet.pt(),weight);
	Hist("pT_s1_tobias_tagged_ly")->Fill(topjet.pt(),weight);  
      }
      }
      //shower deconstruction
      double chi_2= Showerdeconstruction_tagger->Chi(topjet);
      double chi= Showerdeconstruction_tagger->ChiMicro(topjet);
      Hist("SD")->Fill(log(chi),weight);
      Hist("SD2")->Fill(log(chi_2),weight);
      Hist("NMicrojets")->Fill(Showerdeconstruction_tagger->GetNmicrojets(topjet));

      Hist("eta") -> Fill(topjet.eta(), weight);
      Hist("eta_ly") -> Fill(topjet.eta(), weight);
      Hist("phi") -> Fill(topjet.phi(), weight);
      Hist("phi_ly") -> Fill(topjet.phi(), weight);
      
      double mmin=0;
      double mjet=0;
      int nsubjets=0;
      TopTag(topjet,mjet,nsubjets,mmin);
     
      Hist( "MJet" )->Fill( mjet, weight );
      Hist( "MJet_ly" )->Fill( mjet, weight );
      
      if(nsubjets>=3) 
	{
	  Hist( "Mmin" )->Fill( mmin, weight );
	  Hist( "Mmin_ly" )->Fill( mmin, weight );
	  
	}
      Hist( "NSubjets" )->Fill( nsubjets, weight ); 
      Hist( "NSubjets_ly" )->Fill( nsubjets, weight ); 
       Hist("MVA2")->Fill(mva_value,weight);
      Hist("MVA2_ly")->Fill(mva_value,weight);
     
    }

  //some controlplots for CA15
  if (antitag)  for (unsigned int i =0; i<bcc->higgstagjets->size(); ++i){
      TopJet CA15jet =  bcc->higgstagjets->at(i);
	Hist("pT_s2")->Fill(CA15jet.pt(),weight);
	Hist("pT_s2_ly")->Fill(CA15jet.pt(),weight);
	if(HepTopTag(CA15jet)) {
	  Hist("pT_s2_hep2_tagged")->Fill(CA15jet.pt(),weight);
	  Hist("pT_s2_hep2_tagged_ly")->Fill(CA15jet.pt(),weight);
	} 
    }


  
  
  


  sort(bcc->topjets->begin(), bcc->topjets->end(), HigherPt());
  for (unsigned int i =0; i<=3; ++i)
    {
      if (bcc->topjets->size()> i)
	{
	  
	  TopJet topjet =  bcc->topjets->at(i); 
	  TString hname = TString::Format("pT_%d", i+1);
	  Hist(hname)->Fill(topjet.pt(),weight);
	  TString hname_ly = TString::Format("pT_%d_ly", i+1);
	  Hist(hname_ly)->Fill(topjet.pt(),weight);
	  TString hname_eta = TString::Format("eta_%d", i+1);
	  Hist(hname_eta)->Fill(topjet.eta(),weight);
	  TString hname_eta_ly = TString::Format("eta_%d_ly", i+1);
	  Hist(hname_eta_ly)->Fill(topjet.eta(),weight);
	  TString hname_phi = TString::Format("phi_%d", i+1);
	  Hist(hname_phi)->Fill(topjet.phi(),weight);
	  TString hname_phi_ly = TString::Format("phi_%d_ly", i+1);
	  Hist(hname_phi_ly)->Fill(topjet.phi(),weight);
	  
	  double mmin=0;
	  double mjet=0;
	  int nsubjets=0;
	  TopTag(topjet,mjet,nsubjets,mmin);
	  TString hname_MJet = TString::Format("MJet_%d", i+1);
	  Hist(hname_MJet )->Fill( mjet, weight );
	  TString hname_MJet_ly = TString::Format("MJet_%d_ly", i+1);
	  Hist(hname_MJet_ly )->Fill( mjet, weight );
	  if(nsubjets>=3) 
	    {
	      TString hname_Mmin = TString::Format("Mmin_%d", i+1);
	      Hist(hname_Mmin  )->Fill( mmin, weight );
	      TString hname_Mmin_ly = TString::Format("Mmin_%d_ly", i+1);
	      Hist( hname_Mmin_ly )->Fill( mmin, weight );
	    }
	  TString hname_NSubjets = TString::Format("NSubjets_%d", i+1);
	  Hist( hname_NSubjets )->Fill( nsubjets, weight ); 
	  TString hname_NSubjets_ly = TString::Format("NSubjets_%d_ly", i+1);
	  Hist( hname_NSubjets_ly )->Fill( nsubjets, weight ); 


	  if (i==0){
	    
	    Hist("QjetsVol_1")->Fill(topjet.qjets_volatility(), weight);
	    Hist("QjetsVol_1_ly")->Fill(topjet.qjets_volatility(), weight);

	    Hist( "Nsubjettiness1_1")->Fill(topjet.tau1(), weight);
	    Hist( "Nsubjettiness1_1_ly")->Fill(topjet.tau1(), weight);

	    Hist( "Nsubjettiness2_1")->Fill(topjet.tau2(), weight);
	    Hist( "Nsubjettiness2_1_ly")->Fill(topjet.tau2(), weight);

	    Hist( "Nsubjettiness3_1")->Fill(topjet.tau3(), weight);
	    Hist( "Nsubjettiness3_1_ly")->Fill(topjet.tau3(), weight);

	    if (topjet.tau2()>0){
	      Hist( "Nsubjettiness3_2_1")->Fill(topjet.tau3()/topjet.tau2(), weight);
	      Hist( "Nsubjettiness3_2_1_ly")->Fill(topjet.tau3()/topjet.tau2(), weight);
	    }	  
	  }

	  
	}
      
    }
 
}
 
