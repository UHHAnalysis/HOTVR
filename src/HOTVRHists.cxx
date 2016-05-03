#include <time.h>
#include "include/HOTVRHists.h"
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
HOTVRHists::HOTVRHists(const char* name, double mass_min, double mass_max, double mmin_min, double mmin_max, int Nsubjets_min, int Nsubjets_max, double ptfraction_min, double ptfraction_max, double subjet1pt, double subjet2pt) : BaseHists(name)
{
  //  std::cout<<name<<std::endl;
 m_mass_min=mass_min;
  m_mass_max=mass_max;
  m_mmin_min=mmin_min;
  m_mmin_max=mmin_max;
  m_Nsubjets_min=Nsubjets_min;
  m_Nsubjets_max=Nsubjets_max;
  m_ptfraction_min=ptfraction_min;
  m_ptfraction_max=ptfraction_max;
  m_subjet1pt=subjet1pt;
  m_subjet2pt=subjet2pt;
  // named default constructor
   
}

HOTVRHists::~HOTVRHists()
{
  // default destructor, does nothing
}

void HOTVRHists::Init()
{
  // book all histograms here

 
  Book(TH1F("topgen","generated had tops",15,-0.5,14.5));
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
  
  
  
  
  Book( TH1D( "Nsub", ";#tau_{3}/#tau_{2};Events", 100, 0., 1.));
  Book( TH1D( "Nsub_tagged", ";#tau_{3}/#tau_{2};Events", 100, 0., 1.));
  
  
  Book( TH1D( "weight", ";weight;Events", 400, 0., 1000.));
  
  
  
  //Variable R
  Book( TH1D( "Jetmass2", ";Top jet mass (GeV);Events", 100, 0, 500));
  Book( TH1D( "Subjetmass2", ";subjet mass (GeV);Events", 50, 0, 500));
  Book( TH1D( "mmin2", ";Top jet mass (GeV);Events", 50, 0, 200));
  Book( TH1D( "mmin3", ";Top jet mass (GeV);Events", 50, 0, 200));
  Book(TH1D("rho","rho",60,-10,200));
  Book(TH1D("CMS_pt","cms",50,0,2000));
  Book(TH1D("TL_pt","tl",50,0,2000));
  Book(TH2F("pT_njets","njets vs pt",1000,0,2000,10,0,10));
  Book(TH2F("VariableR","VariableR",1000,0,2000,20,0,40));
  Book(TH2F("VariableR2","VariableR",1000,0,2000,20,0,40));
  Book(TH2F("VariableR_area","VariableR area",1000,0,2000,20,0,2));
  Book(TH2F("Area2d","Area vs pt",1000,0,2000,100,0,20));
  Book(TH1D("mmin","mmin",50,0,100));
  Book(TH1D("massjump_mmin200","massjump mmin",50,0,100));
  Book(TH1D("massjump_mmin400","massjump mmin",50,0,100));
  Book(TH1D("massjump_mmin600","massjump mmin",50,0,100));
  Book(TH1D("massjump_mmin800","massjump mmin",50,0,100));
  Book(TH1D("massjump_mmin1000","massjump mmin",50,0,100));
  Book(TH1D("mjet","mjet",100,0,300));
  Book(TH1D("mjet400","mjet",100,0,300));
  Book(TH1D("mjet200","mjet",100,0,300));
  Book(TH1D("mjet600","mjet",100,0,300));
  Book(TH1D("mjet800","mjet",100,0,300));
  Book(TH1D("mjet1000","mjet",100,0,300));

  

  Book(TH1D("mjet_before","mjet",100,0,300));
  Book(TH1D("mjet500_before","mjet",100,0,300));
  Book(TH1D("mjet200_before","mjet",100,0,300));
  Book(TH1D("mjet600_before","mjet",100,0,300));
  Book(TH1D("mjet800_before","mjet",100,0,300));
  Book(TH1D("mjet1000_before","mjet",100,0,300));
   
 
  Book(TH1D("nsubjets","nsubjets",5,-0.5,4.5));
  Book(TH1D("nsubjets_massjump","nsubjets massjump",20,-0.5,19.5));
  Book(TH1D("nsubjets_massjump200","nsubjets massjump",20,-0.5,19.5));
  Book(TH1D("nsubjets_massjump400","nsubjets massjump",20,-0.5,19.5));
  Book(TH1D("nsubjets_massjump600","nsubjets massjump",20,-0.5,19.5));
  Book(TH1D("nsubjets_massjump800","nsubjets massjump",20,-0.5,19.5));
  Book(TH1D("nsubjets_massjump1000","nsubjets massjump",20,-0.5,19.5));
  Book(TH1D("NJets","njets",20,-0.5,19.5));
  Book(TH1D("NJets-1","njets radiation",20,-0.5,19.5));
  Book(TH1D("NJets-2","njets beam jets",20,-0.5,19.5));
  Book(TH1D("NJetstop","NJets top",20,-0.5,19.5));
  Book(TH1D("NJets200","njets",20,-0.5,19.5));
  Book(TH1D("NJetstop200","NJets top",20,-0.5,19.5));
  Book(TH1D("NJets400","njets",20,-0.5,19.5));
  Book(TH1D("NJetstop400","NJets top",20,-0.5,19.5));
  Book(TH1D("NJets600","njets",20,-0.5,19.5));
  Book(TH1D("NJetstop600","NJets top",20,-0.5,19.5));
  Book(TH1F("pT_beam"," p_{T} beam jets",50,0,2000));
  Book(TH1F("pT_radiation"," p_{T} radiation jets",50,0,2000));
  Book(TH1F("pT_fatjets"," p_{T} fatjets",50,0,2000));
 
   
 
  Book(TH1D("pT_jets","pT varRjet",100,0,2000));
  Book(TH1D("pT_topjets","pT varRjet",100,0,2000));
  //  Book(TH1D("Njets_var","pt vs njet",200,0,2000));
  Book(TH1D("Njets_norm","pt vs njet",200,0,2000));
  Book(TH1D("tau1","tau1",100,0,1.0));
  Book(TH1D("tau1_200","tau1 pt>200GeV",100,0,1.0));
  Book(TH1D("tau1_500","tau1 pt>500GeV",100,0,1.0));
  Book(TH1D("tau1_600","tau1 pt>600GeV",100,0,1.0));
  Book(TH1D("tau1_800","tau1 pt>800GeV",100,0,1.0));
  Book(TH1D("tau1_1000","tau1 pt>1000GeV",100,0,1.0));
  Book(TH1D("tau2","tau2",100,0,1.0));
  Book(TH1D("tau2_200","tau2 pt>200GeV",100,0,1.0));
  Book(TH1D("tau2_500","tau2 pt>500GeV",100,0,1.0));
  Book(TH1D("tau2_600","tau2 pt>600GeV",100,0,1.0));
  Book(TH1D("tau2_800","tau2 pt>800GeV",100,0,1.0));
  Book(TH1D("tau2_1000","tau1 pt>1000GeV",100,0,1.0));
  Book(TH1D("tau3","tau3",100,0,1.0));
  Book(TH1D("tau3_200","tau3 pt>200GeV",100,0,1.0));
  Book(TH1D("tau3_500","tau3 pt>500GeV",100,0,1.0));
  Book(TH1D("tau3_600","tau3 pt>600GeV",100,0,1.0));
  Book(TH1D("tau3_800","tau3 pt>800GeV",100,0,1.0));
  Book(TH1D("tau3_1000","tau3 pt>1000GeV",100,0,1.0));
  Book(TH1D("tau4","tau4",100,0,1.0));
  Book(TH1D("tau4_200","tau4 pt>200GeV",100,0,1.0));
  Book(TH1D("tau4_500","tau4 pt>500GeV",100,0,1.0));
  Book(TH1D("tau4_600","tau4 pt>600GeV",100,0,1.0));
  Book(TH1D("tau4_800","tau4 pt>800GeV",100,0,1.0));
  Book(TH1D("tau4_1000","tau2tau1 pt>1000GeV",100,0,1.0));
  Book(TH1D("tau2tau1","tau2tau1",100,0,1.0));
  Book(TH1D("tau2tau1_200","tau2tau1 pt>200GeV",100,0,1.0));
  Book(TH1D("tau2tau1_500","tau2tau1 pt>500GeV",100,0,1.0));
  Book(TH1D("tau2tau1_600","tau2tau1 pt>600GeV",100,0,1.0));
  Book(TH1D("tau2tau1_800","tau2tau1 pt>800GeV",100,0,1.0));
  Book(TH1D("tau2tau1_1000","tau2tau1 pt>1000GeV",100,0,1.0));
  Book(TH1D("tau3tau2","tau2tau1",100,0,1.0));
  Book(TH1D("tau3tau1","tau2tau1",100,0,1.0));
  Book(TH1D("tau3tau2_200","tau3tau2 pt>200GeV",100,0,1.0));
  Book(TH1D("tau3tau2_500","tau3tau2 pt>500GeV",100,0,1.0));
  Book(TH1D("tau3tau2_600","tau3tau2 pt>600GeV",100,0,1.0));
  Book(TH1D("tau3tau2_800","tau3tau2 pt>800GeV",100,0,1.0));
  Book(TH1D("tau3tau2_1000","tau3tau2 pt>1000GeV",100,0,1.0));
  
  //subjetmass
  Book(TH1D("subjet01_massjump","subjet01 mass",100,0,200));
  Book(TH1D("subjet02_massjump","subjet02 mass",100,0,200));
  Book(TH1D("subjet12_massjump","subjet12 mass",100,0,200));
  Book(TH1D("subjet03_massjump","subjet03 mass",100,0,200));
  Book(TH1D("subjet13_massjump","subjet13 mass",100,0,200));
  Book(TH1D("subjet23_massjump","subjet23 mass",100,0,200));
  Book(TH1D("Wmass","reconstructed W-mass",100,0,200));
  Book(TH1D("pT_pf","pT pf",100,0,1000));
  Book(TH1D("sd_mass","softdrop mass",100,0,300));
  
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
 
  Book(TH1D("ptfraction1_200","p_{T,sub1}/p_{T,jet} [GeV]",20,0,1));
  Book(TH1D("ptfraction1_400","p_{T,sub1}/p_{T,jet} [GeV]",20,0,1));
  Book(TH1D("ptfraction1_600","p_{T,sub1}/p_{T,jet} [GeV]",20,0,1));
  Book(TH1D("ptfraction1_800","p_{T,sub1}/p_{T,jet} [GeV]",20,0,1));
  Book(TH1D("ptfraction1_1000","p_{T,sub1}/p_{T,jet} [GeV]",20,0,1));
  
}

void HOTVRHists::SetIdVersion(TString s)
{
  idVersion=s;
}

void HOTVRHists::Fill()
{

  double rho(600.);
  double mu(30.), theta(0.7), mw(0.),mtopLow(140.),mtopHigh(220.), pt_cut(20.);
  double min_r(0.1), max_r(1.5);
  double ptmin(150.);
  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();
  bool IsRealData = calc->IsRealData();
 
  // important: get the event weight
  double weight = calc->GetWeight();
  int bestjetindex2=-1;
 
  if(!IsRealData && (idVersion=="TTbarEff" || idVersion=="TTbarEff_1" || idVersion=="TTbarEff_gen" || idVersion=="TTbarEff_gen2" || idVersion.Contains("Zprime")))
  {
    


    //fill gen particles for clustering
    std::vector<GenParticle>* genparticles = calc->GetGenParticles();
    std::vector<fastjet::PseudoJet> genvector;
    for(unsigned int tx=0;tx<genparticles->size();tx++){
      
      TLorentzVector particle;
      particle.SetPtEtaPhiE(genparticles->at(tx).pt(),genparticles->at(tx).eta(),genparticles->at(tx).phi(),genparticles->at(tx).energy());
      
      fastjet::PseudoJet gen_particle( genparticles->at(tx).pt()*cos(genparticles->at(tx).phi()),
				       genparticles->at(tx).pt()*sin(genparticles->at(tx).phi()),
				       genparticles->at(tx).pt()*sinh(genparticles->at(tx).eta()),
				       genparticles->at(tx).energy() );
     

      if(genparticles->at(tx).status()==1) genvector.push_back(gen_particle);
    }
  
    TTbarGen2* Decay = calc->GetTTbarGen2();//reconstruct TTbar event with Pythia8
    
    
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
   
    
    
    
    
    for (unsigned int j=0; j<Had_Tops.size(); j++)//if (isHadronic)
      {
	


	Hist("pT_s1_all")->Fill(Had_Tops[j].pt(),weight); //Fill denominator hist
	Hist("eta_s1_all")->Fill(Had_Tops[j].eta(),weight);//Fill denominator hist
	
	for(int roc=0;roc<100;roc++) {((TH2D*)Hist("nsub_eff_norm"))->Fill(Had_Tops[j].pt(),roc/100.,weight);}//Fill denominator hist for ROC with nsubjettiness
	EventCalc* calc = EventCalc::Instance();
	fastjet::JetDefinition *JetDef ;
	
	

	//cluster with HOTVR
	HOTVR plugin_hotvr(mu, theta,min_r, max_r,rho,pt_cut, HOTVR::CALIKE);//call HOTVR algorithm
 	fastjet::JetDefinition jet_def(&plugin_hotvr);
	fastjet::ClusterSequence clust_seq(genvector, jet_def);
	

	std::vector<fastjet::PseudoJet> hotvr_jets,rejected_jets,soft_jets ; //vector of hotvr_jets, jets that were rejcted durning the clustering procedure and soft jets
	
	//get vector from the plugin
	hotvr_jets=sorted_by_pt(plugin_hotvr.get_jets());
	rejected_jets=plugin_hotvr.get_rejected_cluster();
	soft_jets=plugin_hotvr.get_soft_cluster();

	
	//Fill mass before matching
	for(unsigned int i=0;i<hotvr_jets.size();i++){
	  double mjet=hotvr_jets[i].m();
	  if(Had_Tops[j].pt()>200)  Hist("mjet200_before")->Fill(mjet,weight);
	  if(Had_Tops[j].pt()>500)  Hist("mjet500_before")->Fill(mjet,weight);
	  if(Had_Tops[j].pt()>600)  Hist("mjet600_before")->Fill(mjet,weight);
	  if(Had_Tops[j].pt()>800)  Hist("mjet800_before")->Fill(mjet,weight);
	  if(Had_Tops[j].pt()>1000)  Hist("mjet1000_before")->Fill(mjet,weight);
	  
	   
	}
	
  	Hist("NJets")->Fill(hotvr_jets.size(),weight); //Fill Number of jets before matching
	
	//matching: Find jet with smallest distance to hadronic top quark
	double minDeltaR2=10000;
	int bestjetindex2=-1;
	for(unsigned int i=0; i<hotvr_jets.size(); ++i)
	{
	  double DeltaR2=sqrt(pow(Had_Tops[j].phi()-hotvr_jets[i].phi_std(),2)+pow(Had_Tops[j].eta()-hotvr_jets[i].pseudorapidity(),2));
	  if(DeltaR2<minDeltaR2 )
	    {
	      minDeltaR2=DeltaR2;
	      bestjetindex2=i;
	    }
	}


	//matching: is smallest distance smaller than radius of the jet?
	if(bestjetindex2>-1) if (minDeltaR2>hotvr_jets[bestjetindex2].user_info<HOTVRinfo>().radius()) bestjetindex2=-1;//protect against large deltaR

	
	//if radius is 0: do not match
	if(bestjetindex2>-1 ) if(hotvr_jets[bestjetindex2].user_info<HOTVRinfo>().radius()==0)  bestjetindex2=-1;

	
	//matching succesfull?
	if(bestjetindex2>-1 )
	  {
	   
	    std::vector<fastjet::PseudoJet> SortedSubJets=sorted_by_pt(hotvr_jets[bestjetindex2].user_info<HOTVRinfo>().subjets());//Get subjets of hotvr jets
	    double jet_radius=hotvr_jets[bestjetindex2].user_info<HOTVRinfo>().radius(); //radius of probejet
	    fastjet::PseudoJet jet=hotvr_jets.at(bestjetindex2); //probejet


	    //calculate the substructure variable mmin
	    double m12=0;
	    double m01=0;
	    double m02=0;
	    double mmin=0;
	    if(SortedSubJets.size()>2){
	      m01 = 0;
	      m01=(SortedSubJets[0]+SortedSubJets[1]).m();
	      m02= 0;
	      m02=(SortedSubJets[0]+SortedSubJets[2]).m();
	      m12 = 0;
	      m12 = (SortedSubJets[1]+SortedSubJets[2]).m();
	      Hist("subjet01_massjump")->Fill(m01,weight);
	      Hist("subjet02_massjump")->Fill(m02,weight);
	      Hist("subjet12_massjump")->Fill(m12,weight);
	    }
	    mmin = std::min(m01,std::min(m02,m12));
	    double jet_mass=jet.m();
	    double jet_mmin=mmin;

	    //Fill subjet pT
	    double jet_subjet1pt=0;
	    double jet_subjet2pt=0;
	    if(SortedSubJets.size()>0) jet_subjet1pt=SortedSubJets.at(0).pt();
	    if(SortedSubJets.size()>1) jet_subjet2pt=SortedSubJets.at(1).pt();
	    
	    if(/*jet_mass>m_mass_min && jet_mass<m_mass_max && jet_mmin>m_mmin_min && jet_mmin<m_mmin_max && jet_nsubjets>m_Nsubjets_min && jet_nsubjets<m_Nsubjets_max && jet_ptfraction>m_ptfraction_min && jet_ptfraction<m_ptfraction_max && jet_subjet1pt>m_subjet1pt && jet_subjet2pt>m_subjet2pt*/true){ 
	      
	      //Fill the number of subjets
	      Hist("nsubjets_massjump")->Fill(SortedSubJets.size(),weight);
	      if(Had_Tops[j].pt()>200 &&Had_Tops[j].pt()<400 ) Hist("nsubjets_massjump200")->Fill(SortedSubJets.size(),weight);
	      if(Had_Tops[j].pt()>400 &&Had_Tops[j].pt()<600) Hist("nsubjets_massjump400")->Fill(SortedSubJets.size(),weight);
	      if(Had_Tops[j].pt()>600&&Had_Tops[j].pt()<800) Hist("nsubjets_massjump600")->Fill(SortedSubJets.size(),weight);
	      if(Had_Tops[j].pt()>800&&Had_Tops[j].pt()<1000) Hist("nsubjets_massjump800")->Fill(SortedSubJets.size(),weight);
	      if(Had_Tops[j].pt()>1000&&Had_Tops[j].pt()<1200) Hist("nsubjets_massjump1000")->Fill(SortedSubJets.size(),weight);
	      
	      //Fill mmin
	      if(Had_Tops[j].pt()>200 &&  Had_Tops[j].pt()<400) Hist("massjump_mmin200")->Fill(mmin,weight);
	      if(Had_Tops[j].pt()>400 &&  Had_Tops[j].pt()<600) Hist("massjump_mmin400")->Fill(mmin,weight);
	      if(Had_Tops[j].pt()>600 &&  Had_Tops[j].pt()<800) Hist("massjump_mmin600")->Fill(mmin,weight);
	      if(Had_Tops[j].pt()>800 &&  Had_Tops[j].pt()<1000) Hist("massjump_mmin800")->Fill(mmin,weight);
	      if(Had_Tops[j].pt()>1000 &&  Had_Tops[j].pt()<1200) Hist("massjump_mmin1000")->Fill(mmin,weight);

	      //Fill pT fraction
	      if(Had_Tops[j].pt()>200 &&  Had_Tops[j].pt()<400) Hist("ptfraction1_200")->Fill( SortedSubJets.at(0).pt()/jet.pt(),weight);
	      if(Had_Tops[j].pt()>400 &&  Had_Tops[j].pt()<600) Hist("ptfraction1_400")->Fill( SortedSubJets.at(0).pt()/jet.pt(),weight);
	      if(Had_Tops[j].pt()>600 &&  Had_Tops[j].pt()<800) Hist("ptfraction1_600")->Fill( SortedSubJets.at(0).pt()/jet.pt(),weight);
	      if(Had_Tops[j].pt()>800 &&  Had_Tops[j].pt()<1000) Hist("ptfraction1_800")->Fill( SortedSubJets.at(0).pt()/jet.pt(),weight);
	      if(Had_Tops[j].pt()>1000 &&  Had_Tops[j].pt()<1200) Hist("ptfraction1_1000")->Fill( SortedSubJets.at(0).pt()/jet.pt(),weight);
	  
	      //Fill jet mass
	      double mjet;
	      mjet=jet.m();
	      if(Had_Tops[j].pt()>200 &&  Had_Tops[j].pt()<400)  Hist("mjet200")->Fill(mjet,weight);
	      if(Had_Tops[j].pt()>400 &&  Had_Tops[j].pt()<600)  Hist("mjet400")->Fill(mjet,weight);
	      if(Had_Tops[j].pt()>600 &&  Had_Tops[j].pt()<800)  Hist("mjet600")->Fill(mjet,weight);
	      if(Had_Tops[j].pt()>800 &&  Had_Tops[j].pt()<1000)  Hist("mjet800")->Fill(mjet,weight);
	      if(Had_Tops[j].pt()>1000 &&  Had_Tops[j].pt()<1200)  Hist("mjet1000")->Fill(mjet,weight);

	      //Fill radius as a function of pT
	      ((TH2D*)Hist("VariableR"))->Fill(jet.perp(),jet_radius*10,weight); 
	  
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
	  	
	      //Fill nominator hists
	      if(SortedSubJets.size()>2)  Hist("pT_s1_hotvr_tagged")->Fill(Had_Tops[j].pt(),weight);
	      if(jet.m()>mtopLow &&  jet.m()<mtopHigh)  Hist("pT_s1_hotvr2_tagged")->Fill(Had_Tops[j].pt(),weight);
	      if(SortedSubJets.size()>2 && jet.m()>mtopLow &&  jet.m()<mtopHigh)  Hist("pT_s1_hotvr3_tagged")->Fill(Had_Tops[j].pt(),weight);
	      if(SortedSubJets.size()>2 && jet.m()>mtopLow &&  jet.m()<mtopHigh  &&SortedSubJets.at(0).pt()>30 && SortedSubJets.at(1).pt()>30 &&  SortedSubJets.at(0).pt()/jet.pt()<0.8 &&mmin>50) {
		Hist("pT_s1_hotvr6_tagged")->Fill(Had_Tops[j].pt(),weight);
		Hist("eta_s1_hotvr6_tagged")->Fill(Had_Tops[j].eta(),weight);
	      }
	    
	     
	     
	      
	      
	      
	      //calculate nsubjettiness
	      if(jet.constituents().size()>0){
		double tau1,tau2,tau3,tau4;
		JetPropsPseudo jp(&jet);
		tau1 = jp.GetNsubjettiness(1, Njettiness::onepass_kt_axes, 1.,jet_radius);
		tau2 = jp.GetNsubjettiness(2, Njettiness::onepass_kt_axes, 1.,jet_radius);
		tau3 = jp.GetNsubjettiness(3, Njettiness::onepass_kt_axes, 1.,jet_radius);
		tau4 = jp.GetNsubjettiness(4, Njettiness::onepass_kt_axes, 1.,jet_radius);
	      
	      //Fill nsubjettiness
	      Hist("tau1")->Fill(tau1,weight);
	      if(Had_Tops[j].pt()>200) Hist("tau1_200")->Fill(tau1,weight);
	      if(Had_Tops[j].pt()>500) Hist("tau1_500")->Fill(tau1,weight);
	      if(Had_Tops[j].pt()>600) Hist("tau1_600")->Fill(tau1,weight);
	      if(Had_Tops[j].pt()>800) Hist("tau1_800")->Fill(tau1,weight);
	      if(Had_Tops[j].pt()>1000) Hist("tau1_1000")->Fill(tau1,weight);
	      
	      Hist("tau2")->Fill(tau2,weight);
	      if(Had_Tops[j].pt()>200) Hist("tau2_200")->Fill(tau2,weight);
	      if(Had_Tops[j].pt()>500) Hist("tau2_500")->Fill(tau2,weight);
	      if(Had_Tops[j].pt()>600) Hist("tau2_600")->Fill(tau2,weight);
	      if(Had_Tops[j].pt()>800) Hist("tau2_800")->Fill(tau2,weight);
	      if(Had_Tops[j].pt()>1000) Hist("tau2_1000")->Fill(tau2,weight);
	      Hist("tau3")->Fill(tau3,weight);
	      if(Had_Tops[j].pt()>200) Hist("tau3_200")->Fill(tau3,weight);
	      if(Had_Tops[j].pt()>500) Hist("tau3_500")->Fill(tau3,weight);
	      if(Had_Tops[j].pt()>600) Hist("tau3_600")->Fill(tau3,weight);
	      if(Had_Tops[j].pt()>800) Hist("tau3_800")->Fill(tau3,weight);
	      if(Had_Tops[j].pt()>1000) Hist("tau3_1000")->Fill(tau3,weight);
	      Hist("tau4")->Fill(tau4,weight);
	      if(Had_Tops[j].pt()>200) Hist("tau4_200")->Fill(tau4,weight);
	      if(Had_Tops[j].pt()>500) Hist("tau4_500")->Fill(tau4,weight);
	      if(Had_Tops[j].pt()>600) Hist("tau4_600")->Fill(tau4,weight);
	      if(Had_Tops[j].pt()>800) Hist("tau4_800")->Fill(tau4,weight);
	      if(Had_Tops[j].pt()>1000) Hist("tau4_1000")->Fill(tau4,weight);
	      Hist("tau2tau1")->Fill(tau2/tau1,weight);
	      if(Had_Tops[j].pt()>200)  Hist("tau2tau1_200")->Fill(tau2/tau1,weight);
	      if(Had_Tops[j].pt()>500)  Hist("tau2tau1_500")->Fill(tau2/tau1,weight);
	      if(Had_Tops[j].pt()>600)  Hist("tau2tau1_600")->Fill(tau2/tau1,weight);
	      if(Had_Tops[j].pt()>800)  Hist("tau2tau1_800")->Fill(tau2/tau1,weight);
	      if(Had_Tops[j].pt()>1000)  Hist("tau2tau1_1000")->Fill(tau2/tau1,weight);
	      Hist("tau3tau2")->Fill(tau3/tau2,weight);
	      if(Had_Tops[j].pt()>200)  Hist("tau3tau2_200")->Fill(tau3/tau2,weight);
	      if(Had_Tops[j].pt()>500)  Hist("tau3tau2_500")->Fill(tau3/tau2,weight);
	      if(Had_Tops[j].pt()>600)  Hist("tau3tau2_600")->Fill(tau3/tau2,weight);
	      if(Had_Tops[j].pt()>800)  Hist("tau3tau2_800")->Fill(tau3/tau2,weight);
	      if(Had_Tops[j].pt()>1000)  Hist("tau3tau2_1000")->Fill(tau3/tau2,weight);
	      
	      //Fill other nomiantors
	      if(SortedSubJets.size()>2 && jet.m()>mtopLow &&  jet.m()<mtopHigh  &&SortedSubJets.at(0).pt()>30 && SortedSubJets.at(1).pt()>30 &&   SortedSubJets.at(0).pt()/jet.pt()<0.8  &&tau3/tau2<0.7 &&mmin>50  )  Hist("pT_s1_hotvr5_tagged")->Fill(Had_Tops[j].pt(),weight);
	      if(SortedSubJets.size()>2 && jet.m()>mtopLow &&  jet.m()<mtopHigh  &&SortedSubJets.at(0).pt()>30 && SortedSubJets.at(1).pt()>30 &&   SortedSubJets.at(0).pt()/jet.pt()<0.8 &&mmin>50  ) for(int roc=0;roc<100;roc++){
		  if(tau3/tau2<roc/100.){
		    ((TH2D*)Hist("nsub_eff"))->Fill(Had_Tops[j].pt(),roc/100.,weight);//Fill the nsubjettiness nominator
		    
		  }
		}}
	      
	    }
	  }
      }
    
    
   
  }
  

  if(!IsRealData && (idVersion.Contains("QCD15to3000") ||idVersion.Contains("QCD")) )
  {
    
    std::vector<GenParticle>* genparticles = calc->GetGenParticles();
    std::vector<fastjet::PseudoJet> genvector;
    //Fill gen particles for clustering
    for(unsigned int tx=0;tx<genparticles->size();tx++){
      //genparticles->at(tx).Print(genparticles);
      TLorentzVector particle;
      particle.SetPtEtaPhiE(genparticles->at(tx).pt(),genparticles->at(tx).eta(),genparticles->at(tx).phi(),genparticles->at(tx).energy());
      fastjet::PseudoJet gen_particle(particle.Px(),particle.Py(),particle.Pz(),particle.E());
      if(genparticles->at(tx).status()==1) genvector.push_back(gen_particle);
    }
    
    EventCalc* calc = EventCalc::Instance();
    fastjet::JetDefinition *JetDef ;
    std::vector<fastjet::PseudoJet> SortedSubJets;
    
    
    HOTVR plugin_hotvr(mu, theta,min_r, max_r,rho,pt_cut, HOTVR::CALIKE); //call hotvr with parameters
    fastjet::JetDefinition jet_def_hotvr(&plugin_hotvr);
       
   
    fastjet::ClusterSequence clust_seq_hotvr(genvector, jet_def_hotvr);
   

    std::vector<fastjet::PseudoJet> hotvr_jets,rejected_jets,soft_jets ; //vector of hotvr_jets, jets that were rejcted durning the clustering procedure and soft jets
	
    //get vector from the plugin
    hotvr_jets=sorted_by_pt(plugin_hotvr.get_jets());
    rejected_jets=plugin_hotvr.get_rejected_cluster();
    soft_jets=plugin_hotvr.get_soft_cluster();

   
    double minDeltaR2=10000;
    int bestjetindex2=-1;
  
  
    for(unsigned int j=0; j<genparticles->size();j++){
      if (abs(genparticles->at(j).pdgId())<=6 || abs(genparticles->at(j).pdgId())==21){
	if(genparticles->at(j).pt()>100) Hist("pT_s1_all")->Fill(genparticles->at(j).pt(),weight); //fill denominator hists for efficiecies
	for(int roc=0;roc<100;roc++)if(genparticles->at(j).pt()>100) {((TH2D*)Hist("nsub_eff_norm"))->Fill(genparticles->at(j).pt(),roc/100.,weight);}//Fill denominator hist for nsubjettiness scan
      }
    }


    //matching ******
    for(unsigned int i=0;i<hotvr_jets.size();i++){
  
   
      double highestPt=-1;
      int  bestjetindex=-1;
     
      for(unsigned int j=0; j<genparticles->size();j++)
	{
	  if (abs(genparticles->at(j).pdgId())<=6 || abs(genparticles->at(j).pdgId())==21)
	    {
	      double DeltaR2=sqrt(pow(genparticles->at(j).phi()-hotvr_jets[i].phi_std(),2)+pow(genparticles->at(j).eta()-hotvr_jets[i].pseudorapidity(),2));
	      if (DeltaR2<hotvr_jets[i].user_info<HOTVRinfo>().radius() && genparticles->at(j).pt()>100)
		{
		  bestjetindex=j;
		  highestPt=genparticles->at(j).pt();
		  
		}
	    }
	}
   
      //******
    
    
 
      if(bestjetindex>-1 ) if(hotvr_jets[i].user_info<HOTVRinfo>().radius()==0)  bestjetindex=-1;//is jet radius=0? Pass not matching
  
      
      //jet matched to generator parton? proceed
      if (bestjetindex>-1)
	{

	  std::vector<fastjet::PseudoJet> SortedSubJets=sorted_by_pt(hotvr_jets[i].user_info<HOTVRinfo>().subjets());//Get subjets of hotvr jets
	  double jet_radius=hotvr_jets[i].user_info<HOTVRinfo>().radius(); //radius of probejet
	  fastjet::PseudoJet jet=hotvr_jets.at(i); //probejet
	  
	  
	  
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
	    Hist("subjet01_massjump")->Fill(m01,weight);
	    Hist("subjet02_massjump")->Fill(m02,weight);
	    Hist("subjet12_massjump")->Fill(m12,weight);
	  }
	  mmin = std::min(m01,std::min(m02,m12));
	  double jet_mass=hotvr_jets[i].m();
	  double jet_mmin=mmin;
	  int jet_nsubjets=SortedSubJets.size();
	  double jet_ptfraction=SortedSubJets.at(0).pt()/hotvr_jets[i].pt();
	  double jet_subjet1pt=0;
	  double jet_subjet2pt=0;
	  if(SortedSubJets.size()>0) jet_subjet1pt=SortedSubJets.at(0).pt();//Fill subjet pT
	  if(SortedSubJets.size()>1) jet_subjet2pt=SortedSubJets.at(1).pt();//Fill subjet pT
      
       if(/*jet_mass>m_mass_min && jet_mass<m_mass_max && jet_mmin>m_mmin_min && jet_mmin<m_mmin_max && jet_nsubjets>m_Nsubjets_min && jet_nsubjets<m_Nsubjets_max && jet_ptfraction>m_ptfraction_min && jet_ptfraction<m_ptfraction_max && jet_subjet1pt>m_subjet1pt && jet_subjet2pt>m_subjet2pt*/true){


      
	 //Fill number of subjets
	 if(genparticles->at(bestjetindex).pt()>200 && genparticles->at(bestjetindex).pt()<400)  Hist("nsubjets_massjump200")->Fill(SortedSubJets.size(),weight);
	 if(genparticles->at(bestjetindex).pt()>500 && genparticles->at(bestjetindex).pt()<600)  Hist("nsubjets_massjump400")->Fill(SortedSubJets.size(),weight);
	 if(genparticles->at(bestjetindex).pt()>600 && genparticles->at(bestjetindex).pt()<800)  Hist("nsubjets_massjump600")->Fill(SortedSubJets.size(),weight);
	 if(genparticles->at(bestjetindex).pt()>800 && genparticles->at(bestjetindex).pt()<1000)  Hist("nsubjets_massjump800")->Fill(SortedSubJets.size(),weight);
	 if(genparticles->at(bestjetindex).pt()>1000 && genparticles->at(bestjetindex).pt()<1200)  Hist("nsubjets_massjump1000")->Fill(SortedSubJets.size(),weight);
     
	 //Fill mmin
	 if(genparticles->at(bestjetindex).pt()>200 && genparticles->at(bestjetindex).pt()<400)  Hist("massjump_mmin200")->Fill(mmin,weight);
	 if(genparticles->at(bestjetindex).pt()>400 && genparticles->at(bestjetindex).pt()<600)  Hist("massjump_mmin400")->Fill(mmin,weight);
	 if(genparticles->at(bestjetindex).pt()>600 && genparticles->at(bestjetindex).pt()<800)  Hist("massjump_mmin600")->Fill(mmin,weight);
	 if(genparticles->at(bestjetindex).pt()>800 && genparticles->at(bestjetindex).pt()<1000)  Hist("massjump_mmin800")->Fill(mmin,weight);
	 if(genparticles->at(bestjetindex).pt()>1000 && genparticles->at(bestjetindex).pt()<1200)  Hist("massjump_mmin1000")->Fill(mmin,weight);
	 
	 //Fill subjet pT
	 for (unsigned int bk =0; bk<=3; ++bk)
	   {
	     if (SortedSubJets.size()> bk)
	       {
	  	 TString hname = TString::Format("pT_subjet%d", bk+1);
		 Hist(hname)->Fill(SortedSubJets.at(bk).pt(),weight);
		 hname = TString::Format("pT_subjetfrac%d", bk+1);
		 Hist(hname)->Fill(SortedSubJets.at(bk).pt()/hotvr_jets[i].pt(),weight);
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
	 mjet=hotvr_jets[i].m();
	 if(genparticles->at(bestjetindex).pt()>200 && genparticles->at(bestjetindex).pt()<400)  Hist("mjet200")->Fill(mjet,weight);
	 if(genparticles->at(bestjetindex).pt()>400 && genparticles->at(bestjetindex).pt()<600)  Hist("mjet400")->Fill(mjet,weight);
	 if(genparticles->at(bestjetindex).pt()>600 && genparticles->at(bestjetindex).pt()<800)  Hist("mjet600")->Fill(mjet,weight);
	 if(genparticles->at(bestjetindex).pt()>800 && genparticles->at(bestjetindex).pt()<1000)  Hist("mjet800")->Fill(mjet,weight);
	 if(genparticles->at(bestjetindex).pt()>1000 && genparticles->at(bestjetindex).pt()<1200)  Hist("mjet1000")->Fill(mjet,weight);


	 //Fill pT fraction
	 if(genparticles->at(bestjetindex).pt()>200 && genparticles->at(bestjetindex).pt()<400)	Hist("ptfraction1_200")->Fill( SortedSubJets.at(0).pt()/hotvr_jets[i].pt(),weight);
	 if(genparticles->at(bestjetindex).pt()>400 && genparticles->at(bestjetindex).pt()<600)	Hist("ptfraction1_400")->Fill( SortedSubJets.at(0).pt()/hotvr_jets[i].pt(),weight);
	 if(genparticles->at(bestjetindex).pt()>600 && genparticles->at(bestjetindex).pt()<800)	Hist("ptfraction1_600")->Fill( SortedSubJets.at(0).pt()/hotvr_jets[i].pt(),weight);
	 if(genparticles->at(bestjetindex).pt()>800 && genparticles->at(bestjetindex).pt()<1000)	Hist("ptfraction1_800")->Fill( SortedSubJets.at(0).pt()/hotvr_jets[i].pt(),weight);
	 if(genparticles->at(bestjetindex).pt()>1000 && genparticles->at(bestjetindex).pt()<1200)	Hist("ptfraction1_1000")->Fill( SortedSubJets.at(0).pt()/hotvr_jets[i].pt(),weight);

	 //Fill nominator hists
	 if(SortedSubJets.size()>2)  Hist("pT_s1_hotvr_tagged")->Fill(genparticles->at(bestjetindex).pt(),weight);
	 if(hotvr_jets[i].m()>mtopLow &&  hotvr_jets[i].m()<mtopHigh)  Hist("pT_s1_hotvr2_tagged")->Fill(genparticles->at(bestjetindex).pt(),weight);
	 if(SortedSubJets.size()>2 && hotvr_jets[i].m()>mtopLow &&  hotvr_jets[i].m()<mtopHigh)  Hist("pT_s1_hotvr3_tagged")->Fill(genparticles->at(bestjetindex).pt(),weight);
	 if(SortedSubJets.size()>2 && hotvr_jets[i].m()>mtopLow &&  hotvr_jets[i].m()<mtopHigh    &&SortedSubJets.at(0).pt()>30 &&SortedSubJets.at(1).pt()>30&&  SortedSubJets.at(0).pt()/hotvr_jets[i].pt()<0.8 /*&&mmin2>50*/)  Hist("pT_s1_hotvr6_tagged")->Fill(genparticles->at(bestjetindex).pt(),weight);
	 if(SortedSubJets.size()>2 && hotvr_jets[i].m()>mtopLow &&  hotvr_jets[i].m()<mtopHigh && mmin>50 && hotvr_jets.size()>2)  Hist("pT_s1_hotvr4_tagged")->Fill(genparticles->at(bestjetindex).pt(),weight);
	
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
	   if(genparticles->at(bestjetindex).pt()>200) Hist("tau1_200")->Fill(tau1,weight);
	   if(genparticles->at(bestjetindex).pt()>500) Hist("tau1_500")->Fill(tau1,weight);
	   if(genparticles->at(bestjetindex).pt()>600) Hist("tau1_600")->Fill(tau1,weight);
	   if(genparticles->at(bestjetindex).pt()>800) Hist("tau1_800")->Fill(tau1,weight);
	   if(genparticles->at(bestjetindex).pt()>1000) Hist("tau1_1000")->Fill(tau1,weight);
	   Hist("tau2")->Fill(tau2,weight);
	   if(genparticles->at(bestjetindex).pt()>200) Hist("tau2_200")->Fill(tau2,weight);
	   if(genparticles->at(bestjetindex).pt()>500) Hist("tau2_500")->Fill(tau2,weight);
	   if(genparticles->at(bestjetindex).pt()>600) Hist("tau2_600")->Fill(tau2,weight);
	   if(genparticles->at(bestjetindex).pt()>800) Hist("tau2_800")->Fill(tau2,weight);
	   if(genparticles->at(bestjetindex).pt()>1000) Hist("tau2_1000")->Fill(tau2,weight);
	   Hist("tau3")->Fill(tau3,weight);
	   if(genparticles->at(bestjetindex).pt()>200) Hist("tau3_200")->Fill(tau3,weight);
	   if(genparticles->at(bestjetindex).pt()>500) Hist("tau3_500")->Fill(tau3,weight);
	   if(genparticles->at(bestjetindex).pt()>600) Hist("tau3_600")->Fill(tau3,weight);
	   if(genparticles->at(bestjetindex).pt()>800) Hist("tau3_800")->Fill(tau3,weight);
	   if(genparticles->at(bestjetindex).pt()>1000) Hist("tau3_1000")->Fill(tau3,weight);
	   Hist("tau4")->Fill(tau4,weight);
	   if(genparticles->at(bestjetindex).pt()>200) Hist("tau4_200")->Fill(tau4,weight);
	   if(genparticles->at(bestjetindex).pt()>500) Hist("tau4_500")->Fill(tau4,weight);
	   if(genparticles->at(bestjetindex).pt()>600) Hist("tau4_600")->Fill(tau4,weight);
	   if(genparticles->at(bestjetindex).pt()>800) Hist("tau4_800")->Fill(tau4,weight);
	   if(genparticles->at(bestjetindex).pt()>1000) Hist("tau4_1000")->Fill(tau4,weight);
	   Hist("tau2tau1")->Fill(tau2/tau1,weight);
	   if(genparticles->at(bestjetindex).pt()>200)  Hist("tau2tau1_200")->Fill(tau2/tau1,weight);
	   if(genparticles->at(bestjetindex).pt()>500)  Hist("tau2tau1_500")->Fill(tau2/tau1,weight);
	   if(genparticles->at(bestjetindex).pt()>600)  Hist("tau2tau1_600")->Fill(tau2/tau1,weight);
	   if(genparticles->at(bestjetindex).pt()>800)  Hist("tau2tau1_800")->Fill(tau2/tau1,weight);
	   if(genparticles->at(bestjetindex).pt()>1000)  Hist("tau2tau1_1000")->Fill(tau2/tau1,weight);
	   Hist("tau3tau2")->Fill(tau3/tau2,weight);
	   if(genparticles->at(bestjetindex).pt()>200)  Hist("tau3tau2_200")->Fill(tau3/tau2,weight);
	   if(genparticles->at(bestjetindex).pt()>500)  Hist("tau3tau2_500")->Fill(tau3/tau2,weight);
	   if(genparticles->at(bestjetindex).pt()>600)  Hist("tau3tau2_600")->Fill(tau3/tau2,weight);
	   if(genparticles->at(bestjetindex).pt()>800)  Hist("tau3tau2_800")->Fill(tau3/tau2,weight);
	   if(genparticles->at(bestjetindex).pt()>1000)  Hist("tau3tau2_1000")->Fill(tau3/tau2,weight);
	 
	 //Fill more nominator hists
	
	   if(SortedSubJets.size()>2 &&  hotvr_jets[i].m()>mtopLow &&   hotvr_jets[i].m()<mtopHigh &&tau3/tau2<0.7  &&SortedSubJets.at(0).pt()>30 &&SortedSubJets.at(1).pt()>30&&  SortedSubJets.at(0).pt()/hotvr_jets[i].pt()<0.8 &&mmin>50)  Hist("pT_s1_hotvr5_tagged")->Fill(genparticles->at(bestjetindex).pt(),weight);
	 if(SortedSubJets.size()>2 &&  hotvr_jets[i].m()>mtopLow &&   hotvr_jets[i].m()<mtopHigh   &&SortedSubJets.at(0).pt()>30 &&SortedSubJets.at(1).pt()>30&&  SortedSubJets.at(0).pt()/hotvr_jets[i].pt()<0.8 &&mmin>50) for(int roc=0;roc<100;roc++){
	     if(tau3/tau2<roc/100.){
	       ((TH2D*)Hist("nsub_eff"))->Fill(genparticles->at(bestjetindex).pt(),roc/100.,weight);
	    
	     }}
	   }
	 
       }
	}
    }
  }
  

    
    

  
  
  
  


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

