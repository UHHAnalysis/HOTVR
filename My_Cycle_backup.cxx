// $Id: MyCycle.cxx,v 1.10 2012/12/07 14:21:51 peiffer Exp $

#include <iostream>

using namespace std;

// Local include(s):
#include "include/MyCycle.h"

ClassImp( MyCycle );

MyCycle::MyCycle()
   : AnalysisCycle() {

  // constructor, declare additional variables that should be 
  // obtained from the steering-xml file
  
  // set the integrated luminosity per bin for the lumi-yield control plots
  SetIntLumiPerBin(500.);

}

MyCycle::~MyCycle() 
{
  // destructor
}

void MyCycle::BeginCycle() throw( SError ) 
{
  // Start of the job, general set-up and definition of 
  // objects are done here

  // Important: first call BeginCycle of base class
  AnalysisCycle::BeginCycle();

  return;

}

void MyCycle::EndCycle() throw( SError ) 
{
  // clean-up, info messages and final calculations after the analysis

  
  // call the base cycle class for all standard methods
  AnalysisCycle::EndCycle();

  return;

}

void MyCycle::BeginInputData( const SInputData& id ) throw( SError ) 
{
  
  // declaration of histograms and selections

  // Important: first call BeginInputData of base class
  AnalysisCycle::BeginInputData( id );

  // -------------------- set up the selections ---------------------------

  Selection* BSel = new Selection( "BSelection");
  BSel->addSelectionModule(new NBTagSelection(1)); //at least one b tag

  Selection* NoBSel = new Selection( "NoBSelection");
  NoBSel->addSelectionModule(new NBTagSelection(0,0)); //no b tags

  Selection* chi2_selection= new Selection("chi2_selection");
  static Chi2Discriminator* m_chi2discr = new Chi2Discriminator();
  chi2_selection->addSelectionModule(new HypothesisDiscriminatorCut( m_chi2discr, -1*double_infinity(), 10));
  chi2_selection->addSelectionModule(new MttbarGenCut(0,700));

  Selection* TopSel = new Selection("TopSelection");
  //DO NOT use trigger selection in PROOF mode at the moment
  //TopSel->addSelectionModule(new TriggerSelection("HLT_PFJet320_v"));
  TopSel->addSelectionModule(new NTopJetSelection(1,int_infinity(),350,2.5));
  TopSel->addSelectionModule(new NTopTagSelection(1,int_infinity()));




  RegisterSelection(BSel);
  RegisterSelection(NoBSel);
  RegisterSelection(TopSel);
  RegisterSelection(chi2_selection);

  // ---------------- set up the histogram collections --------------------

  // histograms without any cuts
  RegisterHistCollection( new MyHists("NoCuts") );
  RegisterHistCollection( new HypothesisHists("Chi2_NoCuts", m_chi2discr ) );

  //histograms with and without b tagging
  RegisterHistCollection( new MyHists("BTag") );
  RegisterHistCollection( new MyHists("NoBTag") );
  RegisterHistCollection( new HypothesisHists("Chi2_BTag", m_chi2discr ) );
  RegisterHistCollection( new HypothesisHists("Chi2_NoBTag", m_chi2discr ) );

  // histograms after the top selection
  RegisterHistCollection( new MyHists("TopSel") );
  RegisterHistCollection( new HypothesisHists("Chi2_TopSel", m_chi2discr ) );



  // important: initialise histogram collections after their definition

  //test Histos
   Book( TH1F( "Nmuon","number of muons",7,-0.5,6.5));
   Book(TH1F("muon_pt","pt muon",1000,0,1000));
   Book(TH1F("chargemuon","charge of muon",5,-2,2));
   Book(TH1F("invarmass","invariant mass",1000,0,5000));
   Book(TH1F("MET","missing Et",1000,0,5000));
   Book(TH1F("met_hadr","missing Et hadr",1000,0,5000));
   Book(TH1F("met_lept","missing Et lept",1000,0,5000));
   Book(TH1F("Ht","Ht",1000,0,5000));
   Book(TH1F("pt_ca_jets","pt ca jets",50,0,1000)); 
   Book(TH1F("pt_ca_jets2","pt ca jets",50,0,1000));
    Book(TH1F("pt_ca_jets3","pt ca jets",50,0,1000));
   Book(TH1F("pt_ca_jet1","pt ca jet1",50,0,1000));
   Book(TH1F("pt_ca_jet2","pt ca jet2",50,0,1000));
   Book(TH1F("pt_ca_jet3","pt ca jet3",50,0,1000));
   Book(TH1F("pt_ca_jet4","pt ca jet4",50,0,1000));
Book(TH1F("controlplot","controlplot",4,0,3));
   Book(TH1F("flavor_ca_jets","flavor ca jets",3,-1.5,1.5));
   Book(TH1F("flavor_ca_jet1","flavor ca jet1",4,0,3));

 Book(TH1F("flavor_ca_jet2","flavor ca jet2",4,0,3));
Book(TH1F("flavor_ca_jet3","flavor ca jet3",4,0,3));
Book(TH1F("flavor_ca_jet4","flavor ca jet4",4,0,3));
Book(TH1F("pt_ca_jets_tagged_cms","pt ca jets tagged cms",50,0,1000));
   Book(TH1F("pt_ca_jets_tagged","pt ca jets tagged",50,0,1000));
   Book(TH1F("pt_ca_jets_tagged2","pt ca jets tagged2",50,0,1000));
   Book(TH1F("pt_ca_lept","pt of leptonic ca",50,0,1000));
   Book(TH1F("pt_ca_hadr","pt of hadronic ca",50,0,1000));
   Book(TH1F("M_ca_lept","mass of leptonic ca",50,0,500));
   Book(TH1F("M_ca_hadr","mass of hadronic ca",50,0,500));
   Book(TH1F("NTopJets","N top-jets",1000,0,20));
   Book(TH1F("NJets","number of Jets",101,0,100));
   Book(TH1F("Nsubjets","N subjets",11,0,10));
   Book(TH1F("Nsubjets2","N subjets2",11,0,10));
   Book(TH1F("Q_subjet1","Q subjet1",1000,-10,10));
   Book(TH1F("Msubjets","Mass subjets",1000,0,1000));
   Book(TH1F("invarsub12","ivariant mass of subjet 1 and 2",1000,0,1000));
   Book(TH1F("invarsub12x","invariant mass of subjet 1 and 2",1000,0,1000));
 Book(TH1F("invarsub13","invariant mass of subjet 1 and 3",1000,0,1000));
 Book(TH1F("invarsub23","invariant mass of subjet 2 and 3",1000,0,1000));
 Book(TH1F("genparticles","particles generated",101,0,100));
 Book(TH1F("pttopgen","pt top generated",50,0,1000));
 Book(TH1F("particle_caused_jet_had","particles that caused jet hadr",3,-1.5,1.5));
 Book(TH1F("particle_caused_jet_lep","particles that caused jet lept",3,-1.5,1.5));
 Book(TH1F("particle_caused_jet_all","particles that caused jet all",3,-1.5,1.5));
 Book(TH1F("topdaughter1","daughter1 of W(top)",101,-50,50));
 Book(TH1F("topdaughter2","daughter2 of W(top)",101,-50,50));
 Book(TH1F("ntoptags","n toptags",1001,0,1000));
 Book(TH1F("ntopgen", "number of generated tops (hadr)",1001,0,1000));
 Book(TH1F("Mtopjet","mass of topjet",1001,0,1000));
 Book(TH1F("Mtopjet_sub","mass of topjet(just 3subs)",1001,0,1000));
 Book(TH1F("ptjettag","pt tagged jet",50,0,1000));
 Book(TH1F("pttoptag","pt tagged alg",50,0,1000));
 Book(TH1F("etajettag","eta jet tagged",50,-5,5));
Book(TH1F("etajettagpt","eta jet tagged",50,-5,5));
 Book(TH1F("etatopgen","eta top generated (hadr)",50,-5,5));
Book(TH1F("etatopgenpt","eta top generated",50,-5,5));
 Book(TH1F("phi","phi jet",300,-4,4));
 Book(TH1F("distance","distance between generated and tagged",110,0,10));
  Book(TH1F("distance_s1_s2","distance between subjet1 and 2",110,0,10));
   Book(TH1F("distance_s1_s3","distance between subjet1 and 2",110,0,10));
    Book(TH1F("distance_s2_s3","distance between subjet1 and 2",110,0,10));
    Book(TH1F("pt_jet_div_subj1","pt_jet/pt_subj1",100,0,5));
    Book(TH1F("match_ca_jet","match ca jet",5,-1,3));
 Book(TH1F("match_ca_jet2","match ca jet",5,-1,3));
 Book(TH1F("match","pt match generated tagged",50,0,1000));
   Book(TH1F("GenTop","number of generated tops",1000,0,1000));\
   Book(TH1F("tageff","eff vs. pt",50,0,1000));
   Book(TH1F("tageff2","eff vs. pt gegeb",50,0,1000));
   Book(TH1F("tageffeta","eff vs. eta",50,-5,5));
 Book(TH1F("tageffetapt","eff vs. eta (pt>650)",50,-5,5));
Book(TH1F("tageff_ca","eff vs. pt (ca jets)",50,0,1000));
 Book(TH1F("tageff_ca2","eff vs. pt2 (ca jets)",50,0,1000));
  InitHistos();

  return;

}

void MyCycle::EndInputData( const SInputData& id ) throw( SError ) 
{
  AnalysisCycle::EndInputData( id );
  return;

}

void MyCycle::BeginInputFile( const SInputData& id ) throw( SError ) 
{
  // Connect all variables from the Ntuple file with the ones needed for the analysis
  // The variables are commonly stored in the BaseCycleContaincer

  // important: call to base function to connect all variables to Ntuples from the input tree
  AnalysisCycle::BeginInputFile( id );

  return;

}

void MyCycle::ExecuteEvent( const SInputData& id, Double_t weight) throw( SError ) 
{
   // this is the most important part: here the full analysis happens
  // user should implement selections, filling of histograms and results

  // first step: call Execute event of base class to perform basic consistency checks
  // also, the good-run selection is performed there and the calculator is reset
  AnalysisCycle::ExecuteEvent( id, weight );
  /* EventCalc* calc = EventCalc::Instance();
  double weight2 = calc -> GetWeight();

    BaseCycleContainer* bcc = calc->GetBaseCycleContainer();

    double tobias,tobias2;
  tobias2=bcc->electrons->size();
  //Electron electron =  bcc->electrons->at(1);
  tobias=electron.pt();
  Hist("number")-> Fill(tobias2,weight2);*/
  

  // get the selections
  static Selection* BSel = GetSelection("BSelection");
  static Selection* NoBSel = GetSelection("NoBSelection");
  static Selection* TopSel = GetSelection("TopSelection");
  static Selection* chi2_selection = GetSelection("chi2_selection");

  
  // get the histogram collections
  BaseHists* HistsNoCuts = GetHistCollection("NoCuts");
  BaseHists* HistsBTag = GetHistCollection("BTag");
  BaseHists* HistsNoBTag = GetHistCollection("NoBTag");
  BaseHists* HistsTopSel = GetHistCollection("TopSel");

  
  BaseHists* Chi2_HistsNoCuts = GetHistCollection("Chi2_NoCuts");
  BaseHists* Chi2_HistsBTag = GetHistCollection("Chi2_BTag");
  BaseHists* Chi2_HistsNoBTag = GetHistCollection("Chi2_NoBTag");
  BaseHists* Chi2_HistsTopSel = GetHistCollection("Chi2_TopSel");
  
  //get test histos
  TH1* genpt = Hist("pttopgen");
  TH1* tagpt = Hist("ptjettag");
  TH1* bla3 = Hist("pttoptag");
  TH1* matched= Hist("match");
  TH1* geneta = Hist("etatopgen");
  TH1* tageta = Hist("etajettag");
  TH1* genetapt = Hist("etatopgenpt");
  TH1* tagetapt = Hist("etajettagpt");
  TH1* pt_ca = Hist("pt_ca_jets");
 TH1* pt_ca2 = Hist("pt_ca_jets2");
  TH1* pt_ca_tagged = Hist("pt_ca_jets_tagged");
   TH1* pt_ca_tagged2 = Hist("pt_ca_jets_tagged2");
   TH1* pt_ca_tagged_cms = Hist("pt_ca_jets_tagged_cms");

  //if(!chi2_selection->passSelection())  throw SError( SError::SkipEvent );

  EventCalc* calc = EventCalc::Instance();
  double weight2 = calc->GetWeight();
   BaseCycleContainer* bcc = calc->GetBaseCycleContainer();
     MET met   = * bcc->met;
   double missingE;
    missingE=met.pt();
    Hist("MET")->Fill(missingE);
   Hist("genparticles")->Fill(calc->GetGenParticles()->size(),weight2);
   Hist("NTopJets")->Fill(bcc->topjets->size(),weight2);
   Hist("NJets")->Fill(bcc->jets->size(),weight2);
   // double weight = calc->GetWeight();

   Hist("Nmuon")->Fill(bcc->muons->size(),weight2);
   /*  for(int t=0;t<bcc->muons->size();++t){
     //sort(bcc->muons.begin(),bcc->muons.end(),HigherPt());
     Muon myon = bcc->muons->at(t);
     Hist("muon_pt")->Fill(myon.pt());
     for(o=0;o<bcc->topjets->size();++o){
       if(sqrt(pow(topjet.phi()-top.phi(),2)+pow(topjet.eta()-top.eta(),2))<0.8) 
     }

     
     }*/



   // TTbarGen* ttbargen= calc->GetTTbarGen();
   GenParticle top;
   for(unsigned int k=0;k<calc->GetGenParticles()->size();++k)//matching generated tagged?
     if((calc->GetGenParticles()->at(k).pdgId())==6) {
       top=calc->GetGenParticles()->at(k);
   Hist("topdaughter1")->Fill(top.daughter(calc->GetGenParticles(),2)->daughter(calc->GetGenParticles(),1)->pdgId(),weight2);//hadronisch oder leptonisch?
   Hist("topdaughter2")->Fill(top.daughter(calc->GetGenParticles(),2)->daughter(calc->GetGenParticles(),2)->pdgId(),weight2);
   if((abs(top.daughter(calc->GetGenParticles(),2)->daughter(calc->GetGenParticles(),1)->pdgId())<7) ) Hist("flavor_ca_jets")->Fill(1);
   else Hist("flavor_ca_jets")->Fill(-1);
   if((abs(top.daughter(calc->GetGenParticles(),2)->daughter(calc->GetGenParticles(),1)->pdgId())<7) ) Hist("met_hadr")->Fill(met.pt());
   else Hist("met_lept")->Fill(met.pt());
   if(abs(top.daughter(calc->GetGenParticles(),2)->daughter(calc->GetGenParticles(),1)->pdgId())<7) {//if hadronic
     Hist("pttopgen")->Fill(top.v4().Pt(),weight2);
     Hist("etatopgen")->Fill(top.v4().eta(),weight2);
     
     if(top.v4().Pt()>650) Hist("etatopgenpt")->Fill(top.v4().eta(),weight2);
   
   }
   }
  
   int ntop=0;
   for(unsigned int k=0;k<calc->GetGenParticles()->size();++k)  if(calc->GetGenParticles()->at(k).pdgId()==6) ntop++;
   if(ntop>0) Hist("ntopgen")->Fill(ntop,weight2);
   //cout<<ntop<<endl;
   int ntoptags = 0;
   std::vector<TopJet> blax(10);
   for(unsigned int t=0; t<bcc->topjets->size(); ++t) {
      TopJet topjet2 = bcc->topjets->at(t);
     blax.at(t)=topjet2;
     
   }
   sort(blax.begin(),blax.end(),HigherPt());
   Hist("pt_ca_jet1")->Fill(blax[0].pt(),weight2); 
   Hist("pt_ca_jet2")->Fill(blax[1].pt(),weight2); 
   Hist("pt_ca_jet3")->Fill(blax[2].pt(),weight2); 
   Hist("pt_ca_jet4")->Fill(blax[3].pt(),weight2); 
    std::vector<Muon> muon(bcc->muons->size());
     for(unsigned int t=0; t<bcc->muons->size();++t) muon.at(t)=bcc->muons->at(t);
     sort(muon.begin(),muon.end(),HigherPt());
     Hist("muon_pt")->Fill(muon[0].pt());
   for(unsigned int t=0; t<bcc->topjets->size();++t){
     
     
   }


   for(unsigned int j=0 ;j<bcc->topjets->size();++j) {
     TopJet topjet = blax.at(j);
     //if(abs(top.daughter(calc->GetGenParticles(),2)->daughter(calc->GetGenParticles(),1)->pdgId())<7) 
      
     bool hadronic_topjet=true;
    
     if(sqrt(pow(topjet.phi()-muon[0].phi(),2)+pow(topjet.eta()-muon[0].eta(),2))<0.8 ) {
       Hist("controlplot")->Fill(2);
        hadronic_topjet=false;
     }
     else {
       Hist("controlplot")->Fill(1);
      
      
       }

     /* hadronic_topjet=true;
      for(unsigned int g=0;g<bcc->jets->size();++g){
       Jet jet = bcc->jets->at(g);
       if(sqrt(pow(topjet.phi()-jet.phi(),2)+pow(topjet.eta()-jet.eta(),2))>0.5) Hist("pt_ca_jets")->Fill(topjet.pt(),weight2);
       else hadronic_topjet=false;
       }*/
    


     bool check=false;
   for(unsigned int k=0;k<calc->GetGenParticles()->size();++k){
	GenParticle top=calc->GetGenParticles()->at(k);
	
	if(abs(calc->GetGenParticles()->at(k).pdgId())==6) {
	  if(sqrt(pow(topjet.phi()-top.phi(),2)+pow(topjet.eta()-top.eta(),2))<0.8){
	    check=true;
	 
	    if((abs(top.daughter(calc->GetGenParticles(),2)->daughter(calc->GetGenParticles(),1)->pdgId())<7)&& hadronic_topjet ) Hist("particle_caused_jet_had")->Fill(1);
	    else if(hadronic_topjet) Hist("particle_caused_jet_had")->Fill(-1);
	    if((abs(top.daughter(calc->GetGenParticles(),2)->daughter(calc->GetGenParticles(),1)->pdgId())<7) && !hadronic_topjet) Hist("particle_caused_jet_lep")->Fill(1);
	    else if(!hadronic_topjet) Hist("particle_caused_jet_lep")->Fill(-1);
	    Hist("pt_ca_jets3")->Fill(topjet.pt(),weight2);
	    if(hadronic_topjet){
	      switch (j){
	      case 0: Hist("match_ca_jet")->Fill(0); break;
	      case 1: Hist("match_ca_jet")->Fill(1); break;
	      case 2: Hist("match_ca_jet")->Fill(2); break;
	      case 3: Hist("match_ca_jet")->Fill(3); break;
	      }
	    }
	    if((abs(top.daughter(calc->GetGenParticles(),2)->daughter(calc->GetGenParticles(),1)->pdgId())<7)) Hist("particle_caused_jet_all")->Fill(1);
	    else Hist("particle_caused_jet_all")->Fill(-1);
	}
	
	}

  
   }
   if(hadronic_topjet & check){ Hist("pt_ca_jets")->Fill(topjet.pt(),weight2);}
      if (check)   Hist("pt_ca_jets2")->Fill(topjet.pt(),weight2);
   // if(!check) hadronic_topjet=false;
   	if(hadronic_topjet){
	      switch (j){
	      case 0: Hist("match_ca_jet2")->Fill(0); break;
	      case 1: Hist("match_ca_jet2")->Fill(1); break;
	      case 2: Hist("match_ca_jet2")->Fill(2); break;
	      case 3: Hist("match_ca_jet2")->Fill(3); break;
	      }
	    }

      



       
    
     //Hist("pt_ca_jets")->Fill(topjet.pt(),weight2);
     std::vector<Particle> subjet = topjet.subjets();
      sort(subjet.begin(),subjet.end(),HigherPt());
      double nsubjets=topjet.numberOfDaughters();
      Hist("Nsubjets")->Fill(nsubjets,weight2);
     
      if(nsubjets>=3){//min 3 subjets in topjet
	//if((subjet[0].v4())!=0); //Hist("pt_jet_div_subj1")->Fill(topjet.v4().M()/subjet[0].v4().M(),weight2);
	// if(IsTagged(topjet,e_JPL)){
    
    
	
	double invarmasssubj12=0;
	invarmasssubj12=(subjet[0].v4()+subjet[1].v4()).M();
	double invarmasssubj13=0;
	invarmasssubj13=(subjet[0].v4()+subjet[2].v4()).M();
	double invarmasssubj23=0;
	invarmasssubj23=(subjet[1].v4()+subjet[2].v4()).M();
	// invarmasssubj12= sqrt(2*subjet[0].pt()*subjet[1].pt()*(cosh(subjet[0].eta()-subjet[1].eta())-cos(subjet[0].phi()-subjet[1].phi())));
	//double invarmasssubj13= sqrt(2*subjet[0].pt()*subjet[2].pt()*(cosh(subjet[0].eta()-subjet[2].eta())-cos(subjet[0].phi()-subjet[2].phi())));  
	//double invarmasssubj23= sqrt(2*subjet[1].pt()*subjet[2].pt()*(cosh(subjet[2].eta()-subjet[2].eta())-cos(subjet[1].phi()-subjet[2].phi())));
	double minsubjet = min(invarmasssubj12,min(invarmasssubj13,invarmasssubj23));
      if(minsubjet>50){ //invariant mass greater than 50Gev?
	Hist("invarsub12")->Fill(invarmasssubj12,weight2);
	Hist("invarsub13")->Fill(invarmasssubj13,weight2);
	Hist("invarsub23")->Fill(invarmasssubj23,weight2);
	Hist("Q_subjet1")->Fill(subjet[2].charge(),weight2);
	Hist("distance_s1_s2")->Fill(sqrt(pow(subjet[1].phi()-subjet[2].phi(),2)+pow(subjet[1].eta()-subjet[2].eta(),2)));
	Hist("distance_s1_s3")->Fill(sqrt(pow(subjet[1].phi()-subjet[3].phi(),2)+pow(subjet[1].eta()-subjet[3].eta(),2)));
	Hist("distance_s2_s3")->Fill(sqrt(pow(subjet[2].phi()-subjet[3].phi(),2)+pow(subjet[2].eta()-subjet[3].eta(),2)));
	double topjetmass=(subjet[0].v4()+subjet[1].v4()+subjet[2].v4()).M();
	if(topjetmass<250 && topjetmass>140){
	Hist("ptjettag")->Fill(topjet.pt(),weight2);//jet mass in right region?
	Hist("Mtopjet")->Fill(topjet.v4().M(),weight2);
	Hist("Mtopjet_sub")->Fill(topjetmass,weight2);
	Hist("phi")->Fill(topjet.phi(),weight2);      
	if(check && hadronic_topjet) Hist("pt_ca_jets_tagged")->Fill(topjet.pt(),weight2);
	if(check) Hist("pt_ca_jets_tagged2")->Fill(topjet.pt(),weight2);
	      for(unsigned int k=0;k<calc->GetGenParticles()->size();++k)//matching generated tagged?
	      if(abs(calc->GetGenParticles()->at(k).pdgId())==6) {
	       GenParticle top=calc->GetGenParticles()->at(k); 
	       double R= sqrt(pow(topjet.phi()-top.phi(),2)+pow(topjet.eta()-top.eta(),2));
	       Hist("distance")->Fill(R);
	    if(R<0.8){ 
	      switch(j){
		case 0:
		  if((abs(top.daughter(calc->GetGenParticles(),2)->daughter(calc->GetGenParticles(),1)->pdgId())<7) ) Hist("flavor_ca_jet1")->Fill(1);
	      else Hist("flavor_ca_jet1")->Fill(2);
	      break;
	    case 1:
	      if((abs(top.daughter(calc->GetGenParticles(),2)->daughter(calc->GetGenParticles(),1)->pdgId())<7) ) Hist("flavor_ca_jet2")->Fill(1);
	      else Hist("flavor_ca_jet2")->Fill(2);
	      break;
	    case 2:
	      if((abs(top.daughter(calc->GetGenParticles(),2)->daughter(calc->GetGenParticles(),1)->pdgId())<7) ) Hist("flavor_ca_jet3")->Fill(1);
	      else Hist("flavor_ca_jet3")->Fill(2);
	      case 3:
		if((abs(top.daughter(calc->GetGenParticles(),2)->daughter(calc->GetGenParticles(),1)->pdgId())<7) ) Hist("flavor_ca_jet4")->Fill(1);
	      else Hist("flavor_ca_jet4")->Fill(2);
	      break;}
	       if(top.daughter(calc->GetGenParticles(),2)->daughter(calc->GetGenParticles(),1)->pdgId()<7) Hist("pt_ca_hadr")->Fill(topjet.pt(),weight2);
	       else Hist("pt_ca_lept")->Fill(topjet.pt(),weight2);
	       if(top.daughter(calc->GetGenParticles(),2)->daughter(calc->GetGenParticles(),1)->pdgId()<7) Hist("M_ca_hadr")->Fill(topjet.v4().M(),weight2);
	       else Hist("M_ca_lept")->Fill(topjet.v4().M(),weight2);
	      Hist("match")->Fill(top.pt(),weight2);
	      Hist("etajettag")->Fill(top.eta(),weight2);
	     
	      HistsTopSel->Fill();
	      if(top.pt()>650) Hist("etajettagpt")->Fill(top.eta(),weight2);
	    }
	// }
	  }
      }}}
     

      double mmin=0;
     double mjet=0;
     int nsubjets2=0;
     if(TopTag(topjet,mjet,nsubjets2,mmin)){
	for(unsigned int k=0;k<calc->GetGenParticles()->size();++k)//matching generated tagged?
	  if((calc->GetGenParticles()->at(k).pdgId())==6) {
	    GenParticle top=calc->GetGenParticles()->at(k); 
	    double R2= sqrt(pow(topjet.phi()-top.phi(),2)+pow(topjet.eta()-top.eta(),2));
	    if(R2<0.8 && hadronic_topjet) Hist("pt_ca_jets_tagged_cms")->Fill(topjet.pt(),weight2);
	  }
       ntoptags++;
       //Hist("invarsub12x")->Fill(sqrt(2*subjet[1].pt()*subjet[2].pt()*(cosh(subjet[1].eta()-subjet[2].eta())-cos(subjet[1].phi()-subjet[2].phi()))));
       Hist("invarsub12x")->Fill((subjet[0].v4()+subjet[1].v4()).M());
       Hist("pttoptag")->Fill(topjet.pt());
	  }

   
    
   }
  
   Hist("ntoptags")->Fill(ntoptags);
     // Hist("tageff")->Divide(tagpt,genpt);
     Hist("tageff")->Divide(matched,genpt);
     Hist("tageff2")->Divide(bla3,genpt);
     Hist("tageffeta")->Divide(tageta,geneta);
      Hist("tageffetapt")->Divide(tagetapt,genetapt);
      Hist("tageff_ca")->Divide(pt_ca_tagged,pt_ca);
       Hist("tageff_ca2")->Divide(pt_ca_tagged2,pt_ca2);

   /*  double tobias,tobias2,missingE,numbertopjets;
  for(unsigned int j=0; j<bcc->topjets->size(); ++j){
     TopJet topjets= bcc->topjets->at(j);
     double sub=topjets.numberOfDaughters();
     Hist("Nsubjets")->Fill(sub);
   }
   tobias2=bcc->muons->size();
     numbertopjets=bcc->topjets->size();
     Hist("NTopJets")->Fill(numbertopjets);
      for(unsigned int i=0; i< bcc->muons->size(); ++i) {
	//	cout<<i<<endl;

    double ht;
    ht=calc->GetHT();
    Hist("Ht")->Fill(ht);
    //double nprimaryVertices;
    //std::vector< PrimaryVertex > primaryVertex = calc->GetPrimaryVertices();
    // nprimaryVertices=primaryVertex.nTracks();
    //Hist("PrimaryVertices")->Fill(nprimaryVertices);
    
   
}
    if(bcc->muons->size()>1){
    Muon muon1 = bcc->muons->at(0);
    Muon muon2 = bcc->muons->at(1);
    if(muon1.charge()+muon2.charge()==0){
      double invarmass =sqrt(2*muon1.pt()*muon2.pt()*(cosh(muon1.eta()-muon2.eta())-cos(muon1.phi()-muon2.phi())));
    Hist("invarmass")->Fill(invarmass);}
    }
    Hist("number")-> Fill(tobias2,weight);
   */
  // start the analysis
   HistsNoCuts->Fill();
   //  Chi2_HistsNoCuts->Fill();
      
  if(BSel->passSelection()){
    HistsBTag->Fill();
    // Chi2_HistsBTag->Fill();
  }
  if(NoBSel->passSelection()){
    HistsNoBTag->Fill();  
    // Chi2_HistsNoBTag->Fill();
  }
  

  if(!TopSel->passSelection())  throw SError( SError::SkipEvent );
  
  
  //Chi2_HistsTopSel->Fill();
  
  return;
  
}


