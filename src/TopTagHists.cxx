#include "include/EventCalc.h"
#include "include/TopTagHists.h"
#include "include/TopFitCalc.h"
#include <iostream>
#include "include/TopTagfunctions.h"

using namespace std;

TopTagHists::TopTagHists(const char* name) : BaseHists(name)
{
  // named default constructor
  
}

TopTagHists::~TopTagHists()
{
  // default destructor, does nothing
}

void TopTagHists::Init()
{
  // book all histograms here
  Book( TH1F( "NTopJets", "number of topjets", 8, -0.5, 7.5 ) );
  Book( TH1F( "NMuons", "number of muons", 8, -0.5, 7.5 ) );
  Book( TH1F( "NElectrons", "number of electrons", 8, -0.5, 7.5 ) );
  Book( TH1F( "NTaus","number of taus", 8, -0.5, 7.5 ) );
  Book( TH1F( "Nbtags","number of b-Tags", 8, -0.5, 7.5 ) );
Book( TH1F( "NPFflow","number of particle flow", 1001, -0.5, 1000.5 ) );
Book( TH1F( "NPFflow2","number of particle flow", 1001, -0.5, 1000.5 ) );
    Book( TH1F( "MJet_out", "m_{jet} outside2", 100,0,2000 ) );
      Book( TH1F( "MJet_out2", "m_{jet} outside2", 100,0,2000 ) );
Book( TH1F( "jetclosetop","number of jets close to topjet", 8, -0.5, 7.5 ) );
 Book( TH1F( "jetnotclosetop","number of jets not close to topjet", 8, -0.5, 7.5 ) );
 Book( TH1F( "jetclosetop2","number of jets close to topjet", 8, -0.5, 7.5 ) );
 Book( TH1F( "jetnotclosetop2","number of jets not close to topjet", 8, -0.5, 7.5 ) );
  Book( TH1F( "NRecoHyps", "number of recohyps", 20, -0.5, 19.5 ) );
  Book( TH1F( "NTopJets_after_selection", "number of topjets  after selection", 8, -0.5, 7.5 ) );
  Book( TH1F( "NTopJets_ly", "number of topjets", 8, -0.5, 7.5 ) );
  Book( TH1F( "flavor_topjet", "decay of topjet", 8, -0.5 ,7.5 ) );
   Book( TH1F( "decay_products0", "decay products in jet (0<p_{T}<400)", 8, -0.5 ,7.5 ) );
  Book( TH1F( "decay_products200", "decay products in jet (400<p_{T}<700)", 8, -0.5 ,7.5 ) );
  Book( TH1F( "decay_products400", "decay products in jet (p_{T}>700)", 8, -0.5 ,7.5 ) );
  /* Book( TH1F( "decay_products600", "decay products in jet (600>p_{T}>800)", 8, -0.5 ,7.5 ) );
   Book( TH1F( "decay_products800", "decay products in jet (800>p_{T}>1000)", 8, -0.5 ,7.5 ) );
   Book( TH1F( "decay_products1000", "decay products in jet (p_{T}>1000)", 8, -0.5 ,7.5 ) );*/
  Book(TH2D("nusbjetsvspt","Nsubjets vs. pt",10,0,2000,6,-0.5,5.5));
  Book( TH1F( "flavor_topjet_after_selection", "decay of topjet after selection", 8, -0.5 ,7.5 ) );
	Book( TH1F( "flavor_topjet_after_selection2", "decay of topjet after selection2", 8, -0.5 ,7.5 ) );
	Book( TH1F( "flavor_topjet_after_selection3", "decay of topjet after selection3", 8, -0.5 ,7.5 ) );
	Book( TH1F( "flavor_topjet_after_selection4", "decay of topjet after selection4", 8, -0.5 ,7.5 ) );
  Book( TH1F( "pT"," p_{T} topjets",100,0,2000));
    Book( TH1F( "HT"," H_{T} ",100,0,2000));
      Book( TH1F( "HT2"," H_{T} 2",100,0,2000));
        Book( TH1F( "HTlep"," p_{Tlep} ",100,0,2000));
	  Book( TH1F( "HTlep2"," p_{Tlep} 2",100,0,2000));
   Book( TH1F( "pT_tag"," p_{T} topjets decay products in jet",100,0,2000));
    Book( TH1F( "pT_mistag"," p_{T} topjets decay products not in jet",100,0,2000));
 Book( TH1F( "pT_toplep"," p_{T} topjets",100,0,2000));
  Book( TH1F( "pT_topleptag"," p_{T} topjets",100,0,2000));
  Book( TH1F( "pT_frac"," p_{T} topjets",100,0,2000));
   Book( TH1F( "pT_frac2"," p_{T} topjets",100,0,2000));
  Book( TH1F( "chi2","chi2",100,0,500));
  Book( TH1F( "pT_after_selection", "p_{T} topjets after selection",100,0,2000));
  Book( TH1F( "pT_theory", "p_{T} topjets theory",100,0,2000));
   Book( TH1F( "pT_theory_tagged", "p_{T} topjets theory tagged",100,0,2000));
  Book( TH1F( "pT_tagged_cms", "p_{T} topjets (tagged)",100,0,2000));
  Book( TH1F( "pT_tagged_hep", "p_{T} topjets (tagged)",100,0,2000));
  Book( TH1F( "pT_tagged_cms_after_selection", "p_{T} topjets (tagged) after selection", 100, 0, 2000));
  Book( TH1F( "pT_tagged_hep_after_selection", "p_{T} topjets (tagged) after selection", 100, 0, 2000));
Book( TH1F( "pT_substract", "p_{T} topjets (tagged) after selection not in jet", 50, 0, 2000));
 Book( TH1F( "pT_substract_cms_tagged", "p_{T} topjets (cms tagged) after selection not in jet tagged", 50, 0, 2000));
 Book( TH1F( "pT_substract_mva_tagged", "p_{T} topjets (mva tagged) after selection not in jet tagged", 50, 0, 2000));
  Book( TH1F( "eff_pT_cms","efficiency in p_{T}(CMS-tagger)",100,0,2000));
  Book( TH1F( "eff_pT_hep","efficiency in p_{T}(HEP-tagger)",100,0,2000));
  Book( TH1F( "eff_pT_cms_theory","efficiency in p_{T} theory  (CMS-tagger)",100,0,2000));
  Book( TH1F( "eff_pT_hep_theory","efficiency in p_{T} theory  (HEP-tagger)",100,0,2000));
  Book( TH1F( "eff_pT_cms_after_selection", "efficiency in p_{T} after selection",100,0,2000));
  Book( TH1F( "eff_pT_hep_after_selection", "efficiency in p_{T} after selection",100,0,2000));
  Book( TH1F ( "lept_ID", "lept ID",3,-0.5,2.5) );
  Book( TH1F( "distance", "distance topjet topjet_lept", 1000,0,10));
  Book( TH1F( "distance2", "distance topjet topjet_lept", 1000,0,10));
  Book( TH1F( "distance3", "distance topjet topjet_lept", 1000,0,10));
   Book( TH1F( "distance4", "distance topjet topjet_lept", 1000,0,10));
   Book( TH1F( "distanceq1", "distance topjet quark1", 1000,0,10));
  Book( TH1F( "distanceq2", "distance topjet quark2", 1000,0,10));
  Book( TH1F( "distanceq3", "distance topjet quark3", 1000,0,10));
   Book( TH1F( "distance_jet1", "distance topjet jet (decay products)", 1000,0,10));
    Book( TH1F( "distance_jet2", "distance topjet jet (!decay products)", 1000,0,10));
  Book( TH1F( "pT_ly"," p_{T} topjets",100,0,2000));
   Book( TH1F( "pT_diff"," p_{T} topjet-top_lep",400,-2000,2000));
    Book( TH1F( "pT_diff2"," p_{T} topjet-top_lep",400,-2000,2000));
  Book( TH1F( "eta"," #eta topjets",100,-3,3));
  Book( TH1F( "eta_ly"," #eta topjets",100,-3,3));
  Book( TH1F( "phi"," #phi topjets",100,-M_PI,M_PI));
  Book( TH1F( "phi_ly"," #phi topjets",100,-M_PI,M_PI));
  Book( TH1F( "MJet", "m_{jet}", 100,0,400 ) );
  Book( TH1F( "MJet_ly", "m_{jet}", 100,0,400 ) );
  Book( TH1F( "Mmin", "m_{min}", 100,0,160 ) );
  Book( TH1F( "Mmin_ly", "m_{min}", 100,0,160 ) );
  Book( TH1F( "NSubjets", "number of subjets", 6,-0.5,5.5) );
   Book( TH1F( "NSubjets_pt0to300", "number of subjets pt 0to300", 6,-0.5,5.5) );
Book( TH1F( "NSubjets_pt300to500", "number of subjets pt 0to300", 6,-0.5,5.5) );
Book( TH1F( "NSubjets_pt500to700", "number of subjets pt 0to300", 6,-0.5,5.5) );
Book( TH1F( "NSubjets_pt700to900", "number of subjets pt 0to300", 6,-0.5,5.5) );
Book( TH1F( "NSubjets_pt900toinf", "number of subjets pt 0to300", 6,-0.5,5.5) );
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

  // Book(TGraphAsymmErrors(Hist("pT"),Hist("pT_tagged"),""));
  //tmva_tagger=TMVA_tagger::Instance();
  tmva_tagger=new TMVA_tagger();
  tmva_tagger->Set_Reader("bestown70");
}


void TopTagHists::Fill()
{
   // important: get the event weight
   EventCalc* calc = EventCalc::Instance();
   TopFitCalc* topfit = TopFitCalc::Instance();
  double weight = calc -> GetWeight();

  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();
  //--------------------------------
   bcc->recoHyps->clear();
  topfit->CalculateSelection(); 
  tagchi2discr = new Chi2Discriminator();
tagchi2discr->FillDiscriminatorValues();
 ReconstructionHypothesis *discr = tagchi2discr->GetBestHypothesis();
 ReconstructionHypothesis *hyp= tagchi2discr->GetBestHypothesis();;
   

  //Hist("NRecoHyps")->Fill(bcc->recoHyps->size(),weight);
  //std::cout<<bcc->recoHyps->size()<<std::endl;
  // for(unsigned int i=0;i<bcc->recoHyps->size();i++){
  //hyp = bcc->recoHyps->at(i);
 
    LorentzVector top_lep2 = hyp->toplep_v4();
    LorentzVector top_had = hyp->tophad_v4();
    LorentzVector top_lep = discr->toplep_v4();
    //  double discr_cut = discr->discriminator("Chi2_tlep");
    // Hist("chi2")->Fill(discr_cut,weight);
    // if(discr_cut<10){
      //if(vec_matches_leptonic_decay(top_lep,0.2)) //Hist("flavor_topjet")->Fill("leptonically2",weight);
      //if(vec_matches_leptonic_decay(top_lep2,0.2)) //Hist("flavor_topjet")->Fill("leptonically3",weight);
    // }
    //for(unsigned i = 0;i<calc->GetGenParticles()->size();++i){
    //GenParticle top=calc->GetGenParticles()->at(i);
    //if(abs(calc->GetGenParticles()->at(i).pdgId())==6) 
    
       Hist("distance3")->Fill(sqrt(pow(top_lep.phi()-top_had.phi(),2)+pow(top_lep.eta()-top_had.eta(),2)), weight);
     //}

    /*if(bcc->recoHyps->size()>0)  hyp = bcc->recoHyps->at(0);
      LorentzVector tobi = hyp.toplep_v4();*/
  



  int NTopJets = bcc-> topjets -> size();
  int NMuon = bcc->muons->size();
  int NElectrons = bcc->electrons->size();
  int NTaus = bcc->taus->size();
  sort(bcc->muons->begin(), bcc->muons->end(), HigherPt());
  Hist("NMuons")->Fill(NMuon,weight);
  Hist("NElectrons")->Fill(NElectrons,weight);
  Hist("NTaus")->Fill(NTaus,weight);
  Hist("NTopJets")->Fill(NTopJets, weight);
  Hist("NTopJets_ly")->Fill(NTopJets, weight);
  NTopJets=0;
  bool trigger = false;
  bool selection_thad;
  sort(bcc->topjets->begin(), bcc->topjets->end(), HigherPt());
  Muon muon = bcc->muons->at(0);
  

  //recoonstruct lep jet
  /* TopFitCalc *bla;
  Particle* lepton = &bcc->muons->at(0);
  std::vector<LorentzVector> neutrinos = bla->NeutrinoFitPolar(lepton->v4(),bcc->met->v4());*/


  // TopJet topjet_lept;
  // //check if there is leptonically top in event
  // for (unsigned int i = 0; i<bcc->topjets->size(); ++i)
  //   {
  //     TopJet topjet =  bcc->topjets->at(i);
      
      
  //     //mass of leptonic jet?
  //          LorentzVector allsubjets(0,0,0,0);
  //     for(int k=0; k<topjet.numberOfDaughters(); ++k) {
  //       allsubjets += topjet.subjets()[k].v4();

  // 	}
  //     double mjet2=0;
  //     double mmin2=0;
  //     int nsubjets2=0;
  //     allsubjets +=muon.v4();
  //     double mjet = allsubjets.M();
  //     if(WTag(topjet,mjet2,nsubjets2,mmin2)){
  // 	//  if(jet_decays_leptonic(topjet,0.8)) Hist("distance2")->Fill(mjet,weight);
  // 	//if(!jet_matches_top(topjet,0.8)) Hist("distance")->Fill(mjet,weight);
  //     if(sqrt(pow(topjet.phi()-muon.phi(),2)+pow(topjet.eta()-muon.eta(),2))<2.0 && mjet<150) {
  // 	trigger = true;
  // 	topjet_lept=topjet;
  // 	if(jet_decays_leptonic(topjet,0.8)) Hist("lept_ID")->Fill(1,weight);
  // 	else Hist("lept_ID")->Fill(2,weight);
 
  //     }
  //     else if(jet_decays_leptonic(topjet,0.8)) Hist("lept_ID")->Fill(0);
  //     }
  //   }
  bool select_event=false;
 


 for(unsigned int i =0; i<bcc->topjets->size();++i){
    TopJet topjet = bcc->topjets->at(i); 
    double mmin=0;
      double mjet=0;
      int nsubjets=0;
      double cone=0.8;
      bool products_in_jet2=(distance_quark(topjet,1,0.8)<0.8 && distance_quark(topjet,2,0.8)<0.8 &&  distance_quark(topjet,3,0.8)<0.8);
      Hist("pT") -> Fill(topjet.pt(),weight);
      if(TopTag(topjet,mjet,nsubjets,mmin)) Hist("pT_tagged_cms") -> Fill(topjet.pt(),weight);
      if(topjet.pt()<300) Hist("NSubjets_pt0to300")->Fill(nsubjets,weight);
      if(topjet.pt()>300 && topjet.pt()<500) Hist("NSubjets_pt300to500")->Fill(nsubjets,weight);
       if(topjet.pt()>500 && topjet.pt()<700) Hist("NSubjets_pt500to700")->Fill(nsubjets,weight);
       if(topjet.pt()>700 && topjet.pt()<900) Hist("NSubjets_pt700to900")->Fill(nsubjets,weight);
       if(topjet.pt()>900 ) Hist("NSubjets_pt900toinf")->Fill(nsubjets,weight);
       Hist("nusbjetsvspt")->Fill(topjet.pt(),nsubjets);
      if(HepTopTag(topjet)) Hist("pT_tagged_hep") -> Fill(topjet.pt(),weight);
      if(products_in_jet2) Hist("pT_theory") ->Fill(topjet.pt(),weight);
      if(products_in_jet2 && TopTag(topjet,mjet,nsubjets,mmin)) Hist("pT_theory_tagged") ->Fill(topjet.pt(),weight);
      // if(decay_products_in_jet_had(topjet,100) && jet_decays_hadronic(topjet,0.8)) Hist("flavor_topjet")->Fill("proinjet_pre",weight);
      //else if(jet_decays_hadronic(topjet,0.8)) Hist("flavor_topjet")->Fill("pronotjet_pre",weight);
	//decay products in jet?
      bool products_in_jet=(distance_quark(topjet,1,0.8)<0.8 && distance_quark(topjet,2,0.8)<0.8 &&  distance_quark(topjet,3,0.8)<0.8);
      if(jet_decays_hadronic(topjet,0.8)){
	if(topjet.pt()>0 && topjet.pt()<400){
	if(products_in_jet) Hist("decay_products0")->Fill("in_jet",weight);
	else Hist("decay_products0")->Fill("not_in_jet",weight);}
	if(topjet.pt()>400  && topjet.pt()<700){
	  if(products_in_jet) Hist("decay_products200")->Fill("in_jet",weight);
	  else Hist("decay_products200")->Fill("not_in_jet",weight);}
	if(topjet.pt()>700 ){
	if(products_in_jet) Hist("decay_products400")->Fill("in_jet",weight);
	else Hist("decay_products400")->Fill("not_in_jet",weight);}}
	/*
	if(topjet.pt()>600 && topjet.pt()<800 )
	if(decay_products_in_jet_had(topjet,cone)) Hist("decay_products600")->Fill("in_jet",weight);
	else Hist("decay_products600")->Fill("not_in_jet",weight);
	if(topjet.pt()>800 && topjet.pt()<1000 )
	if(decay_products_in_jet_had(topjet,cone)) Hist("decay_products800")->Fill("in_jet",weight);
	else Hist("decay_products800")->Fill("not_in_jet",weight);
	if(topjet.pt()>1000)
	if(decay_products_in_jet_had(topjet,cone)) Hist("decay_products1000")->Fill("in_jet",weight);
	else Hist("decay_products1000")->Fill("not_in_jet",weight);
	}*/
	if(jet_decays_hadronic(topjet,cone) && bcc->topjets->size()==2){
	  if(topjet.pt()>0 && topjet.pt()<400 ){
	if(products_in_jet) Hist("decay_products0")->Fill("in_jet2",weight);
	else Hist("decay_products0")->Fill("not_in_jet2",weight);}
	  if(topjet.pt()>400 && topjet.pt()<700 ){
	if(products_in_jet) Hist("decay_products200")->Fill("in_jet2",weight);
	else Hist("decay_products200")->Fill("not_in_jet2",weight);}
	  if(topjet.pt()>700 ){
	  if(products_in_jet) Hist("decay_products400")->Fill("in_jet2",weight);
	  else Hist("decay_products400")->Fill("not_in_jet2",weight);}
	}/*
	if(topjet.pt()>600 && topjet.pt()<800 )
	if(decay_products_in_jet_had(topjet,cone)) Hist("decay_products600")->Fill("in_jet2",weight);
	else Hist("decay_products600")->Fill("not_in_jet2",weight);
	if(topjet.pt()>800 && topjet.pt()<1000 )
	if(decay_products_in_jet_had(topjet,cone)) Hist("decay_products800")->Fill("in_jet2",weight);
	else Hist("decay_products800")->Fill("not_in_jet2",weight);
	if(topjet.pt()>1000 )
	if(decay_products_in_jet_had(topjet,cone)) Hist("decay_products1000")->Fill("in_jet2",weight);
	else Hist("decay_products1000")->Fill("not_in_jet2",weight);
	}
	if(jet_decays_hadronic(topjet,cone)){
	if(decay_products_in_jet_had(topjet,cone)) Hist("pT_tag")->Fill(topjet.pt(),weight);
	else Hist("pT_mistag")->Fill(topjet.pt(),weight);}*/
	if(jet_decays_hadronic(topjet,0.8)){
	Hist("distanceq1")->Fill(distance_quark(topjet,1,0.8));
	Hist("distanceq2")->Fill(distance_quark(topjet,2,0.8));
	Hist("distanceq3")->Fill(distance_quark(topjet,3,0.8)); }
}


 int nbtags=0;
for (unsigned int i =0; i<bcc->jets->size(); ++i)
    {
      Jet jet =  bcc->jets->at(i);
      if(IsTagged(jet,e_CSVT)) nbtags++;
    }
 Hist("Nbtags")->Fill(nbtags,weight);
 if(nbtags>0) select_event=true;
//if(discr_cut>25) select_event=false;
 if(bcc->topjets->size()!=2) ;
  // selection_thad=true;
 if(select_event)
  for (unsigned int i =0; i<bcc->topjets->size(); ++i)
    {
      TopJet topjet =  bcc->topjets->at(i);
       

     





      // Muon muon1 = bcc->muons->at(1);
        selection_thad = true;
      Hist("distance2")->Fill(sqrt(pow(top_had.phi()-topjet.phi(),2)+pow(top_had.eta()-topjet.eta(),2)), weight);
      // Hist("distance")->Fill(sqrt(pow(top_lep.phi()-topjet.phi(),2)+pow(top_lep.eta()-topjet.eta(),2)), weight);
      //  if(sqrt(pow(topjet.phi()-muon.phi(),2)+pow(topjet.eta()-muon.eta(),2))<2.0 ) {
      //	selection_thad = false;
      //	}
       if(sqrt(pow(topjet.phi()-top_lep.phi(),2)+pow(topjet.eta()-top_lep.eta(),2))<2.7 ||  sqrt(pow(topjet.phi()-top_lep.phi(),2)+pow(topjet.eta()-top_lep.eta(),2))>3.5) selection_thad = false;
      //if(sqrt(pow(top_had.phi()-top_lep.phi(),2)+pow(top_had.eta()-top_lep.eta(),2))<2.5 ||  sqrt(pow(top_had.phi()-top_lep.phi(),2)+pow(top_had.eta()-top_lep.eta(),2))>4.5) selection_thad = false;
       //if(top_lep.pt()/topjet.pt()>1) selection_thad = false;
       //if(topjet.pt()<350) selection_thad = false;
       if(topjet.pt()!=0) Hist("pT_toplep") -> Fill(top_lep.pt(),weight);
       if(topjet.pt()!=0 && jet_decays_hadronic(topjet,0.8)) Hist("pT_topleptag") -> Fill(top_lep.pt(),weight);
     
      double cone = 0.8;
       double mmin=0;
      double mjet=0;
      int nsubjets=0;
      if(TopTag(topjet,mjet,nsubjets,mmin)) Hist("distance")->Fill(sqrt(pow(top_lep.phi()-topjet.phi(),2)+pow(top_lep.eta()-topjet.eta(),2)), weight);
      if(jet_decays_hadronic(topjet,0.8)) Hist("distance4")->Fill(sqrt(pow(top_lep.phi()-topjet.phi(),2)+pow(top_lep.eta()-topjet.eta(),2)), weight);
      
       bool products_in_jet=(distance_quark(topjet,1,0.8)<0.8 && distance_quark(topjet,2,0.8)<0.8 &&  distance_quark(topjet,3,0.8)<0.8);
      //check selection
       //if(topjet.pt()>400){
       if(products_in_jet && jet_decays_hadronic(topjet,0.8))  Hist("flavor_topjet")->Fill("products_in_jet",weight);
       else {
	 Hist("flavor_topjet")->Fill("rest",weight);
	 if(jet_decays_hadronic(topjet,0.8) && !products_in_jet) Hist("flavor_topjet")->Fill("not_in_jet(hadronically)",weight);
	 if(jet_decays_leptonic(topjet,0.8))   Hist("flavor_topjet")->Fill("leptonically",weight);
	 if(!jet_decays_hadronic(topjet,0.8) && !jet_decays_leptonic(topjet,0.8)) Hist("flavor_topjet")->Fill("gluon?",weight);}//}
	// else Hist("flavor_topjet")->Fill("not matching",weight);
      
     
      
      if(selection_thad ) {
	

	if(jet_decays_hadronic(topjet,cone) && bcc->topjets->size()==2){
	  if(topjet.pt()>0 && topjet.pt()<400 ){
	if(products_in_jet) Hist("decay_products0")->Fill("in_jet3",weight);
	else Hist("decay_products0")->Fill("not_in_jet3",weight);}
	  if(topjet.pt()>400 && topjet.pt()<700 ){
	if(products_in_jet) Hist("decay_products200")->Fill("in_jet3",weight);
	else Hist("decay_products200")->Fill("not_in_jet3",weight);}
	  if(topjet.pt()>700 ){
	  if(products_in_jet) Hist("decay_products400")->Fill("in_jet3",weight);
	  else Hist("decay_products400")->Fill("not_in_jet3",weight);}
	}
	bool jet_distance=true;
	double jetclosetop=0;
	double jetnotclosetop=0;
	double jetclosetop2=0;
	double jetnotclosetop2=0;
	LorentzVector alljet(0,0,0,0);
	LorentzVector alljet2(0,0,0,0);
	alljet=topjet.v4()+top_lep;
	if(!products_in_jet) alljet2=topjet.v4()+top_lep;
		for(unsigned t=0;t<bcc->jets->size();++t){
		   Jet jet = bcc->jets->at(t);
	if(products_in_jet) Hist("distance_jet1")->Fill(sqrt(pow(jet.phi()-top_lep.phi(),2)+pow(jet.eta()-top_lep.eta(),2)));
	if(!products_in_jet) Hist("distance_jet2")->Fill(sqrt(pow(jet.phi()-top_lep.phi(),2)+pow(jet.eta()-top_lep.eta(),2)));
	//if(sqrt(pow(topjet.phi()-jet.phi(),2)+pow(topjet.eta()-jet.eta(),2))>0.09 && sqrt(pow(topjet.phi()-jet.phi(),2)+pow(topjet.eta()-jet.eta(),2))<0.8) jet_distance=true;
	if((sqrt(pow(topjet.phi()-jet.phi(),2)+pow(topjet.eta()-jet.eta(),2))>0.8 &&sqrt(pow(topjet.phi()-jet.phi(),2)+pow(topjet.eta()-jet.eta(),2))<1.8) ) jet_distance=false;
	//if((sqrt(pow(top_lep.phi()-jet.phi(),2)+pow(top_lep.eta()-jet.eta(),2))>1 &&sqrt(pow(top_lep.phi()-jet.phi(),2)+pow(top_lep.eta()-jet.eta(),2))<2.2)) jet_distance=false;
	//	if((sqrt(pow(top_lep.phi()-jet.phi(),2)+pow(top_lep.eta()-jet.eta(),2))>4)) jet_distance=false;
	if(sqrt(pow(topjet.phi()-jet.phi(),2)+pow(topjet.eta()-jet.eta(),2))<0.8) jetclosetop++; 
	else if(sqrt(pow(topjet.phi()-jet.phi(),2)+pow(topjet.eta()-jet.eta(),2))<1.8) {
	  jetnotclosetop++;
	  }
	
	 if(!products_in_jet) Hist("jetclosetop")->Fill(jetclosetop,weight);
	 if(!products_in_jet)Hist("jetnotclosetop")->Fill(jetnotclosetop,weight);
	 if(products_in_jet) Hist("jetclosetop2")->Fill(jetclosetop,weight);
	 if(products_in_jet) Hist("jetnotclosetop2")->Fill(jetnotclosetop,weight);
	

	 
	if(IsTagged(jet,e_CSVT)){
	if(sqrt(pow(topjet.phi()-jet.phi(),2)+pow(topjet.eta()-jet.eta(),2))<1.3) jetclosetop2++; 
	else if(sqrt(pow(topjet.phi()-jet.phi(),2)+pow(topjet.eta()-jet.eta(),2))<2.5) jetnotclosetop2++;
	 Hist("jetclosetop2")->Fill(jetclosetop2,weight);
	 Hist("jetnotclosetop2")->Fill(jetnotclosetop2,weight);
	}}
	 if(topjet.pt()>400){
	 Hist("MJet_out")->Fill(alljet.M(),weight);
	 Hist("MJet_out2")->Fill(alljet2.M(),weight);}
	 int npfparticle=0;
	 for(unsigned t=0;t<bcc->pfparticles->size(); ++t){
	   PFParticle pfparticle = bcc->pfparticles->at(t);
	   if(sqrt(pow(topjet.phi()-pfparticle.phi(),2)+pow(topjet.eta()-pfparticle.eta(),2))>0.5 && sqrt(pow(topjet.phi()-pfparticle.phi(),2)+pow(topjet.eta()-pfparticle.eta(),2))<2.5 &&topjet.pt()>400)	 npfparticle++;
	   /*  Hist("NPFflow")->Fill(npfparticle,weight);
	       if(products_in_jet) Hist("NPFflow2")->Fill(npfparticle,weight);*/
	 }
	 std::vector<PFParticle> pfparticle2=calc->GetJetPFParticles(&topjet);
	 if(!products_in_jet && topjet.pt()>400) Hist("NPFflow")->Fill(pfparticle2.size(),weight);
	 if(products_in_jet && topjet.pt()>400) Hist("NPFflow2")->Fill(pfparticle2.size(),weight);

	 //if(jetnotclosetop2!=0 && jetclosetop2==1){
	 /*	if(jet_matches_top(topjet,cone)) {
	  if(jet_decays_hadronic(topjet,cone)) {if(products_in_jet) Hist("flavor_topjet_after_selection")->Fill("hadronically_in_jet",weight);
	    else Hist("flavor_topjet_after_selection")->Fill("hadronically",weight);}
	if(jet_decays_leptonic(topjet,cone)) Hist("flavor_topjet_after_selection")->Fill("leptonically",weight);
	//}
	
	else if(products_in_jet) Hist("flavor_topjet_after_selection")->Fill("not matching",weight);
	else Hist("flavor_topjet_after_selection")->Fill("products not in jet",weight);}*/
	 // if(jetnotclosetop2==0 ){

	 double HTlep=calc->GetHTlep();
	 double HT=calc->GetHT();
	 if(!products_in_jet) Hist("HT")->Fill(HT,weight);
	 if(products_in_jet) Hist("HT2")->Fill(HT,weight);
	 if(!products_in_jet) Hist("HTlep")->Fill(HTlep,weight);
	 if(products_in_jet) Hist("HTlep2")->Fill(HTlep,weight);
	 
	 /* if(products_in_jet && jet_decays_hadronic(topjet,0.8)) Hist("flavor_topjet_after_selection")->Fill("products_in_jet",weight);
	 else Hist("flavor_topjet_after_selection")->Fill("products_not_in_jet",weight);
	 if(jet_decays_leptonic(topjet,0.8))  Hist("flavor_topjet_after_selection")->Fill("leptonically",weight);*/

	  if(products_in_jet && jet_decays_hadronic(topjet,0.8))  Hist("flavor_topjet_after_selection")->Fill("products_in_jet",weight);
	  else {
	 Hist("flavor_topjet_after_selection")->Fill("rest",weight);
	 if(jet_decays_hadronic(topjet,0.8) && !products_in_jet) Hist("flavor_topjet_after_selection")->Fill("not_in_jet(hadronically)",weight);
	 if(jet_decays_leptonic(topjet,0.8))   Hist("flavor_topjet_after_selection")->Fill("leptonically",weight);
	 if(!jet_decays_hadronic(topjet,0.8) && !jet_decays_leptonic(topjet,0.8)) Hist("flavor_topjet_after_selection")->Fill("gluon?",weight);}
 if(products_in_jet) Hist("pT_diff")->Fill(topjet.pt()-top_lep.pt(),weight);
	 if(!products_in_jet) Hist("pT_diff2")->Fill(topjet.pt()-top_lep.pt(),weight);

	 if(/*topjet.pt()>400 &&*//*jetnotclosetop==0*/ jet_distance /*&& jetclosetop2==1*/ /*&& alljet.M()<50 && pfparticle2.size()>80*//*&& alljet.M()>950*/){
	   /* if(products_in_jet && jet_decays_hadronic(topjet,0.8)) Hist("flavor_topjet_after_selection2")->Fill("products_in_jet",weight);
	 else Hist("flavor_topjet_after_selection2")->Fill("products_not_in_jet",weight);
	 if(jet_decays_leptonic(topjet,0.8))  Hist("flavor_topjet_after_selection2")->Fill("leptonically",weight);}//}*/
	   
	  if(products_in_jet && jet_decays_hadronic(topjet,0.8))  Hist("flavor_topjet_after_selection2")->Fill("products_in_jet",weight);
	  else {
	 Hist("flavor_topjet_after_selection2")->Fill("rest",weight);
	 if(jet_decays_hadronic(topjet,0.8) && !products_in_jet) Hist("flavor_topjet_after_selection2")->Fill("not_in_jet(hadronically)",weight);
	 if(jet_decays_leptonic(topjet,0.8))   Hist("flavor_topjet_after_selection2")->Fill("leptonically",weight);
	 if(!jet_decays_hadronic(topjet,0.8) && !jet_decays_leptonic(topjet,0.8)) Hist("flavor_topjet_after_selection2")->Fill("gluon?",weight);
	  }


	 if(products_in_jet && jet_decays_hadronic(topjet,0.8))  Hist("flavor_topjet_after_selection3")->Fill("products_in_jet",weight);
	  else {
	 Hist("flavor_topjet_after_selection3")->Fill("rest",weight);
	 if(jet_decays_hadronic(topjet,0.8) && !products_in_jet) Hist("flavor_topjet_after_selection3")->Fill("not_in_jet(hadronically)",weight);
	 if(jet_decays_leptonic(topjet,0.8))   Hist("flavor_topjet_after_selection3")->Fill("leptonically",weight);
	 if(!jet_decays_hadronic(topjet,0.8) && !jet_decays_leptonic(topjet,0.8)) Hist("flavor_topjet_after_selection3")->Fill("gluon?",weight);

	  }
    double mmin=0;
      double mjet=0;
      int nsubjets=0;
      //   tmva_tagger->push_variables(topjet);
      double mva_value;//=tmva_tagger->GetMVA_value();
      	 if(!products_in_jet) Hist("pT_substract")->Fill(topjet.pt(),weight); 
	 if(!products_in_jet && TopTag(topjet,mjet,nsubjets,mmin)) Hist("pT_substract_cms_tagged")->Fill(topjet.pt(),weight);
	 if(!products_in_jet && mva_value>0.165) Hist("pT_substract_mva_tagged")->Fill(topjet.pt(),weight);


	 }


	
	 if(topjet.pt()>400 && /*&&  jetnotclosetop==0 &&*/ jet_distance /*&& alljet.M()<50 && pfparticle2.size()>80*//*&& alljet.M()>950*/){
	   /*if(products_in_jet && jet_decays_hadronic(topjet,0.8)) Hist("flavor_topjet_after_selection3")->Fill("products_in_jet",weight);
	 else Hist("flavor_topjet_after_selection3")->Fill("products_not_in_jet",weight);
	 if(jet_decays_leptonic(topjet,0.8))  Hist("flavor_topjet_after_selection3")->Fill("leptonically",weight);*/
	 
	  if(products_in_jet && jet_decays_hadronic(topjet,0.8))  Hist("flavor_topjet_after_selection4")->Fill("products_in_jet",weight);
	  else {
	 Hist("flavor_topjet_after_selection4")->Fill("rest",weight);
	 if(jet_decays_hadronic(topjet,0.8) && !products_in_jet) Hist("flavor_topjet_after_selection4")->Fill("not_in_jet(hadronically)",weight);
	 if(jet_decays_leptonic(topjet,0.8))   Hist("flavor_topjet_after_selection4")->Fill("leptonically",weight);
	 if(!jet_decays_hadronic(topjet,0.8) && !jet_decays_leptonic(topjet,0.8)) Hist("flavor_topjet_after_selection4")->Fill("gluon?",weight);}
	 
	 
	 }
     
	  if(jet_distance ){
	   Hist("pT_after_selection")->Fill(topjet.pt(),weight);
      double mmin=0;
      double mjet=0;
      int nsubjets=0;
      if(TopTag(topjet,mjet,nsubjets,mmin)) Hist("pT_tagged_cms_after_selection") -> Fill(topjet.pt(),weight);
      if(HepTopTag(topjet)) Hist("pT_tagged_hep_after_selection") -> Fill(topjet.pt(),weight);
   
      NTopJets++;
      }
      }
      Hist("pT_ly") -> Fill(topjet.pt(), weight);
     
      Hist("eta") -> Fill(topjet.eta(), weight);
      Hist("eta_ly") -> Fill(topjet.eta(), weight);
      Hist("phi") -> Fill(topjet.phi(), weight);
      Hist("phi_ly") -> Fill(topjet.phi(), weight);
      
    
     
      Hist( "MJet" )->Fill( mjet, weight );
      Hist( "MJet_ly" )->Fill( mjet, weight );
      
      if(nsubjets>=3) 
	{
	  Hist( "Mmin" )->Fill( mmin, weight );
	  Hist( "Mmin_ly" )->Fill( mmin, weight );
	  
	}
      Hist( "NSubjets" )->Fill( nsubjets, weight ); 
      Hist( "NSubjets_ly" )->Fill( nsubjets, weight ); 
      
     
    }
  Hist("NTopJets_after_selection")->Fill(NTopJets,weight);
  
  

 
 
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

void TopTagHists::Finish()
{
  Hist("eff_pT_cms") -> Divide(Hist("pT_tagged_cms"),Hist("pT"));
  Hist("eff_pT_hep") -> Divide(Hist("pT_tagged_hep"),Hist("pT"));
   Hist("eff_pT_cms_theory") -> Divide(Hist("pT_theory_tagged"),Hist("pT_theory"));
  Hist("eff_pT_hep_theory") -> Divide(Hist("pT_tagged_hep"),Hist("pT_theory"));
  Hist("eff_pT_cms_after_selection") -> Divide(Hist("pT_tagged_cms_after_selection"),Hist("pT_after_selection"));
  Hist("eff_pT_hep_after_selection") -> Divide(Hist("pT_tagged_hep_after_selection"),Hist("pT_after_selection"));
  Hist("pT_frac") -> Divide(Hist("pT_toplep"),Hist("pT"));
  Hist("pT_frac2") -> Divide(Hist("pT_topleptag"),Hist("pT"));
   // final calculations, like division and addition of certain histograms
  
}


