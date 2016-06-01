#include "include/TopTagger.h"
//#include "NtupleWriter/include/JetProps.h"
//#include "NtupleWriter/interface/GenJetProps.h"
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include "SFrameTools/include/EventCalc.h"
#include "fastjet/contrib/HOTVR.hh"
#include "fastjet/contrib/SoftDrop.hh"

using namespace std;
using namespace fastjet;
using namespace contrib;

TopTagger* TopTagger::m_instance = NULL;

TopTagger* TopTagger::Instance()
{
  // Get a pointer to the object handler.
  // This is the only way to access this class, 
  // since it's a singleton. This method is accessible
  // from everywhere.

  if (m_instance == NULL){
    m_instance = new  TopTagger();
  }
  return m_instance; 
   
}


TopTagger::TopTagger() : m_logger("TopTagger")
{
  // constructor: initialise all variables

  
 
     //Reset();

  
}



TopTagger::~TopTagger()
{
  // default destructor
  
}






fastjet::PseudoJet TopTagger::get_softdrop_jet(fastjet::PseudoJet jet, double zcut, double beta){
  fastjet::PseudoJet sd_jet;
  fastjet::contrib::SoftDrop sd(beta, zcut);
  std::cout<<"1"<<std::endl;
  sd_jet = sd(jet);
  std::cout<<"2"<<std::endl;
  return sd_jet;
  

}



bool TopTagger::Is_tagged(TString tagger, fastjet::PseudoJet jet){
  if(tagger== "HOTVR" || tagger=="hotvr") return Is_HOTVR_tagged(jet);
  if(tagger== "CMS" || tagger=="cms")   return Is_CMS_tagged(jet);
  if(tagger== "HEP" || tagger=="hep")    return Is_HEP_tagged(jet);
  if(tagger== "optimalR" || tagger=="optimalr" || tagger=="OPTIMALR" || tagger=="OptimalR")   return Is_OptimalR_tagged(jet);
  if(tagger=="softdrop" || tagger=="SoftDrop" || tagger=="SOFTDROP" || tagger=="Softdrop") return Is_softdrop_tagged(jet);
  else  return false;
  
}





bool TopTagger::Is_softdrop_tagged(fastjet::PseudoJet &jet){
  std::vector<fastjet::PseudoJet> subjets;
  fastjet::PseudoJet sd_jet = get_softdrop_jet(jet,0.2,1);
  jet.set_user_info(new HOTVRinfo(jet,sorted_by_pt(subjets)));
  if(sd_jet.m()>140 && sd_jet.m()<250) return true;
  else return false;
}


bool TopTagger::Is_CMS_tagged(fastjet::PseudoJet &jet){
  double mjet,  mmin;
  int nsubjets;
  std::vector<fastjet::PseudoJet> subjets;
  fastjet::PseudoJet CMSjet;
  bool tagged=CMSTopTagFull_pseudo_CA(jet,3.0,0.05,0.4,0.0004,mjet, nsubjets,mmin,subjets,CMSjet);
  jet.set_user_info(new HOTVRinfo(jet,sorted_by_pt(subjets)));
  //jet=CMSjet;
  return tagged;

 }


bool TopTagger::Is_HOTVR_tagged(fastjet::PseudoJet &jet){
  double mjet=jet.m();
  int nsubjets=jet.user_info<HOTVRinfo>().subjets().size();
  std::vector<fastjet::PseudoJet> subjets= jet.user_info<HOTVRinfo>().subjets();
  double ptfraction=jet.user_info<HOTVRinfo>().ptfraction(0);
  double mmin=jet.user_info<HOTVRinfo>().mmin();
  bool tagged=false;
  if(mjet>140 && mjet<220 && nsubjets>2 && mmin>50 && ptfraction<0.8 /*&& subjets.at(0).pt()>30 && subjets.at(1).pt()>30*/) tagged=true;
  return tagged;
}


bool TopTagger::Is_HEP_tagged(fastjet::PseudoJet &jet){
  double m12,m13,m23,fw,m123,m_pruned,m_unfiltered;
  std::vector<fastjet::PseudoJet> subjets;
  bool tagged=  HepTopTagFull_pseudo(jet, subjets, m12, m13, m23, fw, m123, m_pruned, m_unfiltered);
  jet.set_user_info(new HOTVRinfo(jet,sorted_by_pt(subjets)));
  return tagged;
}


bool TopTagger::Is_OptimalR_tagged(fastjet::PseudoJet &jet){
  double fw,Rmin,mass_Rmin,pt_Rmin,mass_diff,pt_for_exp,Rmin_exp,fw_out;
  std::vector<fastjet::PseudoJet> subjets;
  bool tagged=false;
  if(MultiRTopTag_pseudo(jet, subjets,0.15, Rmin,mass_Rmin,pt_Rmin,mass_diff,pt_for_exp,Rmin_exp,fw_out) && Rmin-Rmin_exp<0.5  && mass_Rmin>140 && mass_Rmin<250) tagged=true; 
// std::cout<<Rmin-Rmin_exp<<" "<<mass_Rmin<<std::endl;
// if(tagged) std::cout<<"this is a tag"<<std::endl;
 return tagged;
}
