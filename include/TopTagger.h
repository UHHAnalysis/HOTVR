#ifndef TopTagger_H
#define TopTagger_H

#include "SFrameTools/include/Objects.h"
#include <fastjet/PseudoJet.hh>
#include "SFrameTools/include/fwd.h"
#include "SFrameTools/include/boost_includes.h" // for shared_array
#include <TMVA/Reader.h>
#include "TVector3.h"
#include <limits>
#include <algorithm>
#include <memory>
#include <TF1.h>
#include "Utils.h"
#include "EventCalc.h"
#include "include/Cleaner.h"
#include "FactorizedJetCorrector.h"
#include "JetCorrectorParameters.h"
#include "fastjet/contrib/ClusteringVetoPlugin.hh"
#include "fastjet/contrib/VariableRPlugin.hh"

//typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVectorXYZE;
class TopTagger{
 private:
  static  TopTagger* m_instance;
  mutable SLogger m_logger;


 public:
TopTagger();
  ~TopTagger();
 static TopTagger* Instance();
 bool Is_tagged(TString tagger, fastjet::PseudoJet jet);
 private:
 bool Is_CMS_tagged(fastjet::PseudoJet &jet);
 bool Is_HEP_tagged(fastjet::PseudoJet &jet);
 bool Is_OptimalR_tagged(fastjet::PseudoJet &jet);
 bool Is_HOTVR_tagged(fastjet::PseudoJet &jet);
 bool Is_softdrop_tagged(fastjet::PseudoJet &jet);
 fastjet::PseudoJet get_softdrop_jet(fastjet::PseudoJet jet, double zcut, double beta);
 
 

};



#endif
