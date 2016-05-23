#ifndef Clustering_H
#define Clustering_H

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
class Clustering{
 private:
  static  Clustering* m_instance;
  mutable SLogger m_logger;
  fastjet::ClusterSequence* _clust_seq2;
  std::vector<fastjet::PseudoJet> _fatjets;
 
 public:
 enum E_algorithm { 
    e_ca, 
    e_akt,
    e_kt,
    e_varR_ca,
    e_varR_akt,
    e_varR_kt,
    e_massjump_ca,
    e_massjump_akt,
    e_massjump_kt,
   };
 E_algorithm m_test;

Clustering();
  ~Clustering();
 static Clustering* Instance();
 void Reset();
  

 //usual squentiell clustering
  std::vector<fastjet::PseudoJet> get_clustered_jets(std::vector<fastjet::PseudoJet> particles,enum Clustering::E_algorithm algorithm, double jet_radius, double ptmin);
  //variableR algorithm
  std::vector<fastjet::PseudoJet> get_clustered_jets(std::vector<fastjet::PseudoJet> particles,enum  Clustering::E_algorithm algorithm, double ptmin, double rho, double min_r, double max_r);
  //  std::vector<fastjet::PseudoJet> get_clustered_jets(std::vector<fastjet::PseudoJet> particles,enum  Clustering::E_algorithm algorithm, double jet_radius,double ptmin, double rho, double min_r, double max_r);
  std::vector<fastjet::PseudoJet> get_clustered_hotvr_jets(std::vector<fastjet::PseudoJet> particles,enum  Clustering::E_algorithm algorithm, double ptmin, double rho, double min_r, double max_r, double mu, double theta, double pt_cut);

};

#endif
