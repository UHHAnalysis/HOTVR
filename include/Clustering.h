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
  fastjet::PseudoJet _setup_jet;
  static  Clustering* m_instance;
  mutable SLogger m_logger;
  fastjet::ClusterSequence* _clust_seq2;
  std::vector<fastjet::PseudoJet> _fatjets;
  double get_value(std::string word);
  void set_modus(std::string word);
  void show_settings();
  std::string get_modus();
  std::vector<fastjet::PseudoJet> set_properties(std::vector<fastjet::PseudoJet> jets);
 public:
 enum E_algorithm { 
    e_ca, 
    e_akt,
    e_kt,
   };
 E_algorithm m_test;
 Clustering();
 Clustering(std::string clustering);
  ~Clustering();
 static Clustering* Instance();
 void Reset();
  
 fastjet::PseudoJet get_settings(){return _setup_jet;};
 //usual squentiell clustering
 std::vector<fastjet::PseudoJet> get_clustered_jets(std::vector<fastjet::PseudoJet> particles);
  std::vector<fastjet::PseudoJet> get_clustered_jets(std::vector<fastjet::PseudoJet> particles,enum Clustering::E_algorithm algorithm, double jet_radius, double ptmin);
  //variableR algorithm
  std::vector<fastjet::PseudoJet> get_clustered_jets(std::vector<fastjet::PseudoJet> particles,enum  Clustering::E_algorithm algorithm, double ptmin, double rho, double min_r, double max_r);
  //  std::vector<fastjet::PseudoJet> get_clustered_jets(std::vector<fastjet::PseudoJet> particles,enum  Clustering::E_algorithm algorithm, double jet_radius,double ptmin, double rho, double min_r, double max_r);
  std::vector<fastjet::PseudoJet> get_clustered_hotvr_jets(std::vector<fastjet::PseudoJet> particles,enum  Clustering::E_algorithm algorithm, double ptmin, double rho, double min_r, double max_r, double mu, double theta, double pt_cut);
  void jets(std::string clustering);

  double _rho, _mu, _theta, _rmin, _rmax, _ptmin,_radius,_pt_cut;
  string _clustering_algorithmus;
 
};

#endif
