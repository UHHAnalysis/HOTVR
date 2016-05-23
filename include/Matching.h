#ifndef Matching_H
#define Matching_H

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
#include <TStopwatch.h>
#include <fastjet/PseudoJet.hh>

//typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVectorXYZE;
class Matching{
 private:
  static  Matching* m_instance;
  mutable SLogger m_logger;



 private:
  std::vector<fastjet::PseudoJet> _hadrons;
  //std::vector<fastjet::PseudoJet> _partons;
  std::vector<fastjet::PseudoJet> _parton_jets;
  fastjet::ClusterSequence*  _clust_seq;
  TStopwatch timer;

  bool IsParton(GenParticle* p);
  GenParticle get_genparticle(GenParticle genparticle,std::vector<GenParticle>* genparticles, std::vector<int> &help);
  bool IsStableHadron(GenParticle* p);
  bool IsTop(GenParticle* p);
  bool FinalStateParton(GenParticle* p, std::vector<GenParticle>* genparticles);
  bool BeforeTopDecay(GenParticle* p, std::vector<GenParticle>* genparticles);
  GenParticle* GetDaughter(GenParticle* p, std::vector<GenParticle>* genparticles, int n);
  void UpdateSkipList(GenParticle* top, std::vector<GenParticle>* genparticles, std::vector<bool>& skip);
  fastjet::PseudoJet convert_particle(GenParticle* genparticle);
  std::vector<fastjet::PseudoJet> get_parton_jets(std::vector<fastjet::PseudoJet> parts);
  bool IsHadronic(GenParticle* p,  std::vector<GenParticle>* genparticles);
 public:
  void Run_matching(std::vector<GenParticle>* genparticles);
  Matching();
  ~Matching();
  static Matching* Instance();
  std::vector<fastjet::PseudoJet> get_hadrons(){return _hadrons;};
  bool IsMatched(fastjet::PseudoJet jet, double matching_distance, fastjet::PseudoJet denominator_jet); 
  fastjet::PseudoJet get_closest_jet(std::vector<fastjet::PseudoJet> jets,fastjet::PseudoJet denominator_jet);
  std::vector<fastjet::PseudoJet> get_denominator_jets(TString idVersion);
 // std::vector<fastjet::PseudoJet> get_genjets();



};

#endif
