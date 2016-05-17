#ifndef Infrared_Saftey_H
#define Infrared_Saftey_H

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


//typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVectorXYZE;
class Infrared_Saftey{
 private:
  static  Infrared_Saftey* m_instance;
  mutable SLogger m_logger;


 public:
Infrared_Saftey();
  ~Infrared_Saftey();
 static Infrared_Saftey* Instance();

 void add_soft_particles(std::vector<fastjet::PseudoJet> &jetpart);
 void add_grid(std::vector<fastjet::PseudoJet> &jetpart,int phi, int eta); 
double create_random(double low, double high);
 void  add_collinear_splitting(std::vector<fastjet::PseudoJet> &jetpart);
};

#endif
