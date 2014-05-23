#ifndef Showerdeconstruction_H
#define Showerdeconstruction_H

#include "SFrameTools/include/Objects.h"
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


//typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVectorXYZE;
class Showerdeconstruction{
 private:
  static Showerdeconstruction* m_instance;
  mutable SLogger m_logger;
 TMVA::Reader *reader=NULL;


 public:
Showerdeconstruction();
  ~Showerdeconstruction();
  double  Chi(TopJet topjet);
  double ChiMicro(TopJet topjet);
  int GetNmicrojets(TopJet topjet);
static Showerdeconstruction* Instance();


};

#endif
