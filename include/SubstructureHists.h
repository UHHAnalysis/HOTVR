#ifndef SubstructureHists_H
#define SubstructureHists_H
#include "SFrameTools/include/BaseHists.h"
#include "SFrameTools/include/Objects.h"
#include "SFrameTools/include/fwd.h"
#include "SFrameTools/include/boost_includes.h"
#include "include/BaseHists.h"
#include "include/BaseCycleContainer.h"
#include "HypothesisDiscriminator.h"
#include "TH2.h"
#include "Utils.h"
#include "EventCalc.h"
//#include "include/TopFitCalc.h"
/**
 *   Class for booking and filling TopJet histograms
 *
 *   
 *   @version $Revision: 1.1 $
 */

class SubstructureHists : public BaseHists {

public:
   /// Named constructor
   SubstructureHists(const char* name);

   /// Default destructor
   ~SubstructureHists();

   void Init();
   void Fill();
   void Fill2(TopJet topjet, double mva_value);
 // class TopTagHists
};
#endif // TopTagcontrol_H


