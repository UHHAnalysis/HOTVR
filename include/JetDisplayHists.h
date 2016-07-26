
#ifndef JetDisplayHists_H
#define JetDisplayHists_H

// ROOT include(s):
#include <TObject.h>
#include <TString.h>
#include "TH2.h"

// Local include(s):
#include "SFrameTools/include/BaseHists.h"
#include "include/Infrared_Saftey.h"
#include "include/Clustering.h"

class JetDisplayHists : public BaseHists {
 

public:
   /// Named constructor
  JetDisplayHists(const char* name);

   /// Default destructor
   ~JetDisplayHists();

   void Init();

   void Fill();
   
 
   void Finish();
   
   void SetIdVersion(TString s);
   bool FillEvent(std::vector<fastjet::PseudoJet> jets,std::vector<fastjet::PseudoJet> parts, std::vector<fastjet::PseudoJet> soft_jets, std::vector<fastjet::PseudoJet> rejected_jets, std::vector<fastjet::PseudoJet> rejected_subjets);

private:
 

   TString idVersion;
   Infrared_Saftey* IR_Saftey;
   Clustering* clustering;
}; // class HOTVRHists


#endif // HOTVRHists_H
