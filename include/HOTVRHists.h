
#ifndef HOTVRHists_H
#define HOTVRHists_H

// ROOT include(s):
#include <TObject.h>
#include <TString.h>
#include "TH2.h"

// Local include(s):
#include "SFrameTools/include/BaseHists.h"


class HOTVRHists : public BaseHists {

public:
   /// Named constructor
  HOTVRHists(const char* name, double mass_min, double mass_max, double mmin_min, double mmin_max, int Nsubjets_min, int Nsubjets_max, double ptfraction_min, double ptfraction_max, double subjet1pt, double subjet2pt);

   /// Default destructor
   ~HOTVRHists();

   void Init();

   void Fill();

   void Finish();
   
   void SetIdVersion(TString s);

private:
   double m_mass_min;
   double m_mass_max;
   double m_mmin_min;
   double m_mmin_max;
   int m_Nsubjets_min;
   int m_Nsubjets_max;
   double m_ptfraction_min;
   double m_ptfraction_max;
   double m_subjet1pt;
   double m_subjet2pt;


   TString idVersion;
 
}; // class HOTVRHists


#endif // HOTVRHists_H
