
#ifndef SubstructureHists_H
#define SubstructureHists_H

// ROOT include(s):
#include <TObject.h>
#include <TString.h>
#include "TH2.h"

// Local include(s):
#include "SFrameTools/include/BaseHists.h"
#include "include/Matching.h"
#include "include/Clustering.h"
 #include "fastjet/PseudoJet.hh"
class SubstructureHists : public BaseHists {

public:
   /// Named constructor
  SubstructureHists(const char* name);

   /// Default destructor
   ~SubstructureHists();

   void Init();
   void Fill(){};
   
   void Fill(fastjet::PseudoJet jet,double jet_radius, fastjet::PseudoJet matched_jet, double weight);

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
   Matching* matching;
   Clustering* clustering;
}; // class SubstructureHists


#endif // Substructure_H
