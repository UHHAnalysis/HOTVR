
#ifndef ClusteringHists_H
#define ClusteringHists_H

// ROOT include(s):
#include <TObject.h>
#include <TString.h>
#include "TH2.h"

// Local include(s):
#include "SFrameTools/include/BaseHists.h"
 #include "fastjet/PseudoJet.hh"
class ClusteringHists : public BaseHists {

public:
   /// Named constructor
  ClusteringHists(const char* name);

   /// Default destructor
   ~ClusteringHists();

   void Init();
   void Fill(){};
   
   void Fill(fastjet::PseudoJet jet, double weight);

   void Finish();
   
  

private:
  
}; // class ClusteringHists


#endif // Clustering_H
