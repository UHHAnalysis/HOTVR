#ifndef MyHists_H
#define MyHists_H

// ROOT include(s):
#include <TObject.h>
#include <TString.h>
#include "TH2.h"

// Local include(s):
#include "include/BaseHists.h"
#include "include/BaseCycleContainer.h"

/**
 *   My class for booking and filling histograms
 *
 *   This class books and fills a collection of histograms.
 *   It should have a unique name, such that the histograms
 *   of multiple instances of this class are ordered in the
 *   output file. 
 *   Always sort your histograms and used methods topically.
 *   This example collection can be used for data and reconstructed
 *   MC events.
 *   
 *   @version $Revision: 1.2 $
 */

class MyHists : public BaseHists {

public:
   /// Named constructor
   MyHists(const char* name);

   /// Default destructor
   ~MyHists();

   void Init();

   void Fill();

   void Finish();

private:

}; // class MyHists


#endif // MyHists_H
