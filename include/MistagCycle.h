// Dear emacs, this is -*- c++ -*-
// $Id: MyCycle.h,v 1.4 2013/06/12 12:40:14 peiffer Exp $
#ifndef MistagCycle_H
#define MistagCycle_H

// SFrame include(s):
//#include "../SFrameTools/include/TMVA_tagger.h"
#include "include/AnalysisCycle.h"
#include "include/SelectionModules.h"
#include "include/MyHists.h"
//#include "include/TopFitCalc.h"
#include "include/TopTagHists.h"
//#include "include/TopTagcontrol.h"
#include "include/Mistagcontrol.h"
#include "include/HypothesisHists.h"
#include "include/TMVA_tagger.h"
#include "EventHists.h"
#include "JetHists.h"
#include "ElectronHists.h"
#include "MuonHists.h"
#include "TauHists.h"
#include "TopJetHists.h"
#include "BTagEffHists.h"
#include "Cleaner.h"
#include "SubstructureHists.h"


/**
 *   @short My of an analysis cycle
 *
 *          This is an example of an analysis cycle. It can be used
 *          as a template for writing your own analysis. Also should
 *          be used for quick cross checks of the system setup.
 *
 *  @author Roman Kogler
 *  @version $Revision: 1.4 $
 */

class MistagCycle : public AnalysisCycle {

public:
  /// Default constructor
  MistagCycle();
  /// Default destructor
  ~MistagCycle();

  /// Function called at the beginning of the cycle
  void BeginCycle() throw( SError );
  /// Function called at the end of the cycle
  void EndCycle() throw( SError );

  /// Function called at the beginning of a new input data
  void BeginInputData( const SInputData& ) throw( SError );
  /// Function called after finishing to process an input data
  void EndInputData  ( const SInputData& ) throw( SError );

  /// Function called after opening each new input file
  void BeginInputFile( const SInputData& ) throw( SError );

  /// Function called for every event
  void ExecuteEvent( const SInputData&, Double_t ) throw( SError );

private:
  //
  // Put all your private variables here
  //

  // Macro adding the functions for dictionary generation
  ClassDef( MistagCycle, 0 );
   bool m_mttgencut;
   TMVA_tagger* tmva_tagger;
     TMVA::Reader* reader;
}; // class MyCycle

#endif // MyCycle_H

