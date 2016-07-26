// Dear emacs, this is -*- c++ -*-
// $Id: RocCycle.h,v 1.4 2013/06/12 12:40:14 peiffer Exp $
#ifndef JetDisplayCycle_H
#define JetDisplayCycle_H

#include "SFrameAnalysis/include/AnalysisCycle.h"
#include "include/Matching.h"
#include "include/Clustering.h"
#include "include/TopTagger.h"
#include "include/Infrared_Saftey.h"
//#include "ZprimeFullHadTools.h"
//#include "include/SubstructureHists.h"

/**
 *   @short ZprimeFullHad of an analysis cycle
 *
 *          This is an example of an analysis cycle. It can be used
 *          as a template for writing your own analysis. Also should
 *          be used for quick cross checks of the system setup.
 *
 *  @author Roman Kogler
 *  @version $Revision: 1.4 $
 */

class JetDisplayCycle : public AnalysisCycle {

public:
  /// Default constructor
  JetDisplayCycle();
  /// Default destructor
  ~JetDisplayCycle();

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
  int m_counter;
  int m_counter2;
  Matching* matching;
  Clustering* clustering, *clustering_ca, *clustering_akt;
  TopTagger* toptagger;
  Infrared_Saftey* IR_Saftey;
  BaseHists* m_Jetdisplay_hists_event[100];
  BaseHists* m_Jetdisplay_akt_hists_event[100];
  BaseHists* m_Jetdisplay_ca_hists_event[100];
  BaseHists* m_Jetdisplay_hists_highpT_event[100];
  BaseHists* m_Jetdisplay_akt_hists_highpT_event[100];
  BaseHists* m_Jetdisplay_ca_hists_highpT_event[100];
   //Selection *TriggerSel;//* BSel, * NoBSel, *TopSel, *chi2_selection;
  int _event_max=50;
  // Macro adding the functions for dictionary generation
  ClassDef( JetDisplayCycle, 0 );

}; // class JetDisplayCycle

#endif // JetDisplayCycle_H

