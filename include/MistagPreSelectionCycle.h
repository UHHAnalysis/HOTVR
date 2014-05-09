// Dear emacs, this is -*- c++ -*-
#ifndef MistagPreSelectionCycle_H
#define MistagPreSelectionCycle_H

// SFrame include(s):
#include "include/AnalysisCycle.h"
#include "Cleaner.h"

/**
 *  @short Selection cycle to perform 
 *         pre-selection for Z'->ttbar analysis
 *  @author Thomas Peiffer
 */


class MistagPreSelectionCycle : public AnalysisCycle {

public:
  /// Default constructor
  MistagPreSelectionCycle();
  /// Default destructor
  ~MistagPreSelectionCycle();

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
  
  // Flg use to reverse electron selection
  bool m_reversed_electron_selection;
  double weight2=0;
  std::string m_Electron_Or_Muon_Selection;

  // Macro adding the functions for dictionary generation
  ClassDef( MistagPreSelectionCycle, 0 );

}; // class ZprimePreSelectionCycle

#endif // ZprimePreSelectionCycle_H

