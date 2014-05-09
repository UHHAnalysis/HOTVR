#include <iostream>

using namespace std;

// Local include(s):
#include "include/MistagPreSelectionCycle.h"
#include "include/SelectionModules.h"
#include "JetCorrectorParameters.h"

ClassImp( MistagPreSelectionCycle );

MistagPreSelectionCycle::MistagPreSelectionCycle()
    : AnalysisCycle()
{

    // constructor, declare additional variables that should be
    // obtained from the steering-xml file

    // set the integrated luminosity per bin for the lumi-yield control plots
    SetIntLumiPerBin(25.);

    DeclareProperty( "Electron_Or_Muon_Selection", m_Electron_Or_Muon_Selection );

    // steering property for data-driven qcd in electron channel
    m_reversed_electron_selection = false;
    DeclareProperty( "ReversedElectronSelection", m_reversed_electron_selection);
}

MistagPreSelectionCycle::~MistagPreSelectionCycle()
{
    // destructor

}

void MistagPreSelectionCycle::BeginCycle() throw( SError )
{
    // Start of the job, general set-up and definition of
    // objects are done here

    // Important: first call BeginCycle of base class
    AnalysisCycle::BeginCycle();

    return;

}

void MistagPreSelectionCycle::EndCycle() throw( SError )
{
    // clean-up, info messages and final calculations after the analysis


    // call the base cycle class for all standard methods
    // and a summary of the made selections
    AnalysisCycle::EndCycle();

    return;

}

void MistagPreSelectionCycle::BeginInputData( const SInputData& id ) throw( SError )
{
    // declaration of histograms and selections

    // Important: first call BeginInputData of base class
    AnalysisCycle::BeginInputData( id );
    // -------------------- set up the selections ---------------------------

    Selection* preselection = new Selection("preselection");
    Book( TH1F( "weight","sum of all weights",3,0,2));
    if(m_Electron_Or_Muon_Selection=="Electrons" || m_Electron_Or_Muon_Selection=="Electron" || m_Electron_Or_Muon_Selection=="Ele" || m_Electron_Or_Muon_Selection=="ELE") {
        preselection->addSelectionModule(new NElectronSelection(0,0));//at least one electron
        preselection->addSelectionModule(new NMuonSelection(0,0));//no muons
    } else if(m_Electron_Or_Muon_Selection=="Muon" || m_Electron_Or_Muon_Selection=="Muons" || m_Electron_Or_Muon_Selection=="Mu" || m_Electron_Or_Muon_Selection=="MU") {
        preselection->addSelectionModule(new NElectronSelection(0,0));//no electron
        preselection->addSelectionModule(new NMuonSelection(0,0));//at least one muon
    } else {
        m_logger << ERROR << "Electron_Or_Muon_Selection is not defined in your xml config file --- should be either `ELE` or `MU`" << SLogger::endmsg;
    }

    preselection->addSelectionModule(new NJetSelection(2,int_infinity(),150));//at least two jets
     preselection->addSelectionModule(new NJetSelection(1,int_infinity(),320));
     preselection->addSelectionModule(new TriggerSelection("HLT_HT750_v"));
      preselection->addSelectionModule(new HThadCut(0,int_infinity(),1000,int_infinity()));
      //cuts from paper
      // preselection->addSelectionModule(new NJetSelection(2,int_infinity(),400));
       //  preselection->addSelectionModule(new NJetdeltaySelection(1.0));
       preselection->addSelectionModule(new NTopJetSelection(2,int_infinity(),400,int_infinity()));
       preselection->addSelectionModule(new NJetdeltaphiSelection(2.1));
    // preselection->addSelectionModule(new NleadingJetSelection(320));//leading Jet pt

    RegisterSelection(preselection);

    return;
}

void MistagPreSelectionCycle::EndInputData( const SInputData& id ) throw( SError )
{
    AnalysisCycle::EndInputData( id );
  std::cout<<weight2<<std::endl;
    return;
}

void MistagPreSelectionCycle::BeginInputFile( const SInputData& id ) throw( SError )
{
    // Connect all variables from the Ntuple file with the ones needed for the analysis
    // The variables are commonly stored in the BaseCycleContaincer

    // important: call to base function to connect all variables to Ntuples from the input tree
    AnalysisCycle::BeginInputFile( id );

    return;
}

void MistagPreSelectionCycle::ExecuteEvent( const SInputData& id, Double_t weight) throw( SError )
{
    // this is the most important part: here the full analysis happens
    // user should implement selections, filling of histograms and results

    // first step: call Execute event of base class to perform basic consistency checks
    // also, the good-run selection is performed there and the calculator is reset
    AnalysisCycle::ExecuteEvent( id, weight);

    Cleaner cleaner;
    static Selection* preselection = GetSelection("preselection");
    EventCalc* calc = EventCalc::Instance();
    //double weight2;
    weight2+=calc->GetWeight();
    Hist("weight")->Fill(1,calc->GetWeight());
    // std::cout<<calc->GetWeight()<<std::endl;
    BaseCycleContainer* bcc = calc->GetBaseCycleContainer();

    //save uncleaned jet collection and MET to be stored in output
    std::vector<Jet> uncleaned_jets;
    for(unsigned int i=0; i<bcc->jets->size(); ++i) {
        uncleaned_jets.push_back(bcc->jets->at(i));
    }
    MET uncleaned_met = *bcc->met;

    //save uncleaned muon collection to be stored in output
    std::vector<Muon> uncleaned_muons;
    for(unsigned int i=0; i<bcc->muons->size(); ++i) {
        uncleaned_muons.push_back(bcc->muons->at(i));
    }

    //save uncleaned electron collection to be stored in output
    std::vector<Electron> uncleaned_electrons;
    for(unsigned int i=0; i<bcc->electrons->size(); ++i) {
        uncleaned_electrons.push_back(bcc->electrons->at(i));
    }


    //clean collections here

    if(bcc->muons) cleaner.MuonCleaner_noIso(45,2.1);
    if(bcc->electrons) cleaner.ElectronCleaner_noIso(35,2.5,m_reversed_electron_selection);
    if(bcc->taus) cleaner.TauCleaner(20,2.1);
    if(bcc->jets) cleaner.JetLeptonSubtractor(m_corrector,false);
    if(!bcc->isRealData && bcc->jets) cleaner.JetEnergyResolutionShifter();
    if(bcc->jets) cleaner.JetCleaner(30,2.5,true);

    //get the jet selections

    if(!preselection->passSelection())  throw SError( SError::SkipEvent );

    //fill the uncleaned collections back to bcc to store them in output tree
    bcc->met->set_pt (uncleaned_met.pt());
    bcc->met->set_phi (uncleaned_met.phi());
    bcc->jets->clear();
    for(unsigned int i=0; i<uncleaned_jets.size(); ++i) {
        bcc->jets->push_back(uncleaned_jets.at(i));
    }
    bcc->muons->clear();
    for(unsigned int i=0; i<uncleaned_muons.size(); ++i) {
        bcc->muons->push_back(uncleaned_muons.at(i));
    }
    bcc->electrons->clear();
    for(unsigned int i=0; i<uncleaned_electrons.size(); ++i) {
        bcc->electrons->push_back(uncleaned_electrons.at(i));
    }

    WriteOutputTree();

    return;
}

