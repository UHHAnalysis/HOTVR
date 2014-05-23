#include "include/Showerdeconstruction.h"
#include "NtupleWriter/include/JetProps.h"

#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include "SFrameTools/include/EventCalc.h"

#include "include/TopTagfunctions.h"
#include "include/topquark.h"

using namespace std;

Showerdeconstruction* Showerdeconstruction::m_instance = NULL;

Showerdeconstruction* Showerdeconstruction::Instance()
{
  // Get a pointer to the object handler.
  // This is the only way to access this class, 
  // since it's a singleton. This method is accessible
  // from everywhere.

  if (m_instance == NULL){
    m_instance = new Showerdeconstruction();
  }
  return m_instance; 
   
}


Showerdeconstruction::Showerdeconstruction() : m_logger("TMVA_tagger")
{
  // constructor: initialise all variables

  
  m_logger << DEBUG << "Constructor called." << SLogger::endmsg;
  
     //Reset();

  
}



Showerdeconstruction::~Showerdeconstruction()
{
  // default destructor
  
}

double Showerdeconstruction::Chi(TopJet topjet){
  double Psignal = 0.0;
  double Pbackground = 0.0;
  DES::fourvector pTop;
  std::vector<Particle> subjets = topjet.subjets();
  sort(subjets.begin(),subjets.end(),HigherPt());
  int nsubjets = subjets.size();
  DES::microjet mJ[DES::maxsize];
  for(unsigned int j=0;j<nsubjets;++j){
    Particle subjet = subjets[j];
    LorentzVector microjet = subjet.v4();
    DES::microjet mjet(microjet.E(),microjet.px(),microjet.py(),microjet.pz());
    /* std::cout<<"E: "<<microjet.E()<<std::endl;
    std::cout<<"px: "<<microjet.Px()<<std::endl;
    std::cout<<"py: "<<microjet.Py()<<std::endl;
    std::cout<<"pz: "<<microjet.Pz()<<std::endl;*/
    mJ[j] = mjet; 
  }
  double chi = DES::Chi(mJ,nsubjets,Psignal,Pbackground, pTop);
  // if(chi!=0) std::cout<<"ca subjets "<<Psignal<<" "<<Pbackground<<"Chi "<<chi<<" "<<nsubjets<<std::endl;
  return chi;
}/*

double Showerdeconstruction::ChiMicro(TopJet topjet){
   EventCalc* calc = EventCalc::Instance();
  JetProps jp(&topjet, calc->GetPFParticles() );
  std::cout<<"test1"<<std::endl;
  std::vector<Particle> microjets = jp.GetMicroSubjet(0.2);
  //same as in Chi
  double Psignal = 0.0;
  double Pbackground = 0.0;
  DES::fourvector pTop;
   std::cout<<"test2"<<std::endl;
  int nsubjets = microjets.size();
  if(nsubjets>9) nsubjets=9;
  DES::microjet mJ[DES::maxsize];
   std::cout<<"test3"<<std::endl;
  for(unsigned int j=0;j<nsubjets;++j){
    Particle subjet = microjets[j];
    LorentzVector microjet = subjet.v4();
    DES::microjet mjet(microjet.E(),microjet.px(),microjet.py(),microjet.pz());
    mJ[j] = mjet; 
     std::cout<<"test4"<<std::endl;
  }
  double chi = DES::Chi(mJ,nsubjets,Psignal,Pbackground, pTop);
  return chi;
  return 0;
  }*/


double Showerdeconstruction::ChiMicro(TopJet topjet){
  double Psignal = 0.0;
  double Pbackground = 0.0;
  DES::fourvector pTop;
  EventCalc* calc = EventCalc::Instance();
  JetProps jp(&topjet, calc->GetPFParticles() );
  std::vector<fastjet::PseudoJet> pseudojets = jp.GetFastJet(0.2);
  int nsubjets=pseudojets.size();
  if(nsubjets>9) nsubjets=9;
  bool test=false;
  if(nsubjets>5) test=true;
   DES::microjet mJ[DES::maxsize];
   for(unsigned int j=0;j<nsubjets;++j){
     fastjet::PseudoJet pseudojet=pseudojets[j];
     DES::microjet mjet(pseudojet.E(),pseudojet.px(),pseudojet.py(),pseudojet.pz());
     if(test){
       /*  std::cout<<"E: "<<pseudojet.E()<<std::endl;
     std::cout<<"px: "<<pseudojet.px()<<std::endl;
     std::cout<<"py: "<<pseudojet.py()<<std::endl;
     std::cout<<"pz: "<<pseudojet.pz()<<std::endl;*/
     }
      mJ[j] = mjet;
   }
double chi = DES::Chi(mJ,nsubjets,Psignal,Pbackground, pTop);
//if(chi!=0) std::cout<<Psignal<<" "<<Pbackground<<"Chi "<<chi<<" "<<nsubjets<<std::endl;
  return chi;
}


int Showerdeconstruction::GetNmicrojets(TopJet topjet){
  EventCalc* calc = EventCalc::Instance();
  JetProps jp(&topjet, calc->GetPFParticles() );
  std::vector<fastjet::PseudoJet> pseudojets = jp.GetFastJet(0.2);
  int nmicrojets=pseudojets.size();
  return nmicrojets;

}
