#ifndef __FASTJET_CONTRIB_JETINFO_HH__
#define __FASTJET_CONTRIB_JETINFO_HH__
#include "fastjet/PseudoJet.hh"

FASTJET_BEGIN_NAMESPACE  

namespace contrib{

  class Jetinfo : public HOTVRinfo {
  public:
    //constructor
    //Jetinfo(int info);
    //Jetinfo(fastjet::PseudoJet jet,std::vector<PseudoJet> subjets): _subjets(subjets),_parent(jet) {}; //construct with the subjets and the the jet itself
    Jetinfo(fastjet::PseudoJet jet,std::vector<PseudoJet> subjets, double mass=-1);
    void test(){std::cout<<"test"<<std::endl;};
    
       
  private:
    //std::vector<PseudoJet> _subjets;
    //fastjet::PseudoJet _parent;
  };
  

}
FASTJET_END_NAMESPACE
#endif // __FASTJET_CONTRIB_JETINFO_HH__
