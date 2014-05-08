#include "include/TopTagfunctions.h"
#include "NtupleWriter/include/JetProps.h"
#include "include/MCDataScaleFactors.h"

#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/GhostedAreaSpec.hh>
namespace external {
#include "include/HEPTopTagger.h"
}
#include "SFrameTools/include/EventCalc.h"

using namespace std;

Double_t HelicityAngle(TLorentzVector j1, TLorentzVector j2)
{
  // calculate helicity angle between j1 and W, assuming j1 and j2 are it's decay products

TLorentzVector W = j1+j2;

TVector3 NegBoostOfW = -(W.BoostVector());

j1.Boost(NegBoostOfW);

 TVector3 j1_3 = j1.Vect();
 TVector3 j2_3 = j2.Vect();
 TVector3 W3 = W.Vect();

Double_t numerator = j1_3.Dot(W3);
Double_t denominator = (j1_3.Mag())*(W3.Mag());
Double_t CosHel = numerator/denominator;

Double_t HelAng = acos(CosHel);
 if (TMath::IsNaN(HelAng)) return -1;
 if (j2_3.Mag() == 0) return -1;
return HelAng;

}

bool jet_matches_top(TopJet topjet, float cone){
  
  EventCalc* calc = EventCalc::Instance();
  for(unsigned int k=0;k<calc->GetGenParticles()->size();++k){
	GenParticle top=calc->GetGenParticles()->at(k);
	if((abs(calc->GetGenParticles()->at(k).pdgId())==6) &&  (sqrt(pow(topjet.phi()-top.phi(),2)+pow(topjet.eta()-top.eta(),2))<cone)) return true;
	
	
  }
  return false;
}

bool jet_decays_hadronic(TopJet topjet, float cone) {

  EventCalc* calc = EventCalc::Instance();
   for(unsigned int k=0;k<calc->GetGenParticles()->size();++k){
     GenParticle top=calc->GetGenParticles()->at(k);
     if((abs(calc->GetGenParticles()->at(k).pdgId())==6) &&  (sqrt(pow(topjet.phi()-top.phi(),2)+pow(topjet.eta()-top.eta(),2))<cone))
       if(abs(top.daughter(calc->GetGenParticles(),1)->pdgId())==24) {
       if((abs(top.daughter(calc->GetGenParticles(),1)->daughter(calc->GetGenParticles(),1)->pdgId())<7)) return true;
       }
       else if((abs(top.daughter(calc->GetGenParticles(),2)->daughter(calc->GetGenParticles(),1)->pdgId())<7)) return true;
	
   }
   return false;
}

   bool jet_decays_leptonic(TopJet topjet, float cone) {

  EventCalc* calc = EventCalc::Instance();
   for(unsigned int k=0;k<calc->GetGenParticles()->size();++k){
     GenParticle top=calc->GetGenParticles()->at(k);
     if((abs(calc->GetGenParticles()->at(k).pdgId())==6) &&  (sqrt(pow(topjet.phi()-top.phi(),2)+pow(topjet.eta()-top.eta(),2))<cone))
       if(abs(top.daughter(calc->GetGenParticles(),1)->pdgId())==24) {
       if((abs(top.daughter(calc->GetGenParticles(),1)->daughter(calc->GetGenParticles(),1)->pdgId())>7)) return true;
       }
       else if((abs(top.daughter(calc->GetGenParticles(),2)->daughter(calc->GetGenParticles(),1)->pdgId())>7)) return true;
	
   }
   return false;
   }


 bool vec_matches_leptonic_decay(LorentzVector vec, float cone) {

  EventCalc* calc = EventCalc::Instance();
   for(unsigned int k=0;k<calc->GetGenParticles()->size();++k){
     GenParticle top=calc->GetGenParticles()->at(k);
     if((abs(calc->GetGenParticles()->at(k).pdgId())==6) &&  (sqrt(pow(vec.Phi()-top.phi(),2)+pow(vec.Eta()-top.eta(),2))<cone))
       if(abs(top.daughter(calc->GetGenParticles(),1)->pdgId())==24) {
       if((abs(top.daughter(calc->GetGenParticles(),1)->daughter(calc->GetGenParticles(),1)->pdgId())>7)) return true;
       }
       else if((abs(top.daughter(calc->GetGenParticles(),2)->daughter(calc->GetGenParticles(),1)->pdgId())>7)) return true;
	
   }
   return false;
   }

bool decay_products_in_jet_had(TopJet topjet, float cone){
      EventCalc* calc = EventCalc::Instance();
      for(unsigned int k=0;k<calc->GetGenParticles()->size();++k){
     GenParticle top=calc->GetGenParticles()->at(k);
     if(sqrt(pow(topjet.phi()-top.phi(),2)+pow(topjet.eta()-top.eta(),2))<0.8){
     if((abs(calc->GetGenParticles()->at(k).pdgId())==6)){
     const GenParticle *quark1;
     const GenParticle *quark2;
     const GenParticle *quark3;
     if(abs(top.daughter(calc->GetGenParticles(),1)->pdgId())==24){
       const GenParticle *quark1=top.daughter(calc->GetGenParticles(),2);
     const GenParticle *wboson=top.daughter(calc->GetGenParticles(),1);
     const GenParticle *quark2=top.daughter(calc->GetGenParticles(),1)->daughter(calc->GetGenParticles(),1);
     const GenParticle *quark3=top.daughter(calc->GetGenParticles(),1)->daughter(calc->GetGenParticles(),2);
    }
      else
       if(abs(top.daughter(calc->GetGenParticles(),2)->pdgId())==24){
	 const GenParticle *quark1=top.daughter(calc->GetGenParticles(),1);
     const GenParticle *wboson=top.daughter(calc->GetGenParticles(),2);
     const GenParticle *quark2=top.daughter(calc->GetGenParticles(),2)->daughter(calc->GetGenParticles(),1);
     const GenParticle *quark3=top.daughter(calc->GetGenParticles(),2)->daughter(calc->GetGenParticles(),2);
       }
     /* const GenParticle *quark2=top.daughter(calc->GetGenParticles(),2)->daughter(calc->GetGenParticles(),1);
	const GenParticle *quark3=wboson.daughter(calc->GetGenParticles(),2);*/
      if(sqrt(pow(topjet.phi()-quark1->phi(),2)+pow(topjet.eta()-quark1->eta(),2))<cone && sqrt(pow(topjet.phi()-quark2->phi(),2)+pow(topjet.eta()-quark2->eta(),2))<cone && sqrt(pow(topjet.phi()-quark3->phi(),2)+pow(topjet.eta()-quark3->eta(),2))<cone) return true;
     }}
	
      }

   
   return false;
}

double distance_quark(TopJet topjet, int nquark, float cone){
  EventCalc* calc = EventCalc::Instance();
      for(unsigned int k=0;k<calc->GetGenParticles()->size();++k){
     GenParticle top=calc->GetGenParticles()->at(k);
     if((abs(calc->GetGenParticles()->at(k).pdgId())==6)){
       if(sqrt(pow(topjet.phi()-top.phi(),2)+pow(topjet.eta()-top.eta(),2))<0.8){
     const GenParticle *quark1;
     const GenParticle *quark2;
     const GenParticle *quark3;
     const GenParticle *wboson;
     if(abs(top.daughter(calc->GetGenParticles(),1)->pdgId())==24){
     quark1=top.daughter(calc->GetGenParticles(),2);
     wboson=top.daughter(calc->GetGenParticles(),1);
     quark2=top.daughter(calc->GetGenParticles(),1)->daughter(calc->GetGenParticles(),1);
     quark3=top.daughter(calc->GetGenParticles(),1)->daughter(calc->GetGenParticles(),2);
    }
      else
       if(abs(top.daughter(calc->GetGenParticles(),2)->pdgId())==24){
     quark1=top.daughter(calc->GetGenParticles(),1);
     wboson=top.daughter(calc->GetGenParticles(),2);
     quark2=top.daughter(calc->GetGenParticles(),2)->daughter(calc->GetGenParticles(),1);
     quark3=top.daughter(calc->GetGenParticles(),2)->daughter(calc->GetGenParticles(),2);
       }
      if(nquark==1) return sqrt(pow(topjet.phi()-quark1->phi(),2)+pow(topjet.eta()-quark1->eta(),2));
     if(nquark==2) return sqrt(pow(topjet.phi()-quark2->phi(),2)+pow(topjet.eta()-quark2->eta(),2));
     if(nquark==3) return sqrt(pow(topjet.phi()-quark3->phi(),2)+pow(topjet.eta()-quark3->eta(),2));
     /*if(nquark==1) return quark1->eta();
     if(nquark==2) return quark2->eta();
     if(nquark==3) return quark3->eta();*/
       }}
      }
}

