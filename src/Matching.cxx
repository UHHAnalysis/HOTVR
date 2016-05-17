#include "include/Matching.h"
//#include "NtupleWriter/include/JetProps.h"
//#include "NtupleWriter/interface/GenJetProps.h"
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include "SFrameTools/include/EventCalc.h"
#include "TRandom3.h"
#include <algorithm> 



using namespace std;

Matching* Matching::m_instance = NULL;

Matching* Matching::Instance()
{
  // Get a pointer to the object handler.
  // This is the only way to access this class, 
  // since it's a singleton. This method is accessible
  // from everywhere.
 
  if (m_instance == NULL){
    m_instance = new  Matching();
  }
  return m_instance; 
   
}


Matching::Matching() : m_logger("Matching")
{
  // constructor: initialise all variables

  
 
     //Reset();

  
}



Matching::~Matching()
{
  // default destructor
  std::cout<<"destructor called"<<std::endl;
  delete m_instance;
}

fastjet::PseudoJet Matching::convert_particle(GenParticle* genparticle){
  TLorentzVector particle;
  particle.SetPtEtaPhiE(genparticle->pt(),genparticle->eta(),genparticle->phi(),genparticle->energy());
  fastjet::PseudoJet gen_particle(particle.Px(),particle.Py(),particle.Pz(),particle.E());
  return gen_particle;

}



bool Matching::IsParton(GenParticle* p)
{
	// keep all quarks, gluons, W, Z bosons
	int fid = std::abs(p->pdgId());
	bool isp = false;
	if (fid<7)   isp = true;
	if (fid==21) isp = true;
	if (fid==23) isp = true;
	if (fid==24) isp = true;
 	return isp;
}

bool Matching::IsStableHadron(GenParticle* p)
{
	int st = p->status();
	if (st==1) return true;
	else return false;
}

bool Matching::IsTop(GenParticle* p)
{
	int id = p->pdgId();
	if (std::abs(id)==6) return true;
	else return false;
}

bool Matching::FinalStateParton(GenParticle* p, std::vector<GenParticle>* genparticles)
{
	// Tobias: what was the definition of a final state parton? 
	// no daughter? please check!!
	GenParticle* d1 = GetDaughter(p, genparticles, 1);
	GenParticle* d2 = GetDaughter(p, genparticles, 2);
	bool result = false;
	// daughters found?
	
	if (!d1 && !d2 && std::abs(p->pdgId())<100) result = true;

	if (d1){
	  int fid = std::abs(d1->pdgId());
	  if (fid>100) result = true;
	}
	if (d2){
	  int fid = std::abs(d2->pdgId());
	  if (fid>100) result = true;
	}
	if(std::abs(p->pdgId())>100) result=false;

	// protect against beam remnant
	LorentzVector v = p->v4();
	if (std::abs(v.Eta()) > 8) result = false;
	if (std::isinf(v.Eta())) result = false;
	if (std::isnan(v.Eta())) result = false;

	return result;
}

bool Matching::BeforeTopDecay(GenParticle* p, std::vector<GenParticle>* genparticles)
{
	// check if this top quark decays into W+b
	// get daughters
	GenParticle* d1 = GetDaughter(p, genparticles, 1);
	GenParticle* d2 = GetDaughter(p, genparticles, 2);

	// daughters found?
	if (!d1) return false;
	if (!d2) return false;

	bool flag = false;
	// check ID of daughters
	int id1 = std::abs(d1->pdgId());
	int id2 = std::abs(d2->pdgId());	

	if (id1==5  && id2==24) flag = true;
	if (id1==24 && id2==5)  flag = true; 

	return flag;
}

GenParticle* Matching::GetDaughter(GenParticle* p, std::vector<GenParticle>* genparticles, int n)
{
	int index = 0;
	GenParticle* d = NULL;
	
	if (n<=1){
		index = p->daughter1();
	} else {
		index = p->daughter2();		
	}

	if (index>0 && index<genparticles->size()){
		d = &genparticles->at(index-1);
	}
	return d;
}

void Matching::UpdateSkipList(GenParticle* top, std::vector<GenParticle>* genparticles, std::vector<bool>& skip)
{
	// go through all decays from the top and flag those to be skipped
	GenParticle* d1 = GetDaughter(top, genparticles, 1);
	GenParticle* d2 = GetDaughter(top, genparticles, 2);

	std::vector<GenParticle*> daughters;
	if (d1){ 
		skip[d1->index()-1] = true;
		daughters.push_back(d1);
	}
	if (d2){ 
		skip[d2->index()-1] = true;
		daughters.push_back(d2);
	}

	// need to check if this works
	int nd = daughters.size();
	for (int i=0; i<nd; ++i){
		GenParticle* p = daughters[i];
		GenParticle* d1 = GetDaughter(p, genparticles, 1);
		GenParticle* d2 = GetDaughter(p, genparticles, 2);
		if (d1){ 
		  if (!IsParton(d1)) continue;
			skip[d1->index()-1] = true;
			daughters.push_back(d1);
		}
		if (d2){ 
		  if (!IsParton(d2)) continue;
			skip[d2->index()-1] = true;
			daughters.push_back(d2);
		}
		nd = daughters.size();
	}
}

std::vector<fastjet::PseudoJet> Matching::get_parton_jets(std::vector<fastjet::PseudoJet> parts){
  std::cout<<"cluster"<<std::endl;
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm,0.4);
  std::vector<fastjet::PseudoJet> fatjets;
  fastjet::ClusterSequence*  clust_seq;
  clust_seq= new fastjet::ClusterSequence(parts, jet_def);
  // _clust_seq=new fastjet::ClusterSequence(parts, jet_def);
  fatjets = sorted_by_pt(clust_seq->inclusive_jets(100.));
  //delete  clust_seq;
  return fatjets;
}

std::vector<fastjet::PseudoJet> Matching::get_denominator_jets(TString idVersion){
  std::vector<fastjet::PseudoJet> denominator_jets;
  if(idVersion.Contains("QCD")) denominator_jets=_parton_jets;
  if(idVersion.Contains("TTbar")) {
    for(int i=0;i<_parton_jets.size();i++) {
      bool candidate=false;
      for(int j=0;j<_parton_jets[i].constituents().size();j++)	if(_parton_jets.at(i).constituents().at(j).user_index()==6) candidate=true;
      if(candidate) {
	denominator_jets.push_back(_parton_jets[i]);
      }
    }
  }
   return denominator_jets;
}


fastjet::PseudoJet Matching::get_closest_jet(std::vector<fastjet::PseudoJet> jets,fastjet::PseudoJet denominator_jet){
  fastjet::PseudoJet matched_jet;   
double minDeltaR=1000;
  double deltaR;
for(int i=0;i<jets.size();i++){
    deltaR=jets[i].delta_R(denominator_jet);
    if(deltaR<minDeltaR) {
      minDeltaR=deltaR;
      matched_jet=jets[i];
    }
 }
 return matched_jet;

}

bool Matching::IsMatched(fastjet::PseudoJet jet, double matching_radius,fastjet::PseudoJet denominator_jet){
 
  //denominator_jets=get_denominator_jets(idVersion);
  double deltaR;
  deltaR=jet.delta_R(denominator_jet);
  if(deltaR<matching_radius) return true;
  else return false;

}


void Matching::Run_matching(std::vector<GenParticle>* genparticles){
  // EventCalc* calc = EventCalc::Instance();
  // std::vector<GenParticle>* genparticles = calc->GetGenParticles();
 
   
  // std::vector<GenParticle*> _hadrons;
  
  // loop over all gen-particles and fill stable hadrons into list
  // (use pointers such that the particles themselves are not copied to save CPU time)
  for (int i=0; i<genparticles->size(); ++i)
    { 
      // get pointer to particle
      GenParticle* part = &(genparticles->at(i));
           // store stable hadrons
	if (IsStableHadron(part)){
	  _hadrons.push_back(convert_particle(part));
	  continue;
	}
	
}
  
  // Tobias: now you can use the array with hadrons for clustering
  
  // now get partonic final state
  
// loop over particles, store only the ones before the hadronisation
// also separate out the ones from the top decays
  
// here are the final state partons
  //std::vector<fastjet::PseudoJet> partons;
  std::vector<GenParticle*> partons;

  // here is a list of partons that should be skipped (daughters from top decay)

  std::vector<bool> skip;
  int nparts = genparticles->size();
  skip.reserve(nparts);

   for (int i=0; i<nparts; ++i) skip[i]=false;
 
   // loop over all particles
   for (int i=0; i<nparts; ++i)
     {
       // check if this parton should be skipped
       if (skip[i]) continue;
	
       GenParticle* p = &(genparticles->at(i));
	
       // is it a parton?
       if (!IsParton(p)) continue;

       // check if it's a top quark
       if (IsTop(p)){
	  	  
	 // is it the one that decays?
	 if (BeforeTopDecay(p,genparticles)){
	    
	   // keep it 
	   partons.push_back(p);

	   // flag all daughter partons to be skipped
	   UpdateSkipList(p, genparticles, skip);
			
	   continue;
	 } else {
	   
	   continue;
	    
	 }
       
       }
	 
       if (FinalStateParton(p, genparticles)){

	 partons.push_back(p);

       } else {

	 // throw away

       }

       continue;

     }

   std::vector<fastjet::PseudoJet> partons_to_cluster;
   LorentzVector s(0,0,0,0);

   for(int i;i<partons.size();i++){
     GenParticle* p = partons[i];
     fastjet::PseudoJet particle=convert_particle(p);
     if(std::abs(partons[i]->pdgId())==6) particle.set_user_index(6);
     else  particle.set_user_index(0);
     partons_to_cluster.push_back(particle); 
     //std::cout<<"partons to cluster "<<p->pdgId()<<std::endl;
     s += p->v4();
   }
   
  


    _parton_jets=get_parton_jets(partons_to_cluster);
    
       
    //cout << "s: px = " << s.px() << " py = " << s.py() << " pz = " << s.pz() << " E = " << s.E() << endl;

  

  

 
 

 
}
//genparticles->at(tx).index()
