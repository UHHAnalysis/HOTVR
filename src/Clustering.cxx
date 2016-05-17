#include "include/Clustering.h"
//#include "NtupleWriter/include/JetProps.h"
//#include "NtupleWriter/interface/GenJetProps.h"
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include "SFrameTools/include/EventCalc.h"
#include "fastjet/contrib/HOTVR.hh"

using namespace std;
using namespace fastjet;
using namespace contrib;

Clustering* Clustering::m_instance = NULL;

Clustering* Clustering::Instance()
{
  // Get a pointer to the object handler.
  // This is the only way to access this class, 
  // since it's a singleton. This method is accessible
  // from everywhere.

  if (m_instance == NULL){
    m_instance = new  Clustering();
  }
  return m_instance; 
   
}


Clustering::Clustering() : m_logger("Clustering")
{
  // constructor: initialise all variables

  
 
     //Reset();

  
}



Clustering::~Clustering()
{
  // default destructor
  
}



std::vector<fastjet::PseudoJet> Clustering::get_clustered_jets(std::vector<fastjet::PseudoJet> particles,enum  Clustering::E_algorithm algorithm, double jet_radius, double ptmin)
{
  if(algorithm==e_akt || algorithm==e_ca || algorithm==e_kt){
    //  IR_Saftey->add_grid(genvector,100,100); 
      fastjet::JetDefinition jet_def(fastjet::antikt_algorithm,jet_radius);
      std::vector<fastjet::PseudoJet> fatjets;
      if(algorithm==e_ca)jet_def.set_jet_algorithm(fastjet::cambridge_algorithm);
      if(algorithm==e_akt)jet_def.set_jet_algorithm(fastjet::antikt_algorithm);
      if(algorithm==e_kt)jet_def.set_jet_algorithm(fastjet::kt_algorithm);
      fastjet::ClusterSequence* clust_seq=new fastjet::ClusterSequence(particles, jet_def);
      fatjets = sorted_by_pt(clust_seq->inclusive_jets(ptmin));
      return fatjets;
  }
}

std::vector<fastjet::PseudoJet> Clustering::get_clustered_jets(std::vector<fastjet::PseudoJet> particles,enum  Clustering::E_algorithm algorithm, double ptmin, double rho, double min_r, double max_r)
{
  std::vector<fastjet::PseudoJet> fatjets;
  VariableRPlugin::ClusterType clustertype;
  if(algorithm==e_akt) clustertype=VariableRPlugin::ClusterType::AKTLIKE;
  if(algorithm==e_ca) clustertype=VariableRPlugin::ClusterType::CALIKE;
  if(algorithm==e_kt) clustertype=VariableRPlugin::ClusterType::KTLIKE;
  VariableRPlugin plugin(rho, min_r, max_r, clustertype);
  fastjet::JetDefinition jet_def(&plugin);
  fastjet::ClusterSequence clust_seq(particles, jet_def);
  fatjets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));
  return fatjets;
}  

std::vector<fastjet::PseudoJet> Clustering::get_clustered_hotvr_jets(std::vector<fastjet::PseudoJet> particles,enum  Clustering::E_algorithm algorithm, double ptmin, double rho, double min_r, double max_r, double mu, double theta, double pt_cut)
{
 HOTVR plugin_hotvr(mu, theta,min_r, max_r,rho,pt_cut, HOTVR::CALIKE);//call HOTVR algorithm
    fastjet::JetDefinition jet_def(&plugin_hotvr);
    fastjet::ClusterSequence* clust_seq=new fastjet::ClusterSequence(particles, jet_def);
    
    
    std::vector<fastjet::PseudoJet> hotvr_jets,rejected_jets,soft_jets ; //vector of hotvr_jets, jets that were rejcted durning the clustering procedure and soft jets
	
    //get vector from the plugin
    hotvr_jets=plugin_hotvr.get_jets();
    rejected_jets=plugin_hotvr.get_rejected_cluster();
    soft_jets=plugin_hotvr.get_soft_cluster();
    return hotvr_jets;
}

/*
std::vector<fastjet::PseudoJet> Clustering::get_clustered_jets(std::vector<fastjet::PseudoJet> particles,enum  Clustering::E_algorithm algorithm, double ptmin, double mu, double theta, double max_r)
{
  std::vector<fastjet::PseudoJet> fatjets;
  VariableRPlugin::ClusterType clustertype;
  if(algorithm==e_akt) clustertype=VariableRPlugin::ClusterType::AKTLIKE;
  if(algorithm==e_ca) clustertype=VariableRPlugin::ClusterType::CALIKE;
  if(algorithm==e_kt) clustertype=VariableRPlugin::ClusterType::KTLIKE;
  ClusteringVetoPlugin pluginAKT(mu, theta, max_r, clustertype);
  fastjet::JetDefinition jet_def(&plugin);
  fastjet::ClusterSequence clust_seq(particles, jet_def);
  fatjets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));
  return fatjets;
}  
*/


/*if(algo=="massjump") {
    ClusteringVetoPlugin pluginAKT(mu, theta, max_r, ClusteringVetoPlugin::AKTLIKE);
     fastjet::JetDefinition jet_defCA(&pluginAKT);
     jet_defCA2=jet_defCA;*/
