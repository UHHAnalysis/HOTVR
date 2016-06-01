#include <iostream>
#include <fstream>
#include "include/Clustering.h"
#include "include/Clusteringinfo.hh"
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

Clustering::Clustering(): m_logger("Clustering")
{
}

Clustering::Clustering(std::string clustering) : m_logger("Clustering")
{
  // constructor: initialise all variables
  ifstream myfile(clustering);
  std::string line;
  if (myfile.is_open()) {
    while ( myfile>>line )
    {
      std::string value;
      myfile>>value;
      if(line=="rho") _rho=get_value(value);
      if(line=="mu") _mu=get_value(value);
      if(line=="theta") _theta=get_value(value);
      if(line=="rmin") _rmin=get_value(value);
      if(line=="rmax") _rmax=get_value(value);
      if(line=="ptmin") _ptmin=get_value(value);
      if(line=="clustering") _clustering_algorithmus=value;
      if(line=="modus") set_modus(value);
      if(line=="radius") _radius=get_value(value);
      if(line=="pt_cut") _pt_cut=get_value(value);
    }
    myfile.close();
  }
  else  std::cout<<"file not found"<<std::endl;
  show_settings();
  _setup_jet.set_user_info(new Clusteringinfo(_radius,_rho,_mu,_theta,_rmin,_rmax,_ptmin,_pt_cut,_clustering_algorithmus));
  
 
     //Reset();

  
}

void Clustering::show_settings(){
  std::cout<<"-----Clustering------"<<std::endl;
  std::cout<<"Algorithmus: "<<_clustering_algorithmus<<std::endl;
  std::cout<<"Modus: "<<get_modus()<<std::endl;
  std::cout<<"Minimum Jet pT: "<<_ptmin<<"GeV"<<std::endl;
  if(_clustering_algorithmus=="squential") std::cout<<"Clustering radius "<<_radius<<std::endl;
  if(_clustering_algorithmus=="hotvr" || _clustering_algorithmus=="variableR") std::cout<<"rho: "<<_rho<<"GeV"<<std::endl;
  if(_clustering_algorithmus=="hotvr" || _clustering_algorithmus=="variableR") std::cout<<"minimum radius: "<<_rmin<<std::endl;
  if(_clustering_algorithmus=="hotvr" || _clustering_algorithmus=="variableR") std::cout<<"maximum radius: "<<_rmax<<std::endl;
  if(_clustering_algorithmus=="hotvr") std::cout<<"mu: "<<_mu<<"GeV"<<std::endl;
  if(_clustering_algorithmus=="hotvr") std::cout<<"theta: "<<_theta<<std::endl;
  std::cout<<"---------------------"<<std::endl;
  
  
}



Clustering::~Clustering()
{
  // default destructor
  // if(_fatjets.size()!=0) delete _fatjets.at(0).associated_cluster_sequence();
}

double Clustering::get_value(std::string word){
  return atof(word.c_str());
}

std::string Clustering::get_modus(){
  if(m_test==e_ca) return "Cambridge/Aachen";
  if(m_test==e_akt) return "Anti-kt";
  if(m_test==e_kt) return "kt";
  else return "NOT SET!";
}


void Clustering::set_modus(std::string word){
  if(word=="cambridge") m_test=e_ca;
  if(word=="antikt") m_test=e_akt;
  if(word=="kt") m_test=e_kt;
}


void Clustering::jets(std::string clustering){
  
}

std::vector<fastjet::PseudoJet> Clustering::get_clustered_jets(std::vector<fastjet::PseudoJet> particles){
  if(_clustering_algorithmus=="hotvr") return get_clustered_hotvr_jets(particles,m_test, _ptmin, _rho, _rmin, _rmax, _mu, _theta, _pt_cut);
  if(_clustering_algorithmus=="sequential") return get_clustered_jets(particles,m_test ,_radius,_ptmin);
  if(_clustering_algorithmus=="variableR") return get_clustered_jets(particles,m_test, _ptmin, _rho, _rmin, _rmax);
  else {
    std::cout<<"Clustering algorithmus not found"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"Please choose: "<<std::endl;
    std::cout<<"squential"<<std::endl;
    std::cout<<"variableR"<<std::endl;
    std::cout<<"hotvr"<<std::endl;
    exit(0);
  }
}


std::vector<fastjet::PseudoJet> Clustering::set_properties(std::vector<fastjet::PseudoJet> jets){
   fastjet::PseudoJet modified_jet;
  modified_jet.set_user_info(new Clusteringinfo(_radius,_rho,_mu,_theta,_rmin,_rmax,_ptmin,_pt_cut,_clustering_algorithmus));
  jets.push_back(modified_jet);
  return jets;
}


std::vector<fastjet::PseudoJet> Clustering::get_clustered_jets(std::vector<fastjet::PseudoJet> particles,enum  Clustering::E_algorithm algorithm, double jet_radius, double ptmin)
{
  if(algorithm==e_akt || algorithm==e_ca || algorithm==e_kt){
  
    fastjet::JetDefinition *jetdef;
    std::vector<fastjet::PseudoJet> fatjets;
  
      if(algorithm==e_ca) jetdef= new fastjet::JetDefinition(fastjet::cambridge_algorithm,jet_radius);
      if(algorithm==e_akt)jetdef= new fastjet::JetDefinition(fastjet::antikt_algorithm,jet_radius);
      if(algorithm==e_kt)jetdef= new fastjet::JetDefinition(fastjet::cambridge_algorithm,jet_radius);
     
      _clust_seq2=new fastjet::ClusterSequence(particles, *jetdef);
    
      fatjets = sorted_by_pt(_clust_seq2->inclusive_jets(ptmin));
      _fatjets=fatjets;

      delete jetdef;
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
   HOTVR::ClusterType clustertype;
  if(algorithm==e_akt) clustertype=HOTVR::ClusterType::AKTLIKE;
  if(algorithm==e_ca) clustertype=HOTVR::ClusterType::CALIKE;
  if(algorithm==e_kt) clustertype=HOTVR::ClusterType::KTLIKE;
 HOTVR plugin_hotvr(mu, theta,min_r, max_r,rho,pt_cut, clustertype);//call HOTVR algorithm
    fastjet::JetDefinition jet_def(&plugin_hotvr);
    _clust_seq2=new fastjet::ClusterSequence(particles, jet_def);
    
    
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

void Clustering::Reset()
{
  delete _clust_seq2;
  // delete _fatjets;
  
}
