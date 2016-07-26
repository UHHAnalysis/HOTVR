#ifndef CLUSTERINGINFO_HH_
#define CLUSTERINGINFO_HH_
#include <string>
#include "fastjet/PseudoJet.hh"

//FASTJET_BEGIN_NAMESPACE  

//namespace contrib{

class Clusteringinfo : public fastjet::PseudoJet::UserInfoBase {
  public:
    double cluster_radius() const{ return _radius;};
    double rho() const{return _rho;};
    double mu() const{return _mu;};
    double theta() const{return _theta;};
    double rmin() const{return _rmin;};
    double rmax() const{return _rmax;};
    double ptmin() const{return _ptmin;};
    double pt_cut() const{return _pt_cut;};
    std::string clustering() const{return _clustering_algorithmus;};
   
    //  HOTVRinfo(fastjet::PseudoJet* jet):_parent(jet){};
    Clusteringinfo(double radius, double rho, double mu, double theta, double rmin, double rmax, double ptmin,double pt_cut, std::string clustering): _rho(rho), _mu(mu), _theta(theta), _rmin(rmin), _rmax(rmax), _ptmin(ptmin),_radius(radius),_pt_cut(pt_cut), _clustering_algorithmus(clustering){};
    ~Clusteringinfo(){};
    
  private:
    double _rho, _mu, _theta, _rmin, _rmax, _ptmin,_radius,_pt_cut;
    std::string _clustering_algorithmus;
  };
  

//}
//FASTJET_END_NAMESPACE
#endif // __FASTJET_CONTRIB_SHAPE_WITH_COMPONENTS_HH__
