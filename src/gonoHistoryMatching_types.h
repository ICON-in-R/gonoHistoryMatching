// gonoHistoryMatching types


#ifndef GONOHISTORYMATCHING_H
#define GONOHISTORYMATCHING_H

#include <Rcpp.h>

using namespace Rcpp;

class parameters {
  public:
    int nRaces;
    int nGenders;
    int nSexualBehs;
    int nSexActs;
    int nAges;
    int nDiseaseStates;
    int timeHorizon;
    double cDisc;
    double uDisc;
    std::vector<std::vector<std::vector<std::vector<double>>>> popDistribution;
    std::vector<std::vector<double>> birthRate;
    std::vector<std::vector<std::vector<double>>> deathRate;
    std::vector<std::vector<std::vector<std::vector<double>>>> newborneDistribution;
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> transmissionRate;
    std::vector<double> latencyRate;
    std::vector<double> propAsymp;
    std::vector<double> tretmentRate;
    std::vector<double> recoveryRate;
    std::vector<std::vector<std::vector<double>>> screeningRate;
    std::vector<std::vector<double>> pRaces;
    std::vector<std::vector<double>> pSexActs;
    std::vector<std::vector<double>> pAges;
    std::vector<double> nmbPartners;
    std::vector<std::vector<std::vector<std::vector<double>>>> sexFrequency;
    
    int Age_min;
    int Age_max;
    int ageGroup_size;
    std::vector<double> ageGroup_bound;
    
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> calibrationTarget;
    
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> vaccA1;
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> vaccA2;
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> vaccA3;
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> vaccA4;
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> vaccB2;
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> vaccB3;
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> vaccB4;
    std::vector<double> waningImmunityA1;
    std::vector<double> waningImmunityA2;
    std::vector<double> waningImmunityA3;
    std::vector<double> waningImmunityA4;
    std::vector<double> waningImmunityB2;
    std::vector<double> waningImmunityB3;
    std::vector<double> waningImmunityB4;
    
    parameters() {
      nRaces = 0;
      nGenders = 0;
      nSexualBehs = 0;
      nSexActs = 0;
      nAges = 0;
      nDiseaseStates = 0;
      timeHorizon = 0;
      cDisc = 0.0;
      uDisc = 0.0;
      Age_min = 0;
      Age_max = 0;
      ageGroup_size = 4;
    }
};

class psa_parameters {
  public:
    std::vector<std::vector<double>> latencyRate;
    std::vector<std::vector<double>> propAsymp;
    std::vector<std::vector<double>> tretmentRate;
    std::vector<std::vector<double>> recoveryRate;
    std::vector<std::vector<double>> recoveryAsympRate;
    std::vector<std::vector<std::vector<std::vector<double>>>> screeningRate;
    
    psa_parameters() {
    }
};

// https://gallery.rcpp.org/articles/mixing-modules-and-export/
RCPP_EXPOSED_AS(parameters)
RCPP_EXPOSED_AS(psa_parameters)
  
#endif