// gonoHistoryMatching types


#ifndef GONOHISTORYMATCHING_JG_H
#define GONOHISTORYMATCHING_JG_H

#include <Rcpp.h>

using namespace Rcpp;

// JG
#define MODEL_FLOAT float

class parameters {
  public:
    int nRaces;
    int nGenders;
    int nSexualBehs;
    int nSexActs;
    int nAges;
    int nDiseaseStates;
    int timeHorizon;
    MODEL_FLOAT cDisc;
    MODEL_FLOAT uDisc;
    std::vector<std::vector<std::vector<std::vector<MODEL_FLOAT>>>> popDistribution;
    std::vector<std::vector<MODEL_FLOAT>> birthRate;
    std::vector<std::vector<std::vector<MODEL_FLOAT>>> deathRate;
    std::vector<std::vector<std::vector<std::vector<MODEL_FLOAT>>>> newborneDistribution;
    std::vector<std::vector<std::vector<std::vector<std::vector<MODEL_FLOAT>>>>> transmissionRate;
    std::vector<MODEL_FLOAT> latencyRate;
    std::vector<MODEL_FLOAT> propAsymp;
    std::vector<MODEL_FLOAT> tretmentRate;
    std::vector<MODEL_FLOAT> recoveryRate;
    std::vector<std::vector<std::vector<MODEL_FLOAT>>> screeningRate;
    std::vector<std::vector<MODEL_FLOAT>> pRaces;
    std::vector<std::vector<MODEL_FLOAT>> pSexActs;
    std::vector<std::vector<MODEL_FLOAT>> pAges;
    std::vector<MODEL_FLOAT> nmbPartners;
    std::vector<std::vector<std::vector<std::vector<MODEL_FLOAT>>>> sexFrequency;
    
    int Age_min;
    int Age_max;
    int ageGroup_size;
    std::vector<MODEL_FLOAT> ageGroup_bound;
    
    std::vector<std::vector<std::vector<std::vector<std::vector<MODEL_FLOAT>>>>> calibrationTarget;
    
    std::vector<std::vector<std::vector<std::vector<std::vector<MODEL_FLOAT>>>>> vaccA1;
    std::vector<std::vector<std::vector<std::vector<std::vector<MODEL_FLOAT>>>>> vaccA2;
    std::vector<std::vector<std::vector<std::vector<std::vector<MODEL_FLOAT>>>>> vaccA3;
    std::vector<std::vector<std::vector<std::vector<std::vector<MODEL_FLOAT>>>>> vaccA4;
    std::vector<std::vector<std::vector<std::vector<std::vector<MODEL_FLOAT>>>>> vaccB2;
    std::vector<std::vector<std::vector<std::vector<std::vector<MODEL_FLOAT>>>>> vaccB3;
    std::vector<std::vector<std::vector<std::vector<std::vector<MODEL_FLOAT>>>>> vaccB4;
    std::vector<MODEL_FLOAT> waningImmunityA1;
    std::vector<MODEL_FLOAT> waningImmunityA2;
    std::vector<MODEL_FLOAT> waningImmunityA3;
    std::vector<MODEL_FLOAT> waningImmunityA4;
    std::vector<MODEL_FLOAT> waningImmunityB2;
    std::vector<MODEL_FLOAT> waningImmunityB3;
    std::vector<MODEL_FLOAT> waningImmunityB4;
    
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
    std::vector<std::vector<MODEL_FLOAT>> latencyRate;
    std::vector<std::vector<MODEL_FLOAT>> propAsymp;
    std::vector<std::vector<MODEL_FLOAT>> tretmentRate;
    std::vector<std::vector<MODEL_FLOAT>> recoveryRate;
    std::vector<std::vector<MODEL_FLOAT>> recoveryAsympRate;
    std::vector<std::vector<std::vector<std::vector<MODEL_FLOAT>>>> screeningRate;
    //fixed_vector<fixed_vector<fixed_vector<fixed_vector<MODEL_FLOAT, 30>, 30>, 30>, 30> screeningRate2;

    psa_parameters() {
    }
};

// https://gallery.rcpp.org/articles/mixing-modules-and-export/
RCPP_EXPOSED_AS(parameters)
RCPP_EXPOSED_AS(psa_parameters)
  
#endif