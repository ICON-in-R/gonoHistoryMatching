
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
    double uDisc ;
    std::vector<std::vector<double>> birthRate;
    std::vector<std::vector<std::vector<double>>> deathRate;
    std::vector<std::vector<std::vector<std::vector<double>>>> newborneDistribution;
    std::vector<double> transmissionRate;
    std::vector<double> latencyRate;
    std::vector<double> propAsymp;
    std::vector<double> tretmentRate;
    std::vector<double> recoveryRate;
    std::vector<std::vector<std::vector<double>>> screeningRate;
    
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
     }
    
    // this operator enables implicit Rcpp::wrap
    operator SEXP();
    
    // this ctor enables implicit Rcpp::as
    parameters(SEXP);
};
