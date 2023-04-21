
class parameters {
  public:
    parameters();
  
  int nRaces = 0;
  int nGenders = 0;
  int nSexualBehs = 0;
  int nSexActs = 0;
  int nAges = 0;
  int nDiseaseStates = 0;
  int timeHorizon = 0;
  double cDisc = 0.0;
  double uDisc = 0.0;
  
  std::vector<std::vector<double>> birthRate;
  std::vector<std::vector<std::vector<double>>> deathRate;
  std::vector<std::vector<std::vector<std::vector<double>>>> newborneDistribution;
  std::vector<double> transmissionRate;
  std::vector<double> latencyRate;
  std::vector<double> propAsymp;
  std::vector<double> tretmentRate;
  std::vector<double> recoveryRate;
  std::vector<std::vector<std::vector<double>>> screeningRate;
  
  // this operator enables implicit Rcpp::wrap
  operator SEXP();
  
  // this ctor enables implicit Rcpp::as
  parameters(SEXP);
};
