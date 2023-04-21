library(Rcpp)

Rcpp::compileAttributes()
sourceCpp("src/GonorrheaDTM.cpp", windowsDebugDLL = FALSE)

inputs <- list(
  nRaces = 3,
  nGenders = 2,
  nSexualBehs = 3,
  nSexActs = 3,
  nAges = 101,
  nDiseaseStates = 10,
  timeHorizon = 100,
  ageGroupSize = 5,
  cDisc = 0.0,
  uDisc = 0.0)

# birthRate =
# deathRate =
# newborneDistribution =
# transmissionRate =
# latencyRate =
# propAsymp =
# tretmentRate =
# recoveryRate =
# screeningRate =

res <- runmodel(inputs)
