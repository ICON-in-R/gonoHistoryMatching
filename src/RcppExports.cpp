// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "gonoHistoryMatching_types.h"
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// loadDemographics
void loadDemographics(std::string inputPath, parameters& Parameters);
RcppExport SEXP _gonoHistoryMatching_loadDemographics(SEXP inputPathSEXP, SEXP ParametersSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type inputPath(inputPathSEXP);
    Rcpp::traits::input_parameter< parameters& >::type Parameters(ParametersSEXP);
    loadDemographics(inputPath, Parameters);
    return R_NilValue;
END_RCPP
}
// loadParameters
void loadParameters(std::string inputPath, parameters& Parameters, psa_parameters* psaParameters, int NumbSim);
RcppExport SEXP _gonoHistoryMatching_loadParameters(SEXP inputPathSEXP, SEXP ParametersSEXP, SEXP psaParametersSEXP, SEXP NumbSimSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type inputPath(inputPathSEXP);
    Rcpp::traits::input_parameter< parameters& >::type Parameters(ParametersSEXP);
    Rcpp::traits::input_parameter< psa_parameters* >::type psaParameters(psaParametersSEXP);
    Rcpp::traits::input_parameter< int >::type NumbSim(NumbSimSEXP);
    loadParameters(inputPath, Parameters, psaParameters, NumbSim);
    return R_NilValue;
END_RCPP
}
// loadCalibratioTargets
void loadCalibratioTargets(std::string inputPath, parameters& Parameters);
RcppExport SEXP _gonoHistoryMatching_loadCalibratioTargets(SEXP inputPathSEXP, SEXP ParametersSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type inputPath(inputPathSEXP);
    Rcpp::traits::input_parameter< parameters& >::type Parameters(ParametersSEXP);
    loadCalibratioTargets(inputPath, Parameters);
    return R_NilValue;
END_RCPP
}
// loadCalibrationParameters
void loadCalibrationParameters(std::string inputPath, parameters& Parameters, int numbParallel);
RcppExport SEXP _gonoHistoryMatching_loadCalibrationParameters(SEXP inputPathSEXP, SEXP ParametersSEXP, SEXP numbParallelSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type inputPath(inputPathSEXP);
    Rcpp::traits::input_parameter< parameters& >::type Parameters(ParametersSEXP);
    Rcpp::traits::input_parameter< int >::type numbParallel(numbParallelSEXP);
    loadCalibrationParameters(inputPath, Parameters, numbParallel);
    return R_NilValue;
END_RCPP
}
// updatedCalibrationParameters
void updatedCalibrationParameters(std::string filename, parameters& Parameters);
RcppExport SEXP _gonoHistoryMatching_updatedCalibrationParameters(SEXP filenameSEXP, SEXP ParametersSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< parameters& >::type Parameters(ParametersSEXP);
    updatedCalibrationParameters(filename, Parameters);
    return R_NilValue;
END_RCPP
}
// runmodel
int runmodel(std::string calibrationPath);
RcppExport SEXP _gonoHistoryMatching_runmodel(SEXP calibrationPathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type calibrationPath(calibrationPathSEXP);
    rcpp_result_gen = Rcpp::wrap(runmodel(calibrationPath));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gonoHistoryMatching_loadDemographics", (DL_FUNC) &_gonoHistoryMatching_loadDemographics, 2},
    {"_gonoHistoryMatching_loadParameters", (DL_FUNC) &_gonoHistoryMatching_loadParameters, 4},
    {"_gonoHistoryMatching_loadCalibratioTargets", (DL_FUNC) &_gonoHistoryMatching_loadCalibratioTargets, 2},
    {"_gonoHistoryMatching_loadCalibrationParameters", (DL_FUNC) &_gonoHistoryMatching_loadCalibrationParameters, 3},
    {"_gonoHistoryMatching_updatedCalibrationParameters", (DL_FUNC) &_gonoHistoryMatching_updatedCalibrationParameters, 2},
    {"_gonoHistoryMatching_runmodel", (DL_FUNC) &_gonoHistoryMatching_runmodel, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_gonoHistoryMatching(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
