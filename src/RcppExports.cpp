// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "gonoHistoryMatching_types.h"
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// sumVector
double sumVector(std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>>>& population, int i_start, int i_end, int j_start, int j_end, int k_start, int k_end, int l_start, int l_end, int m_start, int m_end, int d_start, int d_end, int t_start, int t_end);
RcppExport SEXP _gonoHistoryMatching_sumVector(SEXP populationSEXP, SEXP i_startSEXP, SEXP i_endSEXP, SEXP j_startSEXP, SEXP j_endSEXP, SEXP k_startSEXP, SEXP k_endSEXP, SEXP l_startSEXP, SEXP l_endSEXP, SEXP m_startSEXP, SEXP m_endSEXP, SEXP d_startSEXP, SEXP d_endSEXP, SEXP t_startSEXP, SEXP t_endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>>>& >::type population(populationSEXP);
    Rcpp::traits::input_parameter< int >::type i_start(i_startSEXP);
    Rcpp::traits::input_parameter< int >::type i_end(i_endSEXP);
    Rcpp::traits::input_parameter< int >::type j_start(j_startSEXP);
    Rcpp::traits::input_parameter< int >::type j_end(j_endSEXP);
    Rcpp::traits::input_parameter< int >::type k_start(k_startSEXP);
    Rcpp::traits::input_parameter< int >::type k_end(k_endSEXP);
    Rcpp::traits::input_parameter< int >::type l_start(l_startSEXP);
    Rcpp::traits::input_parameter< int >::type l_end(l_endSEXP);
    Rcpp::traits::input_parameter< int >::type m_start(m_startSEXP);
    Rcpp::traits::input_parameter< int >::type m_end(m_endSEXP);
    Rcpp::traits::input_parameter< int >::type d_start(d_startSEXP);
    Rcpp::traits::input_parameter< int >::type d_end(d_endSEXP);
    Rcpp::traits::input_parameter< int >::type t_start(t_startSEXP);
    Rcpp::traits::input_parameter< int >::type t_end(t_endSEXP);
    rcpp_result_gen = Rcpp::wrap(sumVector(population, i_start, i_end, j_start, j_end, k_start, k_end, l_start, l_end, m_start, m_end, d_start, d_end, t_start, t_end));
    return rcpp_result_gen;
END_RCPP
}
// loadInputParameters
void loadInputParameters(std::string inputPath, std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>>>& Population, parameters& Parameters);
RcppExport SEXP _gonoHistoryMatching_loadInputParameters(SEXP inputPathSEXP, SEXP PopulationSEXP, SEXP ParametersSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type inputPath(inputPathSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>>>& >::type Population(PopulationSEXP);
    Rcpp::traits::input_parameter< parameters& >::type Parameters(ParametersSEXP);
    loadInputParameters(inputPath, Population, Parameters);
    return R_NilValue;
END_RCPP
}
// modelDynamic
void modelDynamic(int year, int month, std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>>>& Population, parameters& Parameters);
RcppExport SEXP _gonoHistoryMatching_modelDynamic(SEXP yearSEXP, SEXP monthSEXP, SEXP PopulationSEXP, SEXP ParametersSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type year(yearSEXP);
    Rcpp::traits::input_parameter< int >::type month(monthSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>>>& >::type Population(PopulationSEXP);
    Rcpp::traits::input_parameter< parameters& >::type Parameters(ParametersSEXP);
    modelDynamic(year, month, Population, Parameters);
    return R_NilValue;
END_RCPP
}
// saveToFile
void saveToFile(std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>>>& population, std::string filename, parameters& Parameters, int runTime, int nAges);
RcppExport SEXP _gonoHistoryMatching_saveToFile(SEXP populationSEXP, SEXP filenameSEXP, SEXP ParametersSEXP, SEXP runTimeSEXP, SEXP nAgesSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>>>& >::type population(populationSEXP);
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< parameters& >::type Parameters(ParametersSEXP);
    Rcpp::traits::input_parameter< int >::type runTime(runTimeSEXP);
    Rcpp::traits::input_parameter< int >::type nAges(nAgesSEXP);
    saveToFile(population, filename, Parameters, runTime, nAges);
    return R_NilValue;
END_RCPP
}
// runmodel
int runmodel(Rcpp::List inputs);
RcppExport SEXP _gonoHistoryMatching_runmodel(SEXP inputsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type inputs(inputsSEXP);
    rcpp_result_gen = Rcpp::wrap(runmodel(inputs));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gonoHistoryMatching_sumVector", (DL_FUNC) &_gonoHistoryMatching_sumVector, 15},
    {"_gonoHistoryMatching_loadInputParameters", (DL_FUNC) &_gonoHistoryMatching_loadInputParameters, 3},
    {"_gonoHistoryMatching_modelDynamic", (DL_FUNC) &_gonoHistoryMatching_modelDynamic, 4},
    {"_gonoHistoryMatching_saveToFile", (DL_FUNC) &_gonoHistoryMatching_saveToFile, 5},
    {"_gonoHistoryMatching_runmodel", (DL_FUNC) &_gonoHistoryMatching_runmodel, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_gonoHistoryMatching(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
