#include <RcppCommon.h>
using namespace Rcpp;

// GonorrheaDTMModel.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <fstream>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>
#include <numeric>
#include <cmath>

using namespace std;

// #include <gonoHistoryMatching_types.h>

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

#include <Rcpp.h>


// [[Rcpp::export]]
double sumVector(std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>>>& population, int i_start, int i_end, int j_start, int j_end, int k_start, int k_end, int l_start, int l_end, int m_start, int m_end, int d_start, int d_end, int t_start, int t_end) {
  double sum = 0;
  for (int i = i_start; i <= i_end; i++) {
    for (int j = j_start; j <= j_end; j++) {
      for (int k = k_start; k <= k_end; k++) {
        for (int l = l_start; l <= l_end; l++) {
          for (int m = m_start; m <= m_end; m++) {
            for (int d = d_start; d <= d_end; d++) {
              for (int t = t_start; t <= t_end; t++)
                sum = sum + population[i][j][k][l][m][d][t];
            }
          }
        }
      }
    }
  }
  return sum;
}

// [[Rcpp::export]]
void loadInputParameters(std::string inputPath, std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>>>& Population, parameters& Parameters) {
  ifstream ifile(inputPath, ios::in);
  
  //check to see that the file was opened correctly:
  if (!ifile.is_open()) {
    cerr << "There was a problem opening the input file!\n";
    exit(1);
  }
  
  std::vector<std::vector<double>> vec;
  string lineData;
  
  //load general parameters
  getline(ifile, lineData); Parameters.cDisc = stod(lineData);
  getline(ifile, lineData); Parameters.uDisc = stod(lineData);
  
  double totalBlack = 0.0;
  double totalHispanic = 0.0;
  double totalWhite = 0.0;
  
  for (int i = 0; i < Parameters.nRaces; i++) {
    for (int j = 0; j < Parameters.nGenders; j++) {
      for (int k = 0; k < Parameters.nSexualBehs; k++) {
        for (int l = 0; l < Parameters.nSexActs; l++) {
          int m = 0;
          double data;
          vector<double> row;
          getline(ifile, lineData);
          stringstream lineStream(lineData);
          while (lineStream >> data) {
            Population[i][j][k][l][m][0][0] = data;
            for (int d = 1; d < Parameters.nDiseaseStates; d++) Population[i][j][k][l][m][d][0] = 0;
            m++;
            if (i == 0) totalBlack = totalBlack + data;
            else
              if (i == 1) totalHispanic = totalHispanic + data;
              else totalWhite = totalWhite + data;
          }
        }
      }
    }
  }
  
  //calculate newborn distribution based on initial population distribution    
  double sumB = 0.0; 
  double sumH = 0.0;
  double sumW = 0.0;
  for (int i = 0; i < Parameters.nRaces; i++) {
    Parameters.newborneDistribution.push_back(std::vector<std::vector<std::vector<double>>>());
    for (int j = 0; j < Parameters.nGenders; j++) {
      Parameters.newborneDistribution[i].push_back(std::vector<std::vector<double>>());
      for (int k = 0; k < Parameters.nSexualBehs; k++) {
        Parameters.newborneDistribution[i][j].push_back(std::vector<double>());
        for (int l = 0; l < Parameters.nSexActs; l++) {
          double pom = sumVector(Population, i, i, j, j, k, k, l, l, 0, Parameters.nAges - 1, 0, 0, 0, 0);
          if (i == 0) { Parameters.newborneDistribution[i][j][k].push_back(pom / totalBlack); sumB = Parameters.newborneDistribution[i][j][k][l] + sumB; }
          else
            if (i == 1) {
              Parameters.newborneDistribution[i][j][k].push_back(pom / totalHispanic); sumH = Parameters.newborneDistribution[i][j][k][l] + sumH;
            }
            else {
              Parameters.newborneDistribution[i][j][k].push_back(pom / totalWhite); sumW = Parameters.newborneDistribution[i][j][k][l] + sumW;
            }
        }
      }
    }
  } //rounding error correction
  Parameters.newborneDistribution[0][Parameters.nGenders - 1][Parameters.nSexualBehs - 1][Parameters.nSexActs - 1] = Parameters.newborneDistribution[0][Parameters.nGenders - 1][Parameters.nSexualBehs - 1][Parameters.nSexActs - 1] + (1 - sumB);
  Parameters.newborneDistribution[1][Parameters.nGenders - 1][Parameters.nSexualBehs - 1][Parameters.nSexActs - 1] = Parameters.newborneDistribution[1][Parameters.nGenders - 1][Parameters.nSexualBehs - 1][Parameters.nSexActs - 1] + (1 - sumH);
  Parameters.newborneDistribution[2][Parameters.nGenders - 1][Parameters.nSexualBehs - 1][Parameters.nSexActs - 1] = Parameters.newborneDistribution[2][Parameters.nGenders - 1][Parameters.nSexualBehs - 1][Parameters.nSexActs - 1] + (1 - sumW);
  
  //load yearly birth rate - per capita
  for (int i = 0; i < Parameters.nRaces; i++) {
    double d;
    std::vector<double> row;
    getline(ifile, lineData);
    stringstream lineStream(lineData);
    while (lineStream >> d) row.push_back(d);
    Parameters.birthRate.push_back(row);
  }
  
  //for (int i = 0; i < Parameters.nRaces; i++)
  //    for (int j = 0; j < 101; j++)
  //        Parameters.birthRate[i][j] = 0;
  
  //load monthly death rate 
  for (int i = 0; i < Parameters.nRaces; i++) {
    Parameters.deathRate.push_back(std::vector<std::vector<double>>());
    for (int j = 0; j < Parameters.nGenders; j++) {
      double d;
      std::vector<double> row;
      getline(ifile, lineData);
      stringstream lineStream(lineData);
      while (lineStream >> d) row.push_back(d);
      Parameters.deathRate[i].push_back(row);
    }
  }
  // for (int i = 0; i < Parameters.nRaces; i++)
  //     for (int j = 0; j < Parameters.nGenders; j++)
  //         for (int m = 0; m < Parameters.nAges; m++)
  //             Parameters.deathRate[i][j][m] = 0.0;
  
  //load transmission rate
  std::vector<double> transmissionRate;
  double d;
  getline(ifile, lineData);
  stringstream lineStream1(lineData);
  while (lineStream1 >> d) transmissionRate.push_back(d);
  
  //load latency rate
  std::vector<double> latencyRate;
  d = 0.0;
  getline(ifile, lineData);
  stringstream lineStream2(lineData);
  while (lineStream2 >> d) latencyRate.push_back(1/d);
  
  
  //load proportion of asymptomatic
  std::vector<double> propAsymp;
  d = 0.0;
  getline(ifile, lineData);
  stringstream lineStream3(lineData);
  while (lineStream3 >> d) propAsymp.push_back(d);
  
  //load treatment rate
  std::vector<double> tretmentRate;
  d = 0.0;
  getline(ifile, lineData);
  stringstream lineStream4(lineData);
  while (lineStream4 >> d) tretmentRate.push_back(1/d);
  
  //load recovery rate
  std::vector<double> recoveryRate;
  d = 0.0;
  getline(ifile, lineData);
  stringstream lineStream5(lineData);
  while (lineStream5 >> d) recoveryRate.push_back(d);
  
  //load screening rate
  std::vector<std::vector<double>> screeningRate;
  for (int i = 0; i < Parameters.nRaces; i++) {
    Parameters.screeningRate.push_back(std::vector<std::vector<double>>());
    for (int j = 0; j < Parameters.nGenders; j++) {
      double d;
      std::vector<double> row;
      getline(ifile, lineData);
      stringstream lineStream(lineData);
      while (lineStream >> d) row.push_back(d);
      Parameters.screeningRate[i].push_back(row);
    }
  }
  
  ifile.close();
}

// [[Rcpp::export]]
void modelDynamic(int year, int month, std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>>>& Population, parameters Parameters)
{
  int index = year * 12 + month;  cout << year << "\n";
  std::vector<double> newbornePopulation(Parameters.nRaces, 0.0);
  std::vector<double> X(Parameters.nDiseaseStates, 0.0);
  
  newbornePopulation[0] = Parameters.birthRate[0][year] * sumVector(Population, 0, 0, 0, Parameters.nGenders - 1, 0, Parameters.nSexualBehs - 1, 0, Parameters.nSexActs - 1, 0, Parameters.nAges - 1, 0, Parameters.nDiseaseStates - 1, index, index);
  newbornePopulation[1] = Parameters.birthRate[1][year] * sumVector(Population, 1, 1, 0, Parameters.nGenders - 1, 0, Parameters.nSexualBehs - 1, 0, Parameters.nSexActs - 1, 0, Parameters.nAges - 1, 0, Parameters.nDiseaseStates - 1, index, index);
  newbornePopulation[2] = Parameters.birthRate[2][year] * sumVector(Population, 2, 2, 0, Parameters.nGenders - 1, 0, Parameters.nSexualBehs - 1, 0, Parameters.nSexActs - 1, 0, Parameters.nAges - 1, 0, Parameters.nDiseaseStates - 1, index, index);
  
  //update indexes 0 and 1, 1->0 1-empty
  for (int i = 0; i < Parameters.nRaces; i++)
    for (int j = 0; j < Parameters.nGenders; j++)
      for (int k = 0; k < Parameters.nSexualBehs; k++)
        for (int l = 0; l < Parameters.nSexActs; l++)
          for (int m = 0; m < Parameters.nAges; m++)
            for (int d = 0; d < Parameters.nDiseaseStates; d++)
              Population[i][j][k][l][m][d][index + 1] = Population[i][j][k][l][m][d][index];
  int index1 = index + 1;
  int index2 = index + 2;
  
  int numberOfSteps = 1;
  double h = 1.0 / numberOfSteps;
  
  for (int step = 0; step < numberOfSteps; step++)
  {
    for (int i = 0; i < Parameters.nRaces; i++) {
      for (int j = 0; j < Parameters.nGenders; j++) {
        for (int k = 0; k < Parameters.nSexualBehs; k++) {
          for (int l = 0; l < Parameters.nSexActs; l++) {
            
            if (step == 0) {
              //monthly discrete dynamic
              for (int m = Parameters.nAges - 1; m > 0; m--)
              {
                for (int d = 0; d < Parameters.nDiseaseStates; d++)
                  if ((month) % 12 == 0) Population[i][j][k][l][m][d][index1] = (1 - Parameters.deathRate[i][j][m - 1]) * Population[i][j][k][l][m - 1][d][index1];
                  else  Population[i][j][k][l][m][d][index1] = (1 - Parameters.deathRate[i][j][m]) * Population[i][j][k][l][m][d][index1];
              }
              //add newbornes as susceptible
              if ((month) % 12 == 0) Population[i][j][k][l][0][0][index1] = newbornePopulation[i] * Parameters.newborneDistribution[i][j][k][l];
              else Population[i][j][k][l][0][0][index1] = (1 - Parameters.deathRate[i][j][0]) * Population[i][j][k][l][0][0][index1] + newbornePopulation[i] * Parameters.newborneDistribution[i][j][k][l];
            }
            //monthly continuous dynamic
            for (int m = 0; m < Parameters.nAges; m++)
            {
              X[0] = 0;//Susceptible
              X[1] = 0;//Exposed
              X[2] = 0;//Simptomatically infected
              X[3] = 0;//Asymptomatically infected
              X[4] = 0;//Treated
              
              X[5] = 0;//VaccinatedA with one dose
              X[6] = 0;//VaccinatedA with two doses
              X[7] = 0;
              X[8] = 0;
              X[9] = 0;
              
              Population[i][j][k][l][m][0][index2] = Population[i][j][k][l][m][0][index1] + h * X[0];
              Population[i][j][k][l][m][1][index2] = Population[i][j][k][l][m][1][index1] + h * X[1];
              Population[i][j][k][l][m][2][index2] = Population[i][j][k][l][m][2][index1] + h * X[2];
              Population[i][j][k][l][m][3][index2] = Population[i][j][k][l][m][3][index1] + h * X[3];
              Population[i][j][k][l][m][4][index2] = Population[i][j][k][l][m][4][index1] + h * X[4];
              Population[i][j][k][l][m][5][index2] = Population[i][j][k][l][m][5][index1] + h * X[5];
              Population[i][j][k][l][m][6][index2] = Population[i][j][k][l][m][6][index1] + h * X[6];
              Population[i][j][k][l][m][7][index2] = Population[i][j][k][l][m][7][index1] + h * X[7];
              Population[i][j][k][l][m][8][index2] = Population[i][j][k][l][m][8][index1] + h * X[8];
              Population[i][j][k][l][m][9][index2] = Population[i][j][k][l][m][9][index1] + h * X[8];
            }
          }}}}
          
          //update indexes 0 and 1, 1->0 1-empty
          for (int i = 0; i < Parameters.nRaces; i++)
            for (int j = 0; j < Parameters.nGenders; j++)
              for (int k = 0; k < Parameters.nSexualBehs; k++)
                for (int l = 0; l < Parameters.nSexActs; l++)
                  for (int m = 0; m < Parameters.nAges; m++)
                    for (int d = 0; d < Parameters.nDiseaseStates; d++)
                      Population[i][j][k][l][m][d][index1] = Population[i][j][k][l][m][d][index2];
  }
}

// [[Rcpp::export]]
void saveToFile(std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>>>& population, std::string filename, parameters& Parameters, int runTime, int nAges) {
  std::ofstream ofile(filename);
  if (ofile.good()) {
    ofile.flags(std::ios::fixed);
    for (int i = 0; i < Parameters.nRaces; i++)
    {
      for (int j = 0; j < Parameters.nGenders; j++)
        for (int k = 0; k < Parameters.nSexualBehs; k++)
          for (int l = 0; l < Parameters.nSexActs; l++)
            for (int m = 0; m < nAges; m++)
              for (int d = 0; d < Parameters.nDiseaseStates; d++)
              {
                for (int t = 0; t < runTime-2; t++) ofile << population[i][j][k][l][m][d][t] << "\t";
                ofile << "\n";
              }
    }
  }
  else {
    cerr << "There was a problem opening the output file!\n";
    exit(1); //exit or do additional error checking
  }
}

// [[Rcpp::export]]
int runmodel(Rcpp::List inputs)
{
  std::string outputPath = "./Outputs/DTM_Outputs.txt";
  std::string inputPath = "./Inputs/DTM_Inputs.txt";
  
  parameters Parameters;
  
  Parameters.nRaces = inputs["nRaces"];
  Parameters.nGenders = inputs["nGenders"];
  Parameters.nSexualBehs = inputs["nSexualBehs"];
  Parameters.nSexActs = inputs["nSexActs"];
  Parameters.nAges = inputs["nAges"];
  Parameters.nDiseaseStates = inputs["nDiseaseStates"];
  Parameters.timeHorizon = inputs["timeHorizon"];
  
  int maxRunTime = Parameters.timeHorizon * 12 + 2;  //time in time units
  
  std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>>> Population(Parameters.nRaces, std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>>(Parameters.nGenders, std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>(Parameters.nSexualBehs, std::vector<std::vector<std::vector<std::vector<double>>>>(Parameters.nSexActs, std::vector<std::vector<std::vector<double>>>(Parameters.nAges, std::vector<std::vector<double>>(Parameters.nDiseaseStates, std::vector<double>(maxRunTime, 0.0)))))));
  
  loadInputParameters(inputPath, Population, Parameters);
  
  for (int year = 0; year < Parameters.timeHorizon; year++) {
    for (int month = 0; month < 12; month++) {
      modelDynamic(year, month, Population, Parameters);
    }
  }
  
  //replace hardcoded value with args[]
  int ageGroupSize = 5;
  
  //aggregate to age groups
  for (int i = 0; i < Parameters.nRaces; i++){
    for (int j = 0; j < Parameters.nGenders; j++){
      for (int k = 0; k < Parameters.nSexualBehs; k++){
        for (int l = 0; l < Parameters.nSexActs; l++){
          for (int d = 0; d < Parameters.nDiseaseStates; d++){
            for (int t = 0; t < maxRunTime; t++)
            {
              int index = 0;
              for (int m1 = 0; m1 < Parameters.nAges / ageGroupSize; m1++) {
                double sum = 0;
                for (int m = 0; m < ageGroupSize; m++) {
                  sum = sum + Population[i][j][k][l][index * ageGroupSize + m][d][t];
                  Population[i][j][k][l][index][d][t] = sum;
                  index++;
                }
              }
              //change this in order to agregate all age groups that are not included - currently works only for 100 years old
              Population[i][j][k][l][Parameters.nAges / ageGroupSize][d][t] = Population[i][j][k][l][Parameters.nAges-1][d][t];
            }
            }}}}}
            
            //save DTM outputs to file  
            saveToFile(Population, outputPath, Parameters, maxRunTime, Parameters.nAges / ageGroupSize + 1);
  
  return 0;
}


