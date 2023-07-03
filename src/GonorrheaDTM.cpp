#include <RcppCommon.h>
using namespace Rcpp;

// GonorrheaDTMModel.cpp : This file contains the 'main' function
// Program execution begins and ends there

#include <fstream>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>
#include <numeric>
#include <cmath>
#include <omp.h>
#include <thread>

using namespace std;

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

#include <math.h>

#include "gonoHistoryMatching_types.h"

#include <Rcpp.h>

#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4

#define NR_END 1
#define FREE_ARG char*


// [[Rcpp::export]]
double SIGN(double a, double b)
{
    if (b >= 0) a = abs(a);
    else a = -abs(a);
    return a;
}

// [[Rcpp::export]]
void nrerror(std::string error_text)
/* Numerical Recipes standard error handler */
{
    //fprintf(stderr, "Numerical Recipes run-time error...\n");
    //fprintf(stderr, "%s\n", error_text);
    //fprintf(stderr, "...now exiting to system...\n");
    exit(1);
}

// [[Rcpp::export]]
double* vector1(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
    double* v;

    v = (double*)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(double)));
    if (!v) nrerror("allocation failure in vector()");
    return v - nl + NR_END;
}

// [[Rcpp::export]]
int* ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
    int* v;

    v = (int*)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(int)));
    if (!v) nrerror("allocation failure in ivector()");
    return v - nl + NR_END;
}

// [[Rcpp::export]]
unsigned char* cvector(long nl, long nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
    unsigned char* v;

    v = (unsigned char*)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(unsigned char)));
    if (!v) nrerror("allocation failure in cvector()");
    return v - nl + NR_END;
}

// [[Rcpp::export]]
unsigned long* lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
    unsigned long* v;

    v = (unsigned long*)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(long)));
    if (!v) nrerror("allocation failure in lvector()");
    return v - nl + NR_END;
}

// [[Rcpp::export]]
double* dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
    double* v;

    v = (double*)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(double)));
    if (!v) nrerror("allocation failure in dvector()");
    return v - nl + NR_END;
}

// [[Rcpp::export]]
double** matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    double** m;

    /* allocate pointers to rows */
    m = (double**)malloc((size_t)((nrow + NR_END) * sizeof(double*)));
    if (!m) nrerror("allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl] = (double*)malloc((size_t)((nrow * ncol + NR_END) * sizeof(double)));
    if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;

    /* return pointer to array of pointers to rows */
    return m;
}

// [[Rcpp::export]]
double** dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    double** m;

    /* allocate pointers to rows */
    m = (double**)malloc((size_t)((nrow + NR_END) * sizeof(double*)));
    if (!m) nrerror("allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl] = (double*)malloc((size_t)((nrow * ncol + NR_END) * sizeof(double)));
    if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;

    /* return pointer to array of pointers to rows */
    return m;
}

// [[Rcpp::export]]
int** imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    int** m;

    /* allocate pointers to rows */
    m = (int**)malloc((size_t)((nrow + NR_END) * sizeof(int*)));
    if (!m) nrerror("allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;


    /* allocate rows and set pointers to them */
    m[nrl] = (int*)malloc((size_t)((nrow * ncol + NR_END) * sizeof(int)));
    if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;

    /* return pointer to array of pointers to rows */
    return m;
}

// [[Rcpp::export]]
double** submatrix(double** a, long oldrl, long oldrh, long oldcl, long oldch, long newrl, long newcl)
{
// point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch]
    long i, j, nrow = oldrh - oldrl + 1, ncol = oldcl - newcl;
    double** m;

    // allocate array of pointers to rows
    m = (double**)malloc((size_t)((nrow + NR_END) * sizeof(double*)));
    if (!m) nrerror("allocation failure in submatrix()");
    m += NR_END;
    m -= newrl;

    // set pointers to rows
    for (i = oldrl, j = newrl; i <= oldrh; i++, j++) m[j] = a[i] + ncol;

    // return pointer to array of pointers to rows
    return m;
}

// [[Rcpp::export]]
double** convert_matrix(double* a, long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
    long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    double** m;

    // allocate pointers to rows
    m = (double**)malloc((size_t)((nrow + NR_END) * sizeof(double*)));
    if (!m) nrerror("allocation failure in convert_matrix()");
    m += NR_END;
    m -= nrl;

    // set pointers to rows
    m[nrl] = a - ncl;
    for (i = 1, j = nrl + 1; i < nrow; i++, j++) m[j] = m[j - 1] + ncol;
    // return pointer to array of pointers to rows
    return m;
}

// [[Rcpp::export]]
double*** f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
// allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh]
    long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
    double*** t;

    // allocate pointers to pointers to rows
    t = (double***)malloc((size_t)((nrow + NR_END) * sizeof(double**)));
    if (!t) nrerror("allocation failure 1 in f3tensor()");
    t += NR_END;
    t -= nrl;

    // allocate pointers to rows and set pointers to them
    t[nrl] = (double**)malloc((size_t)((nrow * ncol + NR_END) * sizeof(double*)));
    if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
    t[nrl] += NR_END;
    t[nrl] -= ncl;

    // allocate rows and set pointers to them
    t[nrl][ncl] = (double*)malloc((size_t)((nrow * ncol * ndep + NR_END) * sizeof(double)));
    if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
    t[nrl][ncl] += NR_END;
    t[nrl][ncl] -= ndl;

    for (j = ncl + 1; j <= nch; j++) t[nrl][j] = t[nrl][j - 1] + ndep;
    for (i = nrl + 1; i <= nrh; i++) {
        t[i] = t[i - 1] + ncol;
        t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
        for (j = ncl + 1; j <= nch; j++) t[i][j] = t[i][j - 1] + ndep;
    }

    // return pointer to array of pointers to rows
    return t;
}

// [[Rcpp::export]]
void free_vector(double* v, long nl, long nh)
{
// free a double vector allocated with vector()
    free((FREE_ARG)(v + nl - NR_END));
}

// [[Rcpp::export]]
void free_ivector(int* v, long nl, long nh)
{
// free an int vector allocated with ivector()
    free((FREE_ARG)(v + nl - NR_END));
}

// [[Rcpp::export]]
void free_cvector(unsigned char* v, long nl, long nh)
{
// free an unsigned char vector allocated with cvector()
    free((FREE_ARG)(v + nl - NR_END));
}

// [[Rcpp::export]]
void free_lvector(unsigned long* v, long nl, long nh)
{
// free an unsigned long vector allocated with lvector()
    free((FREE_ARG)(v + nl - NR_END));
}

// [[Rcpp::export]]
void free_dvector(double* v, long nl, long nh)
{
// free a double vector allocated with dvector()
    free((FREE_ARG)(v + nl - NR_END));
}

// [[Rcpp::export]]
void free_matrix(double** m, long nrl, long nrh, long ncl, long nch)
{
// free a double matrix allocated by matrix()
    free((FREE_ARG)(m[nrl] + ncl - NR_END));
    free((FREE_ARG)(m + nrl - NR_END));
}

// [[Rcpp::export]]
void free_dmatrix(double** m, long nrl, long nrh, long ncl, long nch)
{
// free a double matrix allocated by dmatrix()
    free((FREE_ARG)(m[nrl] + ncl - NR_END));
    free((FREE_ARG)(m + nrl - NR_END));
}

// [[Rcpp::export]]
void free_imatrix(int** m, long nrl, long nrh, long ncl, long nch)
{
// free an int matrix allocated by imatrix()
    free((FREE_ARG)(m[nrl] + ncl - NR_END));
    free((FREE_ARG)(m + nrl - NR_END));
}

// [[Rcpp::export]]
void free_submatrix(double** b, long nrl, long nrh, long ncl, long nch)
{
// free a submatrix allocated by submatrix()
    free((FREE_ARG)(b + nrl - NR_END));
}

// [[Rcpp::export]]
void free_convert_matrix(double** b, long nrl, long nrh, long ncl, long nch)
{
// free a matrix allocated by convert_matrix()
    free((FREE_ARG)(b + nrl - NR_END));
}

// [[Rcpp::export]]
void free_f3tensor(double*** t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
// free a double f3tensor allocated by f3tensor()
    free((FREE_ARG)(t[nrl][ncl] + ndl - NR_END));
    free((FREE_ARG)(t[nrl] + ncl - NR_END));
    free((FREE_ARG)(t + nrl - NR_END));
}


// [[Rcpp::export]]
void rkck(double** PopulationX, parameters& Parameters, psa_parameters* psaParameters, int race, int gender, int sexBehs, int age,
          double y[], double dydx[], int n, double x, double h, double yout[], double yerr[], double year, double month,
          void (*derivs)(double**, parameters&, psa_parameters*, int, int, int, int, double, double[], double[], double, double, int, int), int ns, int np)
{
    int i;
    static double a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875, b21 = 0.2,
        b31 = 3.0 / 40.0, b32 = 9.0 / 40.0, b41 = 0.3, b42 = -0.9, b43 = 1.2,
        b51 = -11.0 / 54.0, b52 = 2.5, b53 = -70.0 / 27.0, b54 = 35.0 / 27.0,
        b61 = 1631.0 / 55296.0, b62 = 175.0 / 512.0, b63 = 575.0 / 13824.0,
        b64 = 44275.0 / 110592.0, b65 = 253.0 / 4096.0, c1 = 37.0 / 378.0,
        c3 = 250.0 / 621.0, c4 = 125.0 / 594.0, c6 = 512.0 / 1771.0,
        dc5 = -277.00 / 14336.0;
    double dc1 = c1 - 2825.0 / 27648.0, dc3 = c3 - 18575.0 / 48384.0,
        dc4 = c4 - 13525.0 / 55296.0, dc6 = c6 - 0.25;
    double* ak2, * ak3, * ak4, * ak5, * ak6, * ytemp;
    ak2 = vector1(1, n);
    ak3 = vector1(1, n);
    ak4 = vector1(1, n);
    ak5 = vector1(1, n);
    ak6 = vector1(1, n);
    ytemp = vector1(1, n);

    for (i = 0; i < n; i++) {
        ytemp[i] = y[i] + b21 * h * dydx[i];
    }
    (*derivs)(PopulationX, Parameters, psaParameters, race, gender, sexBehs, age, x + a2 * h, ytemp, ak2, year, month, ns, np);

    for (i = 0; i < n; i++)
        ytemp[i] = y[i] + h * (b31 * dydx[i] + b32 * ak2[i]);
    (*derivs)(PopulationX, Parameters, psaParameters, race, gender, sexBehs, age, x + a3 * h, ytemp, ak3, year, month, ns, np);

    for (i = 0; i < n; i++)
        ytemp[i] = y[i] + h * (b41 * dydx[i] + b42 * ak2[i] + b43 * ak3[i]);
    (*derivs)(PopulationX, Parameters, psaParameters, race, gender, sexBehs, age, x + a4 * h, ytemp, ak4, year, month, ns, np);

    for (i = 0; i < n; i++)
        ytemp[i] = y[i] + h * (b51 * dydx[i] + b52 * ak2[i] + b53 * ak3[i] + b54 * ak4[i]);
    (*derivs)(PopulationX, Parameters, psaParameters, race, gender, sexBehs, age, x + a5 * h, ytemp, ak5, year, month, ns, np);

    for (i = 0; i < n; i++)
        ytemp[i] = y[i] + h * (b61 * dydx[i] + b62 * ak2[i] + b63 * ak3[i] + b64 * ak4[i] + b65 * ak5[i]);
    (*derivs)(PopulationX, Parameters, psaParameters, race, gender, sexBehs, age, x + a6 * h, ytemp, ak6, year, month, ns, np);

    for (i = 0; i < n; i++) {
        yout[i] = y[i] + h * (c1 * dydx[i] + c3 * ak3[i] + c4 * ak4[i] + c6 * ak6[i]);
    }
    for (i = 0; i < n; i++)
        yerr[i] = h * (dc1 * dydx[i] + dc3 * ak3[i] + dc4 * ak4[i] + dc5 * ak5[i] + dc6 * ak6[i]);

    free_vector(ytemp, 1, n);
    free_vector(ak6, 1, n);
    free_vector(ak5, 1, n);
    free_vector(ak4, 1, n);
    free_vector(ak3, 1, n);
    free_vector(ak2, 1, n);
}

// [[Rcpp::export]]
void rkqs(double** PopulationX, parameters& Parameters, psa_parameters* psaParameters, int race, int gender, int sexBehs, int age,
          double y[], double dydx[], int n, double* x, double htry, double eps, double yscal[], double* hdid, double* hnext, double year, double month,
          void (*derivs)(double**, parameters&, psa_parameters* , int, int, int, int, double, double[], double[], double, double, int, int), int ns, int np)
{
    void rkck(double** PopulationX, parameters & Parameters, psa_parameters * psaParameters, int race, int gender, int sexBehs, int age,
              double y[], double dydx[], int n, double x, double h, double yout[], double yerr[], double year, double month,
              void (*derivs)(double**, parameters&, psa_parameters * , int, int, int, int, double, double[], double[], double, double, int, int), int ns, int np);
    int i;
    double errmax, h, htemp, xnew, * yerr, * ytemp;
    yerr = vector1(1, n);
    ytemp = vector1(1, n);
    h = htry;

    // Set stepsize to the initial trial value.
    for (;;) {
        rkck(PopulationX, Parameters, psaParameters, race, gender, sexBehs, age, y, dydx, n, *x, h, ytemp, yerr, year, month, derivs,ns, np);
        errmax = 0.0;
        for (i = 0; i < n; i++) errmax = fmax(errmax, fabs(yerr[i] / yscal[i]));
        errmax /= eps;
        if (errmax <= 1.0) break;
        htemp = SAFETY * h * pow(errmax, PSHRNK);
        h = (h >= 0.0 ? fmax(htemp, 0.1 * h) : fmin(htemp, 0.1 * h));
        xnew = (*x) + h;
        if (xnew == *x) nrerror("stepsize underflow in rkqs");
    }

    if (errmax > ERRCON) *hnext = SAFETY * h * pow(errmax, PGROW);
    else *hnext = 5.0 * h;
    //No more than a factor of 5 increase.
    *x += (*hdid = h);
    for (i = 0; i < n; i++) {
        y[i] = ytemp[i];
    }
    free_vector(ytemp, 1, n);
    free_vector(yerr, 1, n);
}

//TODO: why are these here? duplication?
#include <math.h>
#define MAXSTP 10000
#define TINY 1.0e-30

// 0.5 for 2 steps
#define tint 1

#define EPS 1.0e-5 
#define H1 1.0e-1
#define HMIN 0

#define KMAXX 200
#define NMAX = 50

// [[Rcpp::export]]
double sumVector(double* Population1, parameters& Parameters, int i_start, int i_end, int j_start, int j_end, int k_start, int k_end, int l_start, int l_end, int m_start, int m_end, int d_start, int d_end) {
    double sum = 0;
    for (int i = i_start; i <= i_end; i++)
        for (int j = j_start; j <= j_end; j++)
            for (int k = k_start; k <= k_end; k++)
                for (int l = l_start; l <= l_end; l++)
                    for (int m = m_start; m <= m_end; m++)
                        for (int d = d_start; d <= d_end; d++)
                            sum = sum + Population1[d + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))))];
    return sum;
}

// [[Rcpp::export]]
double newInfections(double* y, parameters& Parameters, int race, int gender, int sexBeh, int sexAct, int age, int np) {

    //HARDCODED 
    //double beta = 0.5;

    int i = 0.0;
    int j = 0.0;
    int k = 0.0;
    int l = 0.0;
    int m = 0.0;

    // limit only to the population of interest
    if ((age < Parameters.Age_min) || (age > Parameters.Age_max) || (sexAct == 0) || ((gender == 1) & (sexBeh > 0))) return 0.0;
    else
    {
        double susceptible = y[0 + Parameters.nDiseaseStates * (age + Parameters.nAges * (sexAct + Parameters.nSexActs * (sexBeh + Parameters.nSexualBehs * (gender + Parameters.nGenders * race))))];
        double infectedRatio = 0.0;

        switch (gender * Parameters.nSexualBehs + sexBeh) {
        case 0:  //MSW -> WSM
            j = 1; k = 0;
            for (int i = 0; i < Parameters.nRaces; i++)
                for (int l = 1; l < Parameters.nSexActs; l++)
                    for (int m = Parameters.Age_min; m <= Parameters.Age_max; m++) {
                        double totalSubPopulation = 0.0;
                        for (int d = 0; d < Parameters.nDiseaseStates - 2; d++) //minus 2 to remove incidence
                            totalSubPopulation = totalSubPopulation + y[d + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))))];
                        double contact_ratio = 1;
                        infectedRatio = infectedRatio + Parameters.nmbPartners[sexAct] * Parameters.sexFrequency[race][gender][sexBeh][age] * contact_ratio * Parameters.pAges[age][m] * Parameters.pRaces[race][i] * Parameters.pSexActs[sexAct][l] * (y[2 + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))))] + y[3 + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))))]) / totalSubPopulation;
                    }
            break;
        case 1:  //MSM -> MSM, MSMW
            // with MSM
            //HARDCODED 2 - shoul be proportional to homo vs bi
            j = 0; k = 1;
            for (int i = 0; i < Parameters.nRaces; i++)
                for (int l = 1; l < Parameters.nSexActs; l++)
                    for (int m = Parameters.Age_min; m <= Parameters.Age_max; m++) {
                        double totalSubPopulation = 0.0;
                        for (int d = 0; d < Parameters.nDiseaseStates - 2; d++) //minus 2 to remove incidence
                            totalSubPopulation = totalSubPopulation + y[d + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))))];
                        double contact_ratio = (Parameters.popDistribution[i][0][1][1] + Parameters.popDistribution[i][0][1][2]) / (Parameters.popDistribution[i][0][1][1] + Parameters.popDistribution[i][0][2][1] + Parameters.popDistribution[i][0][1][2] + Parameters.popDistribution[i][0][2][2]); //HARDCODED
                        infectedRatio = infectedRatio + Parameters.nmbPartners[sexAct] * Parameters.sexFrequency[race][gender][sexBeh][age] * contact_ratio * Parameters.pAges[age][m] * Parameters.pRaces[race][i] * Parameters.pSexActs[sexAct][l] * (y[2 + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))))] + y[3 + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))))]) / totalSubPopulation;

                    }
            // with MSMW
            j = 0; k = 2;
            for (int i = 0; i < Parameters.nRaces; i++)
                for (int l = 1; l < Parameters.nSexActs; l++)
                    for (int m = Parameters.Age_min; m <= Parameters.Age_max; m++) {
                        double totalSubPopulation = 0.0;
                        for (int d = 0; d < Parameters.nDiseaseStates - 2; d++) //minus 2 to remove incidence
                            totalSubPopulation = totalSubPopulation + y[d + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))))];
                        double contact_ratio = (Parameters.popDistribution[i][0][2][1] + Parameters.popDistribution[i][0][2][2]) / (Parameters.popDistribution[i][0][1][1] + Parameters.popDistribution[i][0][2][1] + Parameters.popDistribution[i][0][1][2] + Parameters.popDistribution[i][0][2][2]); //HARDCODED
                        infectedRatio = infectedRatio + Parameters.nmbPartners[sexAct] * Parameters.sexFrequency[race][gender][sexBeh][age] * contact_ratio * Parameters.pAges[age][m] * Parameters.pRaces[race][i] * Parameters.pSexActs[sexAct][l] * (y[2 + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))))] + y[3 + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))))]) / totalSubPopulation;
                    }
            break;
        case 2:  //MSMW -> MSM, MSMW, WSM
            //HARDCODED 3 - should be proportional to hetero
            // with MSMW
            j = 0; k = 2;
            for (int i = 0; i < Parameters.nRaces; i++)
                for (int l = 1; l < Parameters.nSexActs; l++)
                    for (int m = Parameters.Age_min; m <= Parameters.Age_max; m++) {
                        double totalSubPopulation = 0.0;
                        for (int d = 0; d < Parameters.nDiseaseStates - 2; d++) //minus 2 to remove incidence
                            totalSubPopulation = totalSubPopulation + y[d + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))))];
                        double contact_ratio = (Parameters.popDistribution[i][0][2][1] + Parameters.popDistribution[i][0][2][2]) / (Parameters.popDistribution[i][0][2][1] + Parameters.popDistribution[i][1][0][1] + Parameters.popDistribution[i][0][1][1] + Parameters.popDistribution[i][0][2][2] + Parameters.popDistribution[i][1][0][2] + Parameters.popDistribution[i][0][1][2]); //HARDCODED
                        infectedRatio = infectedRatio + Parameters.nmbPartners[sexAct] * Parameters.sexFrequency[race][gender][sexBeh][age] * contact_ratio * Parameters.pAges[age][m] * Parameters.pRaces[race][i] * Parameters.pSexActs[sexAct][l] * (y[2 + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))))] + y[3 + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))))]) / totalSubPopulation;
                    }
            // with WSM
            j = 1; k = 0;
            for (int i = 0; i < Parameters.nRaces; i++)
                for (int l = 1; l < Parameters.nSexActs; l++)
                    for (int m = Parameters.Age_min; m <= Parameters.Age_max; m++) {
                        double totalSubPopulation = 0.0;
                        for (int d = 0; d < Parameters.nDiseaseStates - 2; d++) //minus 2 to remove incidence
                            totalSubPopulation = totalSubPopulation + y[d + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))))];
                        double contact_ratio = (Parameters.popDistribution[i][1][0][1] + Parameters.popDistribution[i][1][0][2]) / (Parameters.popDistribution[i][0][2][1] + Parameters.popDistribution[i][1][0][1] + Parameters.popDistribution[i][0][1][1] + Parameters.popDistribution[i][0][2][2] + Parameters.popDistribution[i][1][0][2] + Parameters.popDistribution[i][0][1][2]); //HARDCODED
                        infectedRatio = infectedRatio + Parameters.nmbPartners[sexAct] * Parameters.sexFrequency[race][gender][sexBeh][age] * contact_ratio * Parameters.pAges[age][m] * Parameters.pRaces[race][i] * Parameters.pSexActs[sexAct][l] * (y[2 + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))))] + y[3 + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))))]) / totalSubPopulation;
                    }
            // with MSM
            j = 0; k = 1;
            for (int i = 0; i < Parameters.nRaces; i++)
                for (int l = 1; l < Parameters.nSexActs; l++)
                    for (int m = Parameters.Age_min; m <= Parameters.Age_max; m++) {
                        double totalSubPopulation = 0.0;
                        for (int d = 0; d < Parameters.nDiseaseStates - 2; d++) //minus 2 to remove incidence
                            totalSubPopulation = totalSubPopulation + y[d + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))))];
                        double contact_ratio = (Parameters.popDistribution[i][0][1][1] + Parameters.popDistribution[i][0][1][2]) / (Parameters.popDistribution[i][0][2][1] + Parameters.popDistribution[i][1][0][1] + Parameters.popDistribution[i][0][1][1] + Parameters.popDistribution[i][0][2][2] + Parameters.popDistribution[i][1][0][2] + Parameters.popDistribution[i][0][1][2]); //HARDCODED
                        infectedRatio = infectedRatio + Parameters.nmbPartners[sexAct] * Parameters.sexFrequency[race][gender][sexBeh][age] * contact_ratio * Parameters.pAges[age][m] * Parameters.pRaces[race][i] * Parameters.pSexActs[sexAct][l] * (y[2 + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))))] + y[3 + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))))]) / totalSubPopulation;
                    }
            break;
        case 3:  //WSM -> MSW, MSMW
            // with MSW
            j = 0; k = 1;
            for (int i = 0; i < Parameters.nRaces; i++)
                for (int l = 1; l < Parameters.nSexActs; l++)
                    for (int m = Parameters.Age_min; m <= Parameters.Age_max; m++) {
                        double totalSubPopulation = 0.0;
                        for (int d = 0; d < Parameters.nDiseaseStates - 2; d++) //minus 2 to remove incidence
                            totalSubPopulation = totalSubPopulation + y[d + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))))];
                        double contact_ratio = (Parameters.popDistribution[i][0][1][1] + Parameters.popDistribution[i][0][1][2]) / (Parameters.popDistribution[i][0][1][1] + Parameters.popDistribution[i][0][2][1] + Parameters.popDistribution[i][0][1][2] + Parameters.popDistribution[i][0][2][2]); //HARDCODED
                        infectedRatio = infectedRatio + Parameters.nmbPartners[sexAct] * Parameters.sexFrequency[race][gender][sexBeh][age] * contact_ratio * Parameters.pAges[age][m] * Parameters.pRaces[race][i] * Parameters.pSexActs[sexAct][l] * (y[2 + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))))] + y[3 + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))))]) / totalSubPopulation;
                    }
            // with MSMW
            j = 0; k = 2;
            for (int i = 0; i < Parameters.nRaces; i++)
                for (int l = 1; l < Parameters.nSexActs; l++)
                    for (int m = Parameters.Age_min; m <= Parameters.Age_max; m++) {
                        double totalSubPopulation = 0.0;
                        for (int d = 0; d < Parameters.nDiseaseStates - 2; d++) //minus 2 to remove incidence
                            totalSubPopulation = totalSubPopulation + y[d + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))))];
                        double contact_ratio = (Parameters.popDistribution[i][0][2][1] + Parameters.popDistribution[i][0][2][2]) / (Parameters.popDistribution[i][0][1][1] + Parameters.popDistribution[i][0][2][1] + Parameters.popDistribution[i][0][1][2] + Parameters.popDistribution[i][0][2][2]); //HARDCODED
                        infectedRatio = infectedRatio + Parameters.nmbPartners[sexAct] * Parameters.sexFrequency[race][gender][sexBeh][age] * contact_ratio * Parameters.pAges[age][m] * Parameters.pRaces[race][i] * Parameters.pSexActs[sexAct][l] * (y[2 + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))))] + y[3 + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))))]) / totalSubPopulation;
                    }
            break;
        }
        //HARDCODED

       //determine the calibration age group
       int age_group = 0;
       for (int m1 = 0; m1 < Parameters.ageGroup_size; m1++)
           if (age >= Parameters.ageGroup_bound[m1]) age_group = m1;
       
        return Parameters.transmissionRate[np][race][gender][sexBeh][age_group] * susceptible * (infectedRatio);
       //return 0.022 * susceptible * (infectedRatio);
    }
}

// [[Rcpp::export]]
void derivs(double** PopulationX, parameters& Parameters, psa_parameters* psaParameters, int race, int gender, int sexBehs, int age, double x, double* y, double* dydx, double year, double month, int ns, int np)
{
    //HARDCODED
    double aging_rate = 1.0 / 12.0;
    int time = year * 12 + month;
    int max_tau = (time <= 20) ? time : 20;
    
    for (int i = 0; i < Parameters.nRaces; i++) {
        double newbornePopulation = Parameters.birthRate[i][year] * sumVector(y, Parameters, i, i, 0, Parameters.nGenders - 1, 0, Parameters.nSexualBehs - 1, 0, Parameters.nSexActs - 1, 0, Parameters.nAges - 1, 0, Parameters.nDiseaseStates - 3);//to remove incidence

        for (int j = 0; j < Parameters.nGenders; j++)
            for (int k = 0; k < Parameters.nSexualBehs; k++)
                for (int l = 0; l < Parameters.nSexActs; l++)
                    for (int m = 0; m < Parameters.nAges; m++) {

                        double Incidence = newInfections(y, Parameters, i, j, k, l, m, np);

                        for (int d = 0; d < Parameters.nDiseaseStates; d++) {
                            int index_a1 = d + Parameters.nDiseaseStates * (m - 1 + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));
                            int index = d + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));
                            
                            double totalWaning = 0.0;
     
                            switch (d) {
                            case 0:
                                //susceptible
                                totalWaning = 0.0;
                                //for (int tau = 1; tau <= max_tau; tau++)
                                //    totalWaning = totalWaning + wanA1[tau] * PopulationX[time - tau][index + 5] + wanA2[tau] * PopulationX[time - tau][index + 6] + wanA3[tau] * PopulationX[time - tau][index + 7] + wanA4[tau] * PopulationX[time - tau][index + 8]
                                //                              + wanB2[tau] * PopulationX[time - tau][index + 6] + wanB3[tau] * PopulationX[time - tau][index + 6] + wanB4[tau] * PopulationX[time - tau][index + 6];
                                
                                dydx[index] = ((m == 0) ? newbornePopulation * Parameters.popDistribution[i][j][k][l] * (1 - Parameters.deathRate[i][j][m]) : aging_rate * y[index_a1]) 
                                             + psaParameters[ns].recoveryRate[j][k] * y[index + 4] + psaParameters[ns].recoveryAsympRate[j][k] * y[index + 2] 
                                             //+ totalWaning
                                             - Incidence 
                                             //- (vaccA1[i][j][k][l][m] + vaccA2[i][j][k][l][m] + vaccA3[i][j][k][l][m]+ vaccA4[i][j][k][l][m]+ vaccB2[i][j][k][l][m]+ vaccB3[i][j][k][l][m]+ vaccB4[m]) * y[index]
                                             - (aging_rate + Parameters.deathRate[i][j][m]) * y[index];
                                break;
                            case 1:
                                //exposed
                                dydx[index] = ((m == 0) ? 0 : aging_rate * y[index_a1]) 
                                            + Incidence 
                                            //- (vaccA1[i][j][k][l][m] + vaccA2[i][j][k][l][m] + vaccA3[i][j][k][l][m]+ vaccA4[i][j][k][l][m]+ vaccB2[i][j][k][l][m]+ vaccB3[i][j][k][l][m]+ vaccB4[m]) * y[index]
                                            - (aging_rate + psaParameters[ns].latencyRate[j][k] + Parameters.deathRate[i][j][m]) * y[index];
                                break;
                            case 2:
                                //asymptomatic
                                dydx[index] = ((m == 0) ? 0 : aging_rate * y[index_a1]) 
                                            + psaParameters[ns].propAsymp[j][k] * psaParameters[ns].latencyRate[j][k] * y[index - 1] 
                                            //- (vaccA1[i][j][k][l][m] + vaccA2[i][j][k][l][m] + vaccA3[i][j][k][l][m]+ vaccA4[i][j][k][l][m]+ vaccB2[i][j][k][l][m]+ vaccB3[i][j][k][l][m]+ vaccB4[m]) * y[index]
                                            - (aging_rate + psaParameters[ns].screeningRate[i][j][k][m] + psaParameters[ns].recoveryAsympRate[j][k] + Parameters.deathRate[i][j][m]) * y[index];
                                break;
                            case 3:
                                //symptomatic
                                dydx[index] = ((m == 0) ? 0 : aging_rate * y[index_a1]) + (1 - psaParameters[ns].propAsymp[j][k]) * psaParameters[ns].latencyRate[j][k] * y[index - 2] - (aging_rate + psaParameters[ns].tretmentRate[j][k] + Parameters.deathRate[i][j][m]) * y[index];
                                break;
                            case 4:
                                //treated
                                dydx[index] = ((m == 0) ? 0 : aging_rate * y[index_a1]) + psaParameters[ns].screeningRate[i][j][k][m] * y[index - 2] + psaParameters[ns].tretmentRate[j][k] * y[index - 1] - (aging_rate + psaParameters[ns].recoveryRate[j][k] + Parameters.deathRate[i][j][m]) * y[index];
                                break;
                            case 5:
                                //vacc schedule A dose 1
                                totalWaning = 0.0;
                                //for (int tau = 1; tau <= max_tau; tau++) totalWaning = totalWaning + wanA1[tau] * PopulationX[time - tau][index];
                                dydx[index] = ((m == 0) ? 0 : aging_rate * y[index_a1])
                                    //+ vaccA1[i][j][k][l][m] * y[index-5] + vaccA1[i][j][k][l][m] * y[index - 4] + vaccA1[i][j][k][l][m] * y[index - 3]
                                    //- vaccA2[i][j][k][l][m]*y[index]
                                    //- vaccB2[i][j][k][l][m]*y[index]
                                    //- totalWaning;
                                    - (aging_rate + Parameters.deathRate[i][j][m]) * y[index];
                                break;
                            case 6:
                                //vacc schedule A dose 2
                                totalWaning = 0.0;
                                //for (int tau = 1; tau <= max_tau; tau++) totalWaning = totalWaning + wanA2[tau] * PopulationX[time - tau][index];
                                dydx[index] = ((m == 0) ? 0 : aging_rate * y[index_a1])
                                    //+ vaccA2[i][j][k][l][m] * y[index-5] + vaccA2[i][j][k][l][m] * y[index - 4] + vaccA2[i][j][k][l][m] * y[index - 3]
                                    //+ vaccA2[i][j][k][l][m] * y[index-1]
                                    //- vaccA3[i][j][k][l][m] * y[index]
                                    //- totalWaning;
                                    - (aging_rate + Parameters.deathRate[i][j][m]) * y[index];
                                break;
                            case 7:
                                //vacc schedule A dose 3
                                totalWaning = 0.0;
                                //for (int tau = 1; tau <= max_tau; tau++) totalWaning = totalWaning + wanA3[tau] * PopulationX[time - tau][index];
                                dydx[index] = ((m == 0) ? 0 : aging_rate * y[index_a1])
                                    //+ vaccA3[i][j][k][l][m] * y[index-5] + vaccA3[i][j][k][l][m] * y[index - 4] + vaccA3[i][j][k][l][m] * y[index - 3]
                                    //+ vaccA3[i][j][k][l][m] * y[index-1]
                                    //- vaccA4[i][j][k][l][m] * y[index]
                                    //- totalWaning;
                                    - (aging_rate + Parameters.deathRate[i][j][m]) * y[index];
                                break;
                            case 8://vacc schedule A dose 4
                                totalWaning = 0.0;
                                //for (int tau = 1; tau <= max_tau; tau++) totalWaning = totalWaning + wanA4[tau] * PopulationX[time - tau][index];
                                dydx[index] = ((m == 0) ? 0 : aging_rate * y[index_a1])
                                    //+ vaccA4[i][j][k][l][m] * y[index-5] + vaccA4[i][j][k][l][m] * y[index - 4] + vaccA4[i][j][k][l][m] * y[index - 3]
                                    //+ vaccA4[i][j][k][l][m] * y[index-1]
                                    //- totalWaning;
                                    - (aging_rate + Parameters.deathRate[i][j][m]) * y[index];
                                break;
                            case 9:
                                //vacc schedule B dose 2
                                totalWaning = 0.0;
                                //for (int tau = 1; tau <= max_tau; tau++) totalWaning = totalWaning + wanB2[tau] * PopulationX[time - tau][index];
                                dydx[index] = ((m == 0) ? 0 : aging_rate * y[index_a1])
                                    //+ vaccB2[i][j][k][l][m] * y[index-5] + vaccB2[i][j][k][l][m] * y[index - 4] + vaccB2[i][j][k][l][m] * y[index - 3]
                                    //+ vaccB2[i][j][k][l][m]*y[index-4]
                                    //- vaccB3[i][j][k][l][m] * y[index]
                                    //- totalWaning;
                                    - (aging_rate + Parameters.deathRate[i][j][m]) * y[index];
                                break;
                            case 10:
                                //vacc schedule b dose 3
                                totalWaning = 0.0;
                                //for (int tau = 1; tau <= max_tau; tau++) totalWaning = totalWaning + wanB3[tau] * PopulationX[time - tau][index];
                                dydx[index] = ((m == 0) ? 0 : aging_rate * y[index_a1])
                                    //+ vaccB3[i][j][k][l][m] * y[index-5] + vaccB3[i][j][k][l][m] * y[index - 4] + vaccB3[i][j][k][l][m] * y[index - 3]
                                    //+ vaccB3[i][j][k][l][m] * y[index-1]
                                    //- vaccB4[i][j][k][l][m] * y[index]
                                    - (aging_rate + Parameters.deathRate[i][j][m]) * y[index];
                                break;
                            case 11:
                                //vacc schedule B dose 4
                                totalWaning = 0.0;
                                //for (int tau = 1; tau <= max_tau; tau++) totalWaning = totalWaning + wanB4[tau] * PopulationX[time - tau][index];
                                dydx[index] = ((m == 0) ? 0 : aging_rate * y[index_a1])
                                    //+ vaccB4[i][j][k][l][m] * y[index-5] + vaccB4[i][j][k][l][m] * y[index - 4] + vaccB4[i][j][k][l][m] * y[index - 3]
                                    //+ vaccB4[i][j][k][l][m] * y[index-1]
                                    - (aging_rate + Parameters.deathRate[i][j][m]) * y[index];
                                break;
                            case 12: // placeholder
                                dydx[index] = 0;
                                break;
                            case 13: // Reported incidence
                                dydx[index] = psaParameters[ns].tretmentRate[j][k] * y[index - 10] + psaParameters[ns].screeningRate[i][j][k][m] * y[index - 11];;
                                break;
                            }
                        }
                    }
    }
}

// [[Rcpp::export]]
void odeint(double** PopulationX, parameters& Parameters, psa_parameters* psaParameters, int race, int gender, int sexBehs, int age, double* ystart, double year, double month, int size, int ns, int np) {

    /* Runge - Kutta driver with adaptive stepsize control.Integrates the starting values ystart(1:NVAR) from x1 to x2 _
       with accuracy eps, storing intermediate results in the global variables xp(1:KMAXX), yp(1:NMAX, 1 : KMAXX).h1 should _
       be set as a guessed first step, hmin as the minimum allowed stepsize(can be zero).On output nokand nbad are the _
       number of goodand bad(but retriedand fixed) steps taken, and ystart(1:NVAR) is replaced by values at the end of _
       the integration interval.derivs is the user supplied subroutine for calculating the right - hand side derivative, _
       while rkqs is the name of the stepper routine to be used. */

    double eps = EPS;
    double h1 = H1;
    double hmin = HMIN;
    double kmax = KMAXX;
    int nvar = size;

    double x1 = 0;
    double x2 = x1 + tint;
    double dxsav = (x2 - x1) / 100.0;

    vector<double> xp(KMAXX, 0.0);
    vector<vector<double>> yp(nvar, vector<double>(KMAXX, 0.0));  //???

    int nstp, i;
    int kount = 0;
    double xsav, x, hnext, hdid, h;
    double* yscal, * y, * dydx;
    yscal = vector1(1, nvar);
    y = vector1(1, nvar);
    dydx = vector1(1, nvar);
    x = x1;
    h = SIGN(h1, x2 - x1);
    for (i = 0; i < nvar; i++)
        y[i] = ystart[i];
    if (kmax > 0) xsav = x - dxsav * 2.0;

    for (nstp = 1; nstp <= MAXSTP; nstp++) {
        (*derivs)(PopulationX, Parameters, psaParameters, race, gender, sexBehs, age, x, y, dydx, year, month, ns, np);
        for (i = 0; i < nvar; i++) {
            yscal[i] = fabs(y[i]) + fabs(dydx[i] * h) + TINY;
        }

        if (kmax > 0 && kount < kmax - 1 && fabs(x - xsav) > fabs(dxsav)) {
            xp[++kount] = x;
            for (i = 0; i < nvar; i++) yp[i][kount] = y[i];
            xsav = x;
        }
        if ((x + h - x2) * (x + h - x1) > 0.0) h = x2 - x;

        (*rkqs)(PopulationX, Parameters, psaParameters, race, gender, sexBehs, age, y, dydx, nvar, &x, h, eps, yscal, &hdid, &hnext, year, month, derivs, ns, np);
        if ((x - x2) * (x2 - x1) >= 0.0) {

            for (i = 0; i < nvar; i++) ystart[i] = y[i];
            if (kmax) {
                xp[++kount] = x;
                for (i = 0; i < nvar; i++) yp[i][kount] = y[i];
            }
            free_vector(dydx, 1, nvar);
            free_vector(y, 1, nvar);
            free_vector(yscal, 1, nvar);
            return;

        }
        if (fabs(hnext) <= hmin) nrerror("Step size too small in odeint");
        h = hnext;
    }
    nrerror("Too many steps in routine odeint");
}

// [[Rcpp::export]]
void loadInitialPopulation(std::string inputPath, double*** PopulationX, parameters& Parameters, int NumbParallel, double initialInfection) {
    ifstream ifile(inputPath, ios::in);

    //check to see that the file was opened correctly
    if (!ifile.is_open()) {
        cerr << "There was a problem opening the input file!-dem\n";
        exit(1);
    }

    std::vector<std::vector<double>> vec;
    std::string lineData;

    //load initial population distribution
    double totalBlack = 0.0;
    double totalHispanic = 0.0;
    double totalWhite = 0.0;

    double pop;
    int idx = 0;
    for (int d = 0; d < Parameters.nDiseaseStates; d++)
        for (int i = 0; i < Parameters.nRaces; i++)
            for (int j = 0; j < Parameters.nGenders; j++)
                for (int k = 0; k < Parameters.nSexualBehs; k++)
                    for (int l = 0; l < Parameters.nSexActs; l++)
                        for (int m = 0; m < Parameters.nAges; m++) {
                            if (d >= Parameters.nDiseaseStates - 2) {
                                int index = d + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));
                                PopulationX[0][0][d + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))))] = 0.0;
                                for (int np = 1; np < NumbParallel; np++) PopulationX[np][0][index] = PopulationX[0][0][index];
                            }
                            else {
                                idx++;
                                std::getline(ifile, lineData);
                                pop = stod(lineData);
                                int index = d + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));

                                PopulationX[0][0][d + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))))] = pop;
                                for (int np = 1; np < NumbParallel; np++) PopulationX[np][0][index] = PopulationX[0][0][index];

                                if (i == 0) totalBlack = totalBlack + pop;
                                else
                                    if (i == 1) totalHispanic = totalHispanic + pop;
                                    else totalWhite = totalWhite + pop;
                            }
                        }
    
    //calculate newborn and sexual contact distribution  - static    
    double sumB = 0.0;
    double sumH = 0.0;
    double sumW = 0.0;
    for (int i = 0; i < Parameters.nRaces; i++) {
        Parameters.popDistribution.push_back(vector<vector<vector<double>>>());
        for (int j = 0; j < Parameters.nGenders; j++) {
            Parameters.popDistribution[i].push_back(vector<vector<double>>());
            for (int k = 0; k < Parameters.nSexualBehs; k++) {
                Parameters.popDistribution[i][j].push_back(vector<double>());
                for (int l = 0; l < Parameters.nSexActs; l++) {
                    double pom = sumVector(PopulationX[0][0], Parameters, i, i, j, j, k, k, l, l, 0, Parameters.nAges - 1, 0, Parameters.nDiseaseStates - 3);  // - 3 to remove incidence
                    if (i == 0) { Parameters.popDistribution[i][j][k].push_back(pom / totalBlack); sumB = Parameters.popDistribution[i][j][k][l] + sumB; }
                    else
                        if (i == 1) {
                            Parameters.popDistribution[i][j][k].push_back(pom / totalHispanic); sumH = Parameters.popDistribution[i][j][k][l] + sumH;
                        }
                        else {
                            Parameters.popDistribution[i][j][k].push_back(pom / totalWhite); sumW = Parameters.popDistribution[i][j][k][l] + sumW;
                        }
                }
            }
        }
    }//rounding error correction
    Parameters.popDistribution[0][Parameters.nGenders - 1][Parameters.nSexualBehs - 1][Parameters.nSexActs - 1] = Parameters.popDistribution[0][Parameters.nGenders - 1][Parameters.nSexualBehs - 1][Parameters.nSexActs - 1] + (1 - sumB);
    Parameters.popDistribution[1][Parameters.nGenders - 1][Parameters.nSexualBehs - 1][Parameters.nSexActs - 1] = Parameters.popDistribution[1][Parameters.nGenders - 1][Parameters.nSexualBehs - 1][Parameters.nSexActs - 1] + (1 - sumH);
    Parameters.popDistribution[2][Parameters.nGenders - 1][Parameters.nSexualBehs - 1][Parameters.nSexActs - 1] = Parameters.popDistribution[2][Parameters.nGenders - 1][Parameters.nSexualBehs - 1][Parameters.nSexActs - 1] + (1 - sumW);
  ifile.close();
}

// [[Rcpp::export]]
void loadDemographics(std::string inputPath, parameters& Parameters) {
    ifstream ifile(inputPath, ios::in);
    
    //check to see that the file was opened correctly:
    if (!ifile.is_open()) {
        cerr << "There was a problem opening the input file!-dem\n";
        exit(1);
    }
    
    std::vector<std::vector<double>> vec;
    std::string lineData;

    //load age group structure for calibration
    getline(ifile, lineData);  Parameters.Age_min = stoi(lineData);
    getline(ifile, lineData);  Parameters.Age_max = stoi(lineData);

    getline(ifile, lineData);  Parameters.ageGroup_bound.push_back(stoi(lineData));
    getline(ifile, lineData);  Parameters.ageGroup_bound.push_back(stoi(lineData));
    getline(ifile, lineData);  Parameters.ageGroup_bound.push_back(stoi(lineData));
    getline(ifile, lineData);  Parameters.ageGroup_bound.push_back(stoi(lineData));

    //load yearly birth rate - per capita
    for (int i = 0; i < Parameters.nRaces; i++) {
        double d;
        vector<double> row;
        std::getline(ifile, lineData);
        stringstream lineStream(lineData);
        while (lineStream >> d) 
            row.push_back(d);
        Parameters.birthRate.push_back(row);
    }
    
    //load monthly death rate 
    for (int i = 0; i < Parameters.nRaces; i++) {
        Parameters.deathRate.push_back(vector<vector<double>>());
        for (int j = 0; j < Parameters.nGenders; j++) {
            double d;
            vector<double> row;
            std::getline(ifile, lineData);
            stringstream lineStream(lineData);
            while (lineStream >> d) row.push_back(d);
            Parameters.deathRate[i].push_back(row);
        }
    }
   
    //load assortative mixing for race/ethnicity
    double d;
    vector<double> row;
    std::getline(ifile, lineData);
    stringstream lineStream(lineData);
    
    //load assortative mixing for race/ethnicity
    lineStream >> d;
    for (int i = 0; i < Parameters.nRaces; i++) {
        Parameters.pRaces.push_back(vector<double>());
        for (int i1 = 0; i1 < Parameters.nRaces; i1++)
            if (i == i1) Parameters.pRaces[i].push_back(d);
            else Parameters.pRaces[i].push_back((1 - d) / 2);  //HARDCODED 2 to equally split proportion of contact among other race/ethnic groups
    }
    
    //load assortative mixing for sexual activity
    lineStream >> d;
    for (int l = 0; l < Parameters.nSexActs; l++) {
        Parameters.pSexActs.push_back(vector<double>());
        for (int l1 = 0; l1 < Parameters.nSexActs; l1++)
            if (l == 0 || l1 == 0)Parameters.pSexActs[l].push_back(0);
            else
                if (l == l1) Parameters.pSexActs[l].push_back(d);
                else Parameters.pSexActs[l].push_back((1-d));
    }

    //load assortative mixing for age group between Amin and Amax years of age with defined ageBand
    lineStream >> d;

    for (int m = 0; m < Parameters.nAges; m++) {
        Parameters.pAges.push_back(vector<double>());
        for (int m1 = 0; m1 < Parameters.nAges; m1++) {
            if ((m >= Parameters.Age_min) & (m <= Parameters.Age_max) & (m1 >= Parameters.Age_min) & (m1 <= Parameters.Age_max))
            {
                //determine the calibration age group
                int age_group = 0;
                int age_group1 = 0;

                for (int m2 = 0; m2 < Parameters.ageGroup_size; m2++) {
                    if (m >= Parameters.ageGroup_bound[m2]) age_group = m2;
                    if (m1 >= Parameters.ageGroup_bound[m2]) age_group1 = m2;
                }

                if (age_group == age_group1) Parameters.pAges[m].push_back(d);
                else Parameters.pAges[m].push_back(1 - d);
            }
            else Parameters.pAges[m].push_back(0.0);
        }
    }

    // load average number of partners
    // HARDCODED for testing
    Parameters.nmbPartners.push_back(0.0);
    Parameters.nmbPartners.push_back(1.0);
    Parameters.nmbPartners.push_back(1.0);

    //load mean number of sexual partners per month by age group
    for (int i = 0; i < Parameters.nRaces; i++) {
        Parameters.sexFrequency.push_back(vector<vector<vector<double>>>());
        for (int j = 0; j < Parameters.nGenders; j++) {
            Parameters.sexFrequency[i].push_back(vector<vector<double>>());
            for (int k = 0; k < Parameters.nSexualBehs; k++) {
                double d;
                vector<double> row;
                std::getline(ifile, lineData);
                stringstream lineStream(lineData);
                while (lineStream >> d) row.push_back(d/12);
                Parameters.sexFrequency[i][j].push_back(row);
            }
        }
    }
    ifile.close();
}

// [[Rcpp::export]]
void loadParameters(std::string inputPath, parameters& Parameters, psa_parameters* psaParameters, int NumbSim) {
    ifstream ifile(inputPath, ios::in);

    //check to see that the file was opened correctly:
    if (!ifile.is_open()) {
        cerr << "There was a problem opening the input file - Parameters!\n";
        exit(1);
    }
        
    std::string lineData;
    
    //load vaccination.rate.A1 
    for (int i = 0; i < Parameters.nRaces; i++) {
        Parameters.vaccA1.push_back(vector<vector<vector<vector<double>>>>());
        for (int j = 0; j < Parameters.nGenders; j++) {
            Parameters.vaccA1[i].push_back(vector<vector<vector<double>>>());
            for (int k = 0; k < Parameters.nSexualBehs; k++) {
                Parameters.vaccA1[i][j].push_back(vector<vector<double>>());
                for (int l = 0; l < Parameters.nSexActs; l++) {
                    double d;
                    vector<double> row;
                    getline(ifile, lineData);
                    stringstream lineStream(lineData);
                    while (lineStream >> d) row.push_back(d);
                    Parameters.vaccA1[i][j][k].push_back(row);
                }
            }
        }
    }

    //load vaccination.rate.A2
    for (int i = 0; i < Parameters.nRaces; i++) {
        Parameters.vaccA2.push_back(vector<vector<vector<vector<double>>>>());
        for (int j = 0; j < Parameters.nGenders; j++) {
            Parameters.vaccA2[i].push_back(vector<vector<vector<double>>>());
            for (int k = 0; k < Parameters.nSexualBehs; k++) {
                Parameters.vaccA2[i][j].push_back(vector<vector<double>>());
                for (int l = 0; l < Parameters.nSexActs; l++) {
                    double d;
                    vector<double> row;
                    getline(ifile, lineData);
                    stringstream lineStream(lineData);
                    while (lineStream >> d) row.push_back(d);
                    Parameters.vaccA2[i][j][k].push_back(row);
                }
            }
        }
    }

    //load vaccination.rate.A3 
    for (int i = 0; i < Parameters.nRaces; i++) {
        Parameters.vaccA3.push_back(vector<vector<vector<vector<double>>>>());
        for (int j = 0; j < Parameters.nGenders; j++) {
            Parameters.vaccA3[i].push_back(vector<vector<vector<double>>>());
            for (int k = 0; k < Parameters.nSexualBehs; k++) {
                Parameters.vaccA3[i][j].push_back(vector<vector<double>>());
                for (int l = 0; l < Parameters.nSexActs; l++) {
                    double d;
                    vector<double> row;
                    getline(ifile, lineData);
                    stringstream lineStream(lineData);
                    while (lineStream >> d) row.push_back(d);
                    Parameters.vaccA3[i][j][k].push_back(row);
                }
            }
        }
    }

    //load vaccination.rate.A4 
    for (int i = 0; i < Parameters.nRaces; i++) {
        Parameters.vaccA4.push_back(vector<vector<vector<vector<double>>>>());
        for (int j = 0; j < Parameters.nGenders; j++) {
            Parameters.vaccA4[i].push_back(vector<vector<vector<double>>>());
            for (int k = 0; k < Parameters.nSexualBehs; k++) {
                Parameters.vaccA4[i][j].push_back(vector<vector<double>>());
                for (int l = 0; l < Parameters.nSexActs; l++) {
                    double d;
                    vector<double> row;
                    getline(ifile, lineData);
                    stringstream lineStream(lineData);
                    while (lineStream >> d) row.push_back(d);
                    Parameters.vaccA4[i][j][k].push_back(row);
                }
            }
        }
    }

    //load vaccination.rate.B2 
    for (int i = 0; i < Parameters.nRaces; i++) {
        Parameters.vaccB2.push_back(vector<vector<vector<vector<double>>>>());
        for (int j = 0; j < Parameters.nGenders; j++) {
            Parameters.vaccB2[i].push_back(vector<vector<vector<double>>>());
            for (int k = 0; k < Parameters.nSexualBehs; k++) {
                Parameters.vaccB2[i][j].push_back(vector<vector<double>>());
                for (int l = 0; l < Parameters.nSexActs; l++) {
                    double d;
                    vector<double> row;
                    getline(ifile, lineData);
                    stringstream lineStream(lineData);
                    while (lineStream >> d) row.push_back(d);
                    Parameters.vaccB2[i][j][k].push_back(row);
                }
            }
        }
    }

    //load vaccination.rate.B3 
    for (int i = 0; i < Parameters.nRaces; i++) {
        Parameters.vaccB3.push_back(vector<vector<vector<vector<double>>>>());
        for (int j = 0; j < Parameters.nGenders; j++) {
            Parameters.vaccB3[i].push_back(vector<vector<vector<double>>>());
            for (int k = 0; k < Parameters.nSexualBehs; k++) {
                Parameters.vaccB3[i][j].push_back(vector<vector<double>>());
                for (int l = 0; l < Parameters.nSexActs; l++) {
                    double d;
                    vector<double> row;
                    getline(ifile, lineData);
                    stringstream lineStream(lineData);
                    while (lineStream >> d) row.push_back(d);
                    Parameters.vaccB3[i][j][k].push_back(row);
                }
            }
        }
    }

    //load vaccination.rate.B4
    for (int i = 0; i < Parameters.nRaces; i++) {
        Parameters.vaccB4.push_back(vector<vector<vector<vector<double>>>>());
        for (int j = 0; j < Parameters.nGenders; j++) {
            Parameters.vaccB4[i].push_back(vector<vector<vector<double>>>());
            for (int k = 0; k < Parameters.nSexualBehs; k++) {
                Parameters.vaccB4[i][j].push_back(vector<vector<double>>());
                for (int l = 0; l < Parameters.nSexActs; l++) {
                    double d;
                    vector<double> row;
                    getline(ifile, lineData);
                    stringstream lineStream(lineData);
                    while (lineStream >> d) row.push_back(d);
                    Parameters.vaccB4[i][j][k].push_back(row);
                }
            }
        }
    }

    //load waning.immunity
    double d;
    getline(ifile, lineData); stringstream lineStreamA1(lineData); while (lineStreamA1 >> d) Parameters.waningImmunityA1.push_back(d);
    getline(ifile, lineData); stringstream lineStreamA2(lineData); while (lineStreamA2 >> d) Parameters.waningImmunityA2.push_back(d);
    getline(ifile, lineData); stringstream lineStreamA3(lineData); while (lineStreamA3 >> d) Parameters.waningImmunityA3.push_back(d);
    getline(ifile, lineData); stringstream lineStreamA4(lineData); while (lineStreamA4 >> d) Parameters.waningImmunityA4.push_back(d);
    getline(ifile, lineData); stringstream lineStreamB2(lineData); while (lineStreamB2 >> d) Parameters.waningImmunityB2.push_back(d);
    getline(ifile, lineData); stringstream lineStreamB3(lineData); while (lineStreamB3 >> d) Parameters.waningImmunityB3.push_back(d);
    getline(ifile, lineData); stringstream lineStreamB4(lineData); while (lineStreamB4 >> d) Parameters.waningImmunityB4.push_back(d);
    
    //load other parameters
    for (int ns = 0; ns < NumbSim; ns++)
    {
        //load latency rate
        for (int j = 0; j < Parameters.nGenders; j++) {
            psaParameters[ns].latencyRate.push_back(vector<double>());
            for (int k = 0; k < Parameters.nSexualBehs; k++) {
                getline(ifile, lineData); psaParameters[ns].latencyRate[j].push_back(1 / stod(lineData));
            }
        }

        //load proportion of asymptomatic
        for (int j = 0; j < Parameters.nGenders; j++) {
            psaParameters[ns].propAsymp.push_back(vector<double>());
            for (int k = 0; k < Parameters.nSexualBehs; k++) {
                getline(ifile, lineData); psaParameters[ns].propAsymp[j].push_back(stod(lineData));
            }
        }

        //load treatment rate
        for (int j = 0; j < Parameters.nGenders; j++) {
            psaParameters[ns].tretmentRate.push_back(vector<double>());
            for (int k = 0; k < Parameters.nSexualBehs; k++) {
                getline(ifile, lineData); psaParameters[ns].tretmentRate[j].push_back(1 / stod(lineData));
            }
        }

        //load recovery rate
        for (int j = 0; j < Parameters.nGenders; j++) {
            psaParameters[ns].recoveryRate.push_back(vector<double>());
            for (int k = 0; k < Parameters.nSexualBehs; k++) {
                getline(ifile, lineData); psaParameters[ns].recoveryRate[j].push_back(1 / stod(lineData));
            }
        }

        //load recovery asymptgomatic rate
        for (int j = 0; j < Parameters.nGenders; j++) {
            psaParameters[ns].recoveryAsympRate.push_back(vector<double>());
            for (int k = 0; k < Parameters.nSexualBehs; k++) {
                getline(ifile, lineData); psaParameters[ns].recoveryAsympRate[j].push_back(1 / stod(lineData));
            }
        }
        
        //load screening rate
        for (int i = 0; i < Parameters.nRaces; i++) {
            psaParameters[ns].screeningRate.push_back(vector<vector<vector<double>>>());
            for (int j = 0; j < Parameters.nGenders; j++) {
                psaParameters[ns].screeningRate[i].push_back(vector<vector<double>>());
                for (int k = 0; k < Parameters.nSexualBehs; k++) {
                    double d;
                    vector<double> row;
                    getline(ifile, lineData);
                    stringstream lineStream(lineData);
                    while (lineStream >> d) row.push_back(d);
                    psaParameters[ns].screeningRate[i][j].push_back(row);
                }
            }
        }
    }
    ifile.close();
}

// [[Rcpp::export]]
void loadCalibratioTargets(std::string inputPath, parameters& Parameters) {
    ifstream ifile(inputPath, ios::in);

    //check to see that the file was opened correctly:
    if (!ifile.is_open()) {
        cerr << "There was a problem opening the input file!\n";
        exit(1);
    }

    std::vector<std::vector<double>> vec;
    std::string lineData;

    //load calibration targets: race, gender, sex.beh, age group, time 
    for (int i = 0; i < Parameters.nRaces; i++) {
        Parameters.calibrationTarget.push_back(vector<vector<vector<vector<double>>>>());
        for (int j = 0; j < Parameters.nGenders; j++) {
            Parameters.calibrationTarget[i].push_back(vector<vector<vector<double>>>());
            for (int k = 0; k < Parameters.nSexualBehs; k++) {
                Parameters.calibrationTarget[i][j].push_back(vector<vector<double>>());
                for (int m1 = 0; m1 < Parameters.ageGroup_size; m1++) { 
                    double d;
                    vector<double> row;
                    getline(ifile, lineData);
                    stringstream lineStream(lineData);
                    while (lineStream >> d) row.push_back(d);
                    Parameters.calibrationTarget[i][j][k].push_back(row);
                }
            }
        }
    }
}

// [[Rcpp::export]]
void loadCalibrationParameters(std::string inputPath, parameters& Parameters, int numbParallel) {
    ifstream ifile(inputPath, ios::in);

    //check file opened correctly
    if (!ifile.is_open()) {
        cerr << "There was a problem opening the input file!\n";
        exit(1);
    }

    std::vector<std::vector<double>> vec;
    std::string lineData;

    //load calibrated contact frequency: race, gender, sex.beh 
    for (int np = 0; np < numbParallel; np++) {
        Parameters.transmissionRate.push_back(vector<vector<vector<vector<double>>>>());
        for (int i = 0; i < Parameters.nRaces; i++) {
            Parameters.transmissionRate[np].push_back(vector<vector<vector<double>>>());
            for (int j = 0; j < Parameters.nGenders; j++) {
                Parameters.transmissionRate[np][i].push_back(vector<vector<double>>());
                for (int k = 0; k < Parameters.nSexualBehs; k++) {
                    Parameters.transmissionRate[np][i][j].push_back(vector<double>());
                    for (int m1 = 0; m1 < Parameters.ageGroup_size; m1++) {  //HARDCODED
                        if (np == 0) {
                            double d;
                            getline(ifile, lineData);
                            stringstream lineStream1(lineData);
                            lineStream1 >> d;
                            Parameters.transmissionRate[np][i][j][k].push_back(d);
                        }
                        else
                            Parameters.transmissionRate[np][i][j][k].push_back(Parameters.transmissionRate[0][i][j][k][m1]);
                    }
                }
            }
        }
    }
}

// [[Rcpp::export]]
void runModelDynamics(double** PopulationX, parameters& Parameters, psa_parameters* psaParameters, int iter, int np) {

    int size = Parameters.nRaces * Parameters.nGenders * Parameters.nSexualBehs * Parameters.nSexActs * Parameters.nAges * Parameters.nDiseaseStates;

    for (int year = 0; year < Parameters.timeHorizon; year++) {  //0  to 3*12+2 = 37-1 = 36
        for (int month = 0; month < 12; month++) {
            int time = year * 12 + month;

            //update starting population
            for (int index = 0; index < size; index++) PopulationX[time + 1][index] = PopulationX[time][index];

            //run RK integrator change 1 - HARDCODED
            for (int rk = 0; rk < 1; rk++) odeint(PopulationX, Parameters, psaParameters, 0, 0, 0, 0, PopulationX[time + 1], year, month, size, iter, np);
        }
    }
}

// [[Rcpp::export]]
void updatedCalibrationParameters(std::string filename, parameters& Parameters){
    std::ofstream ofile(filename);
    if (ofile.good()) {
        for (int i = 0; i < Parameters.nRaces; i++)
            for (int j = 0; j < Parameters.nGenders; j++)
                for (int k = 0; k < Parameters.nSexualBehs; k++)
                    for (int m1 = 0; m1 < Parameters.ageGroup_size; m1++)
                        ofile << Parameters.transmissionRate[0][i][j][k][m1]<< "\n";  
    }
    else {
        cerr << "There was a problem opening the output file!\n";
        exit(1);
    }
}

// [[Rcpp::export]]
void saveCalibratedIncidence(double** Population, std::string filename, parameters& Parameters) {
    
    //aggregate and save incidence
    int dim = Parameters.nRaces * Parameters.nGenders * Parameters.nSexualBehs * Parameters.ageGroup_size;

    double** Incidence = new double* [5];  //HARDCODED for 5 years
        for (int t = 0; t < 5; t++)
            Incidence[t] = new double[dim];
    
    std::ofstream ofile(filename);
    if (ofile.good()) {
        int d = 13;   //HARDCODED 13 for incidence

        for (int index = 0; index < dim; index++) {
            //for (int t = 4; t <= 4; t++){
            for (int t = 0; t < 5; t++) {

                Incidence[t][index] = 0.0;

                int m1 = index % Parameters.ageGroup_size;
                int k = ((index - m1) / Parameters.ageGroup_size) % Parameters.nSexualBehs;
                int j = ((((index - m1) / Parameters.ageGroup_size) - k) / Parameters.nSexualBehs) % Parameters.nGenders;
                int i = ((((index - m1) / Parameters.ageGroup_size) - k) / Parameters.nSexualBehs - j) / Parameters.nGenders;

                for (int month = 0; month < 12; month++) {

                    int year = Parameters.timeHorizon - 5 + t;   //HARDCODED 5

                    int time = year * 12 + month;   //HARDCODED 12 for months

                    for (int l = 0; l < Parameters.nSexActs; l++)
                    {
                        if (m1 == 0)
                            for (int m = Parameters.ageGroup_bound[0]; m < Parameters.ageGroup_bound[1]; m++) {
                                int index_13 = 13 + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));
                                if (time == 0) Incidence[t][index] = Incidence[t][index] + Population[time][index_13];
                                else Incidence[t][index] = Incidence[t][index] + (Population[time][index_13] - Population[time - 1][index_13]);
                            }
                        else
                            if (m1 == 1)
                                for (int m = Parameters.ageGroup_bound[1]; m < Parameters.ageGroup_bound[2]; m++) {
                                    int index_13 = 13 + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));
                                    if (time == 0) Incidence[t][index] = Incidence[t][index] + Population[time][index_13];
                                    else Incidence[t][index] = Incidence[t][index] + (Population[time][index_13] - Population[time - 1][index_13]);
                                }
                            else
                                if (m1 == 2)
                                    for (int m = Parameters.ageGroup_bound[2]; m < Parameters.ageGroup_bound[3]; m++) {
                                        int index_13 = 13 + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));
                                        if (time == 0) Incidence[t][index] = Incidence[t][index] + Population[time][index_13];
                                        else Incidence[t][index] = Incidence[t][index] + (Population[time][index_13] - Population[time - 1][index_13]);
                                    }
                                else
                                    for (int m = Parameters.ageGroup_bound[3]; m <= Parameters.Age_max; m++) {
                                        int index_13 = 13 + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));
                                        if (time == 0) Incidence[t][index] = Incidence[t][index] + Population[time][index_13];
                                        else Incidence[t][index] = Incidence[t][index] + (Population[time][index_13] - Population[time - 1][index_13]);
                                    }
                    }
                }
                ofile << Incidence[t][index] << "\t";
            }
            ofile << "\n";
        }
    }
    else {
        cerr << "There was a problem opening the output file!\n";
        exit(1);  
    }
}

// [[Rcpp::export]]
void saveIncidence(double** population, std::string filename, parameters& Parameters, int runTime) {

    std::ofstream ofile(filename);
    if (ofile.good()) {
        ofile.flags(std::ios::fixed);
        int j = 0;
        int k = 0;

        //write header
        ofile << " " << "," << " " << "," << " " << ",";
        for (int k = 0; k <= 2; k++) //MSW,MSM,WSM
            for (int m = 0; m < Parameters.nAges; m++)
                ofile << k << ",";
        ofile << "\n";
        ofile << "race" << "," << "state" << "," << "time" << ",";
        for (int k = 0; k <= 4; k++)
            for (int m = 0; m < Parameters.nAges; m++)
                ofile << m << ",";
        ofile << "\n";
        for (int i = 0; i < Parameters.nRaces; i++) {

            //incidence - asymptomatically infected
            int d = 12;
            for (int t = 0; t < runTime - 1; t++) {
                ofile << i << "," << 0 << "," << t << ", ";  //HARDCODED 0 for asymptomatic infections

                //MSW
                j = 0;
                k = 0;
                for (int m = 0; m < Parameters.nAges; m++) {
                    int index_0 = d + Parameters.nDiseaseStates * (m + Parameters.nAges * (0 + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));
                    int index_1 = d + Parameters.nDiseaseStates * (m + Parameters.nAges * (1 + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));
                    int index_2 = d + Parameters.nDiseaseStates * (m + Parameters.nAges * (2 + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));

                    if (t == 0) ofile << (population[t][index_0] + population[t][index_1] + population[t][index_2]) << ",";
                    else ofile << (population[t][index_0] + population[t][index_1] + population[t][index_2]) - (population[t - 1][index_0] + population[t - 1][index_1] + population[t - 1][index_2]) << ",";
                }

                //MSM
                j = 0;
                for (int m = 0; m < Parameters.nAges; m++) {
                    int index_0 = d + Parameters.nDiseaseStates * (m + Parameters.nAges * (0 + Parameters.nSexActs * (1 + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));
                    int index_1 = d + Parameters.nDiseaseStates * (m + Parameters.nAges * (1 + Parameters.nSexActs * (1 + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));
                    int index_2 = d + Parameters.nDiseaseStates * (m + Parameters.nAges * (2 + Parameters.nSexActs * (1 + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));
                    int index_3 = d + Parameters.nDiseaseStates * (m + Parameters.nAges * (0 + Parameters.nSexActs * (2 + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));
                    int index_4 = d + Parameters.nDiseaseStates * (m + Parameters.nAges * (1 + Parameters.nSexActs * (2 + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));
                    int index_5 = d + Parameters.nDiseaseStates * (m + Parameters.nAges * (2 + Parameters.nSexActs * (2 + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));

                    if (t == 0) ofile << (population[t][index_0] + population[t][index_1] + population[t][index_2] + population[t][index_3] + population[t][index_4] + population[t][index_5]) << ",";
                    else ofile << (population[t][index_0] + population[t][index_1] + population[t][index_2] + population[t][index_3] + population[t][index_4] + population[t][index_5]) - (population[t - 1][index_0] + population[t - 1][index_1] + population[t - 1][index_2] + population[t - 1][index_3] + population[t - 1][index_4] + population[t - 1][index_5]) << ",";
                }

                //WSM
                j = 1;
                k = 0;
                for (int m = 0; m < Parameters.nAges; m++) {
                    int index_0 = d + Parameters.nDiseaseStates * (m + Parameters.nAges * (0 + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));
                    int index_1 = d + Parameters.nDiseaseStates * (m + Parameters.nAges * (1 + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));
                    int index_2 = d + Parameters.nDiseaseStates * (m + Parameters.nAges * (2 + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));

                    if (t == 0) ofile << (population[t][index_0] + population[t][index_1] + population[t][index_2]) << ",";
                    else ofile << (population[t][index_0] + population[t][index_1] + population[t][index_2]) - (population[t - 1][index_0] + population[t - 1][index_1] + population[t - 1][index_2]) << ",";
                }
                ofile << "\n";
            }

            //incidence - symptomatic infections
            d = 13;
            for (int t = 0; t < runTime - 1; t++) {
                ofile << i << "," << 0 << "," << t << ", ";  //HARDCODED 0 for asymptomatic infections

                //MSW
                j = 0;
                k = 0;
                for (int m = 0; m < Parameters.nAges; m++) {
                    int index_0 = d + Parameters.nDiseaseStates * (m + Parameters.nAges * (0 + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));
                    int index_1 = d + Parameters.nDiseaseStates * (m + Parameters.nAges * (1 + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));
                    int index_2 = d + Parameters.nDiseaseStates * (m + Parameters.nAges * (2 + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));

                    if (t == 0) ofile << (population[t][index_0] + population[t][index_1] + population[t][index_2]) << ",";
                    else ofile << (population[t][index_0] + population[t][index_1] + population[t][index_2]) - (population[t - 1][index_0] + population[t - 1][index_1] + population[t - 1][index_2]) << ",";
                }

                //MSM
                j = 0;
                for (int m = 0; m < Parameters.nAges; m++) {
                    int index_0 = d + Parameters.nDiseaseStates * (m + Parameters.nAges * (0 + Parameters.nSexActs * (1 + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));
                    int index_1 = d + Parameters.nDiseaseStates * (m + Parameters.nAges * (1 + Parameters.nSexActs * (1 + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));
                    int index_2 = d + Parameters.nDiseaseStates * (m + Parameters.nAges * (2 + Parameters.nSexActs * (1 + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));
                    int index_3 = d + Parameters.nDiseaseStates * (m + Parameters.nAges * (0 + Parameters.nSexActs * (2 + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));
                    int index_4 = d + Parameters.nDiseaseStates * (m + Parameters.nAges * (1 + Parameters.nSexActs * (2 + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));
                    int index_5 = d + Parameters.nDiseaseStates * (m + Parameters.nAges * (2 + Parameters.nSexActs * (2 + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));

                    if (t == 0) ofile << (population[t][index_0] + population[t][index_1] + population[t][index_2] + population[t][index_3] + population[t][index_4] + population[t][index_5]) << ",";
                    else ofile << (population[t][index_0] + population[t][index_1] + population[t][index_2] + population[t][index_3] + population[t][index_4] + population[t][index_5]) - (population[t - 1][index_0] + population[t - 1][index_1] + population[t - 1][index_2] + population[t - 1][index_3] + population[t - 1][index_4] + population[t - 1][index_5]) << ",";
                }

                //WSM
                j = 1;
                k = 0;
                for (int m = 0; m < Parameters.nAges; m++) {
                    int index_0 = d + Parameters.nDiseaseStates * (m + Parameters.nAges * (0 + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));
                    int index_1 = d + Parameters.nDiseaseStates * (m + Parameters.nAges * (1 + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));
                    int index_2 = d + Parameters.nDiseaseStates * (m + Parameters.nAges * (2 + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));

                    if (t == 0) ofile << (population[t][index_0] + population[t][index_1] + population[t][index_2]) << ",";
                    else ofile << (population[t][index_0] + population[t][index_1] + population[t][index_2]) - (population[t - 1][index_0] + population[t - 1][index_1] + population[t - 1][index_2]) << ",";
                }

                ofile << "\n";
            }
        }
    }
    else {
        cerr << "There was a problem opening the output file!\n";
        exit(1);  //exit or do additional error checking
    }
}

// [[Rcpp::export]]
void saveTrajectories(double** population, std::string filename, parameters& Parameters, int runTime) {

    std::ofstream ofile(filename);
    if (ofile.good()) {
        ofile.flags(std::ios::fixed);
        for (int d = 0; d < Parameters.nDiseaseStates - 2; d++) //to remove incidence (symptomatic and asymptomatic)
            for (int i = 0; i < Parameters.nRaces; i++)
                for (int j = 0; j < Parameters.nGenders; j++)
                    for (int k = 0; k < Parameters.nSexualBehs; k++)
                        for (int l = 0; l < Parameters.nSexActs; l++)
                            for (int m = 0; m < Parameters.nAges; m++)
                            {
                                int index_d = d + Parameters.nDiseaseStates * (m + Parameters.nAges * (l + Parameters.nSexActs * (k + Parameters.nSexualBehs * (j + Parameters.nGenders * i))));
                                for (int t = 0; t < runTime - 2; t++) ofile << population[t][index_d] << "\t";
                                ofile << "\n";
                            }
    }
    else {
        cerr << "There was a problem opening the output file!\n";
        exit(1);
    }
}

// [[Rcpp::export]]
int runmodel(Rcpp::List inputs)
{
  int NumbSim = 2;
  
  parameters Parameters;
  
  Parameters.nRaces = 3;
  Parameters.nGenders = 2;
  Parameters.nSexualBehs = 3;
  Parameters.nSexActs = 3;
  Parameters.nAges = 101;
  Parameters.nDiseaseStates = 14;  //12 disease related stated plus two more state for incidence to be corrected
  Parameters.timeHorizon = 7;      // atoi(argv[10]);
  
  std::string Path = "./Inputs/";
  
  psa_parameters psaParameters[1000];
  
  int maxRunTime = Parameters.timeHorizon * 12 + 2;  //time in time units
  
  int NumbParallel = 1; //HARDCODED placeholder for multithreading 
  
  double*** PopulationX = new double** [NumbParallel];
  
  for (int np = 0; np < NumbParallel; np++) {
    PopulationX[np] = new double* [maxRunTime];
    
    for (int t = 0; t < maxRunTime; t++)
      PopulationX[np][t] = new double[Parameters.nRaces * Parameters.nGenders * Parameters.nSexualBehs * Parameters.nSexActs * Parameters.nAges * Parameters.nDiseaseStates];
  }
  
  double initialInfection = 0.1;
  
  loadInitialPopulation(Path + "Initial population.txt", PopulationX, Parameters, NumbParallel, initialInfection);
  loadDemographics(Path + "Demographic parameters.txt", Parameters);
  loadParameters(Path + "Parameters.txt", Parameters, psaParameters, NumbSim);
  loadCalibrationParameters(Path + "Calibration parameters.txt", Parameters, NumbParallel);
  
  //run model dynamics
  runModelDynamics(PopulationX[0], Parameters, psaParameters, 0,0);
  
  //save 
  updatedCalibrationParameters(Path + "Calibration parameters.txt", Parameters);
  saveCalibratedIncidence(PopulationX[0], Path + "Calibrated incidence.txt", Parameters);
  
  //save outputs
  saveIncidence(PopulationX[0], Path + "Incidence.csv", Parameters, maxRunTime);
  saveTrajectories(PopulationX[0], Path + "Trajectories.txt", Parameters, maxRunTime);
  
  //release memory
  for (int np = 0; np < NumbParallel; np++) {
    for (int j = 0; j < maxRunTime; j++) {
      delete[] PopulationX[np][j];
    }
    delete[] PopulationX[np];
  }
  delete[] PopulationX;
  
  return 0;
}

