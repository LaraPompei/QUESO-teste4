#ifndef COMPUTE_H
#define COMPUTE_H

#include <queso/Environment.h>
#include <cmath>
#include <queso/GslMatrix.h>
#include <queso/GenericVectorFunction.h>
#include <queso/GaussianVectorRV.h>
#include <queso/LogNormalVectorRV.h>
#include <queso/GenericVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/StatisticalForwardProblem.h>
#include <queso/VectorSet.h>

#include <sys/time.h>
#include <random>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include <cmath>

using namespace QUESO;
using namespace std;

void compute(const FullEnvironment& env);
void filling_matrix(double* t, double* values_a, double** data);
void save_data(double* model, double** baseModel, double* data, double* values_of_a);
#endif
