// get_results.h

#ifndef GET_RESULTS_H
#define GET_RESULTS_H

#include <vector>
#include <random>
#include "config.h"
#include <cstdint>
#include <cmath>

#include<iostream>
#include<ctime>
#include<fstream>
#include<string>
#include<chrono>

using namespace std;

//Get the results for g(d), nu(d), etc for the models
void get_results_FHPI(double rho, double& d, double& cs, double& gd, double& nud, double& etad, double& re);
void get_results_FHPII(double rho, double& d, double& cs, double& gd, double& nud, double& etad, double& re);
void get_results_FHPIII(double rho, double& d, double& cs, double& gd, double& nud, double& etad, double& re);

#endif // GET_RESULTS_H