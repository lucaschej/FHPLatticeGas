// promediate.h
#ifndef PROMEDIATE_H
#define PROMEDIATE_H

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

//Measures the most important values of the distribution
//void measure(vector< vector<uint64_t> > cell[7], double ci[7][2], vector< vector<double> > Ni[7], vector< vector<double> >&  rho, vector< vector<double> >  u[2], double mean_data[NMEANDATA]);
void promediate(double mean_data[NMEANDATA], vector<double> ni_to_write[7], vector<double> u_to_write[2], double& rho_eq, double u_eq[2], double Ni_eq[7], int fraction_its, int n_its);

#endif // PROMEDIATE_H