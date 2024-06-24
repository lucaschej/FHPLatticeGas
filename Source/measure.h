// measure.h
#ifndef MEASURE_H
#define MEASURE_H_H

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

void measure(vector< vector<uint64_t> > cell[7], double ci[7][2], vector< vector<double> > Ni[7], vector< vector<double> >&  rho, vector< vector<double> >  u[2], double mean_data[NMEANDATA]);

#endif // MEASURE_H_H