// propagation.h
#ifndef PROPAGATION_H
#define PROPAGATION_H

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

//Particle propagation for any system.
void propagation(vector< vector<uint64_t> > result_cell[7], vector< vector<uint64_t> > cell[7], vector< vector<uint64_t> > nsbit);

#endif // PROPAGATION_H