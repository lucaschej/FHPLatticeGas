// collisionFHP_III.h
#ifndef COLLISION_FHP_III_H
#define COLLISION_FHP_III_H

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

void collisionFHP_III(vector< vector<uint64_t> > cell[7],  vector< vector<uint64_t> > result_cell[7], vector< vector<uint64_t> > nsbit, mt19937 generator);

#endif // COLLISION_FHP_III_H