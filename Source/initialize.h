// initialize.h
#ifndef INITIALIZE_H
#define INITIALIZE_H

#include <vector>
#include <random>
#include <cstdint>

using namespace std;

int initialize(vector< vector<uint64_t> > cell[7],  vector< vector<uint64_t> > result_cell[7], vector< vector<uint64_t> >& nsbit, mt19937 generator, double prob);

#endif // INITIALIZE_H
