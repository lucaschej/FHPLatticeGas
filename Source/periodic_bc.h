// periodic_bc.h
#ifndef PERIODIC_BC_H
#define PERIODIC_BC_H

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

//Periodic boundary conditions
int periodic_bc(int k, int maxk);
int periodic_bc(int b, int maxb, int& x);

#endif // PERIODIC_BC_H