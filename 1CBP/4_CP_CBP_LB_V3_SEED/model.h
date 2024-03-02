#ifndef MODEL_H
#define MODEL_H

using namespace std;
#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <algorithm>
#include <math.h> 
#include "gurobi_c++.h"
#include "helper_functions.h"

void CPM(const Instance& inst, Solution& sol, const int& seed);
int CPMUB(const Instance& inst, const int& seed);
int LB(const Instance& inst);
#endif 