#ifndef LB_H
#define LB_H

using namespace std;
#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <algorithm>
#include <math.h> 
#include "gurobi_c++.h"
#include "helper_functions.h"

void LB(const Instance& inst, Solution& sol);

#endif 