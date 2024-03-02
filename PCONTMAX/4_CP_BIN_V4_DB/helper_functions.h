#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H

using namespace std;
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <sstream> 
#include <time.h>
#include <sys/time.h>
#include <ilcp/cp.h>

const double EPSILON = 0.00001; // small constant
const int M = 1000000;			// big constant

struct Instance
{
	int n; 						// number of items
	int W, Wo; 					// number of bins
	vector<vector<int> > items;	// items
	double contR;
	void print();
};

struct Solution
{
	int opt, LB, UB, Nvar, Nconstr, Ncoeff;
	double timeP, timeT, contR;
    vector<vector<int>> assignedBin;
};

double getCPUTime();
bool sortItems(const vector<int>& v1, const vector<int>& v2);
Instance readInstance(const string& filename, Solution& sol);
int CP(const Instance& inst, const int& heightL, vector<vector<int> >& assignedBin);
void printInfo(const string& pathAndFileout, const Solution& sol, const string& filein);
void printSInfo(const Solution& sol, const Instance& inst);

#endif 