#include "main.h"

int main(int argc, char **argv){          

	// Read input and output paths
	double start = getCPUTime(); 	
	string path = argv[1];	
	string filein = argv[2];
	string pathAndFileout = argv[3];
	int seed = atoi(argv[4]);
	Solution sol; sol.LB = 0; sol.UB = 0; sol.Nvar = 0; sol.Nconstr = 0; sol.Ncoeff = 0; sol.contR = 0;

    // initialize the input variables from a file
    Instance inst = readInstance(path + filein, sol);
	sol.timeP = getCPUTime() - start;
	inst.print();

    // find the CPM solution
    CPM(inst,sol,seed);
	sol.timeT = getCPUTime() - start;	
	printInfo(pathAndFileout, sol, filein);
	printSInfo(sol, inst);
}

