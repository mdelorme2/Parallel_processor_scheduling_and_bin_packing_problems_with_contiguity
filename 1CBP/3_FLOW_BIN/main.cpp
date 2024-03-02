#include "main.h"

int main(int argc, char **argv){          

	// Read input and output paths
	double start = getCPUTime(); 	
	string path = argv[1];	
	string filein = argv[2];
	string pathAndFileout = argv[3];
	Solution sol; sol.LB = 0; sol.UB = 0; 
	
    // initialize the input variables from a file
    Instance inst = readInstance(path + filein, sol);
	sol.timeP = getCPUTime() - start;
	inst.print();

    // find the FLOW solution
    FLOW(inst,sol);
	sol.timeT = getCPUTime() - start;	
	printInfo(pathAndFileout, sol, filein);
	printSInfo(sol, inst);
}

