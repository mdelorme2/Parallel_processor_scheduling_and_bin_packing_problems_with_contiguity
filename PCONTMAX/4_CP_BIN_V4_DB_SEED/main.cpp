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
	BMLP(inst);
	int DB = ceil(inst.contR-EPSILON);
	int TL = 3600 - (getCPUTime() - start);
	
	// Compute normal patterns for DB
	vector<bool> DBNP(DB*2,false);
	for (int j = 0; j < inst.n; j++) {
		if(inst.items[j][2] == 0) continue; 
		DBNP[0] = true;
		for (int jp = 0; jp<inst.n;jp++){
			int lim = inst.items[jp][2];
			for(int l = 0; l<lim;l++){
				for (int i = DBNP.size()-inst.items[jp][1];i >= 0; i--){
					if (DBNP[i]) DBNP[i+inst.items[jp][1]] = true;
				}
			}
		}
	}

	while(TL >= 1 && CPM(inst,sol,DB,TL,seed) == -1){
		TL = 3600 - (getCPUTime() - start);
		DB++;
		while(!DBNP[DB]) DB++; 
	}
	sol.contR = inst.contR;
	sol.timeT = getCPUTime() - start;	
	printInfo(pathAndFileout, sol, filein);
	printSInfo(sol, inst);
}

