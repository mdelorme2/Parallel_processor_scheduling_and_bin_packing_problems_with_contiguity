#include "helper_functions.h"

void Instance::print(){
	cout << "n = " << n << " " << "W = " << W << endl;
	for(int j=0; j<n;j++)
		cout << "Item " << j << ": " << items[j][0] << " " << items[j][1] << " " << items[j][2]  << endl;
}

double getCPUTime(){
	return (double)clock() / CLOCKS_PER_SEC;
}

bool sortItems(const vector<int>& v1, const vector<int>& v2) {
    return M * v1[0] +  v1[1] >  M * v2[0] + v2[1];
}

Instance readInstance(const string& filename, Solution& sol) {
    // define variables
	Instance inst;

    // open the file
    ifstream file(filename); 

    // read the file
    if (file.is_open()) { //if the file is open
        string line;

        // first line contains number of items
        getline(file, line, '\n');
        inst.n = stoi(line); 
		sol.assignedBin.resize(inst.n);

        // reshape the array that will hold the items
        inst.items.resize(inst.n, vector<int>(3, 0));

        // second line contains the number of bins
        getline(file, line, '\n');
        inst.W = stoi(line); inst.Wo = inst.W;

        // the remaining lines contain the width, the height, and the demand of the items
        for (int j = 0; j < inst.n; j++) {
            getline(file, line, ' ');
            inst.items[j][0] = stoi(line);
            getline(file, line, ' ');
            inst.items[j][1] = stoi(line);
            getline(file, line, '\n');
            inst.items[j][2] = stoi(line);
        }

        // close the file
        file.close(); 
    }
    else {
        // if the file cannot be opened: print error and return default Inst
        cout << "Unable to open file"; 
    }

	// Preprocessing 1 -- prepack items	
	/*sort(inst.items.begin(), inst.items.end(), sortItems);
	Instance inst2 = inst;
	int sumAreaL = 0; int sumAreaR = 0;
	int heightL = 0; int heightR = 0;
	int idx2 = inst.n;
	for(int j = 0; j < inst.n; j++) inst2.items[j][2] = 0;
	
	for (int j = 0; j < inst.n; j++) {
		if (inst.items[j][0] <= inst.W / 2.0) break;
		
		// The left
		heightL += inst.items[j][1] * inst.items[j][2];
		sumAreaL += inst.items[j][0] * inst.items[j][1] * inst.items[j][2];
		inst2.items[j][2] = inst.items[j][2];
		
		// The right
		while (idx2 >= j+2 && inst.items[idx2-1][0] + inst.items[j][0] <= inst.W){
			idx2--;
			heightR = max(heightR,inst.items[idx2][1]);
			sumAreaR += inst.items[idx2][0] * inst.items[idx2][1] * inst.items[idx2][2];
			inst2.items[idx2][2] = inst.items[idx2][2];			
		}

		// Feasability test
		cout << "Test " << j << " " << idx2;
		if(sumAreaL + sumAreaR <= inst.W * heightL && heightL >= heightR){
			if (CP(inst2,heightL,sol.assignedBin) == 0){
				sol.LB += heightL; sol.UB += heightL; sol.contR += heightL; 
				sumAreaL = 0; sumAreaR = 0; heightL = 0; heightR =0;
				for(int k = 0; k < inst.n; k++) inst2.items[k][2] = 0;		
				for(int k = 0; k <= j; k++) inst.items[k][2] = 0;
				for(int k = idx2; k < inst.n; k++) inst.items[k][2] = 0;
				cout << " prepacked!" << endl;
			}
			else{
				cout << " not prepacked" << endl;
			}
		}
		else{
			cout << " no hope" << endl;
		}
	}

	sort(inst.items.begin(), inst.items.end(), sortItems);
	
	// Preprocessing 2 -- decrease W
	vector<bool> isR(inst.W+1,false); isR[0] = true;
	for (int j = 0; j < inst.n; j++) {
		for(int l = 0; l<inst.items[j][2];l++){
			for (int i = inst.W-inst.items[j][0];i >= 0; i--){
				if (isR[i]) isR[i+inst.items[j][0]] = true;
			}
		}
	}
	cout << "W from " << inst.W << " to ";
	while (!isR[inst.W] && inst.W >= 2) inst.W--;
	cout << inst.W << endl;

	// Preprocessing 3 -- inrease w
	for (int j = 0; j < inst.n; j++) {
		if(inst.items[j][2] == 0) continue;
		isR.resize(0); isR.resize(inst.W - inst.items[j][0] +1,false); isR[0] = true;
		for (int jp = 0; jp<inst.n;jp++){
			int lim = inst.items[jp][2];
			if (jp == j) lim--; 
			for(int l = 0; l<lim;l++){
				for (int i = inst.W- inst.items[j][0]-inst.items[jp][0];i >= 0; i--){
					if (isR[i]) isR[i+inst.items[jp][0]] = true;
				}
			}
		}
		int tw = inst.W - inst.items[j][0];
		cout << "Item " << j << " from " << inst.items[j][0] << " to ";
		while(!isR[tw]) tw--;
		if(2 * inst.items[j][0] > inst.W) inst.items[j][0] += (inst.W - inst.items[j][0] - tw);
		else inst.items[j][0] += (inst.W - inst.items[j][0] - tw)/inst.items[j][2];
		cout << inst.items[j][0] << endl;
	}*/
	
    return inst;
}

int CP(const Instance& inst, const int& heightL, vector<vector<int> >& assignedBin) {
	
    // create a model
	IloEnv env;     
	IloModel model(env);
	
	// from CSP to BPP
	vector<vector<int> > items2;
	int n2 = 0;
    for (int j = 0; j < inst.n; j++) { 
		for (int k = 0; k < inst.items[j][2]; k++){
			items2.push_back({inst.items[j][0],inst.items[j][1],j});
			n2++;
		}			
	}
	
    // declaration of the variables for the model
	IloIntervalVarArray x(env, n2);
	IloCumulFunctionExpr z (env);
	
    // initizalization of the variables for the model
    for (int j = 0; j < n2; j++) { 
		x[j] = IloIntervalVar(env); 
		x[j].setStartMin(0);
		if(items2[j][0] > inst.W / 2) x[j].setStartMax(0);
		else x[j].setStartMax(inst.W - items2[j][0]);
		if(j > 0 && items2[j][2] == items2[j-1][2])  model.add(IloStartBeforeStart(env, x[j-1], x[j]));
	    x[j].setSizeMin(items2[j][0]);
		x[j].setSizeMax(items2[j][0]);
		z += IloPulse(x[j], items2[j][1]);
    }

    // set the objective: minimize z
	model.add(z <= heightL);

    // change some settings
    IloCP cp(model);
    cp.setParameter(IloCP::TimeLimit, 3600);
    cp.setParameter(IloCP::Workers, 1);
	cp.setParameter(IloCP::LogPeriod, 500000);
	cp.setParameter(IloCP::LogVerbosity, IloCP::Quiet); 
	cp.solve();	

    // find the optimal solution    
	if(cp.solve()){
		// get bin for each item
        for (int j = 0; j < n2; j++){                                 
			// cout << "Task "<< j << " starts at time " << cp.getStart(x[j]) << endl;
			assignedBin[items2[j][2]].push_back(cp.getStart(x[j]));
		}
		return 0;
	}
	
    return -1;
}

void printInfo(const string& pathAndFileout, const Solution& sol, const string& filein){
	string nameFile = pathAndFileout;
	std::ofstream file(nameFile.c_str(), std::ios::out | std::ios::app);
	file << filein << "\t"  << sol.opt << "\t" << sol.timeT << "\t" << sol.timeP << "\t" << sol.LB << "\t" << sol.UB << "\t" << sol.contR  << "\t" << sol.Nvar << "\t" << sol.Nconstr << "\t" << sol.Ncoeff  << endl;
	file.close();
}

void printSInfo(const Solution& sol, const Instance& inst){
	vector<vector<int> > bins(inst.Wo);
	vector<int>	load(inst.Wo,0);
	for (int j = 0; j < sol.assignedBin.size();j++){
		for (int i = 0; i < sol.assignedBin[j].size();i++){
			for (int ip = sol.assignedBin[j][i]; ip < sol.assignedBin[j][i] + inst.items[j][0];ip++){
				load[ip] += inst.items[j][1];
				bins[ip].push_back(j);
			}
		}
	}
	for (int i = 0; i < inst.Wo;i++){
		cout << "Bin " << i << "\t" << load[i] << "\t [";
		for (int j = 0; j <  bins[i].size(); j++) 
			cout << bins[i][j] << " ";
		cout << "]" << endl;
	}
}
