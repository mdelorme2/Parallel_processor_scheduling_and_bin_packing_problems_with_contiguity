#include "model.h"

void FLOW(const Instance& inst, Solution& sol) {
	
	// local variable
	vector<bool> isWActive(inst.W,false);
	vector<int> isWActiveV;

	// compute normal patterms
	vector<vector<bool> > NPs (inst.n,vector<bool>(inst.W,false));
	for (int j = 0; j < inst.n; j++) {
		if(inst.items[j][2] == 0) continue; 
		NPs[j][0] = true;
		for (int jp = 0; jp<inst.n;jp++){
			int lim = inst.items[jp][2];
			if (jp == j) lim--; 
			for(int l = 0; l<lim;l++){
				for (int i = inst.W-inst.items[jp][0] - 1;i >= 0; i--){
					if (NPs[j][i]) NPs[j][i+inst.items[jp][0]] = true;
				}
			}
		}
	}
	
	// create item arcs
	vector<vector<int> > arcs;
	for (int j = 0; j < inst.n; j++) {
		if(inst.items[j][2] == 0) continue; 
		for (int i = 0; i <= inst.W - inst.items[j][0]; i++){ 
			if(NPs[j][i]){ 
				arcs.push_back({i,i+inst.items[j][0],inst.items[j][1],j});
				isWActive[i] = true;
				isWActive[i+inst.items[j][0]] = true;
			}
		}
	}
	cout << "There are " << arcs.size() << " arcs, loss arcs excluded" << endl;

	// create loss arcs
	for(int i = 0; i<inst.W;i++)
		if(isWActive[i]) isWActiveV.push_back(i);
	isWActiveV.push_back(inst.W);	
	for (int i = 0; i < isWActiveV.size()-1; i++)
		arcs.push_back({isWActiveV[i],isWActiveV[i+1],1,-1});
	cout << "There are " << arcs.size() << " arcs, loss arcs included" << endl;
	
    // create a model
    GRBEnv env = GRBEnv();              	
    GRBModel model = GRBModel(env);         

    // declaration of the variables for the model
    vector<GRBVar> x (arcs.size());
	
    // initizalization of the variables for the model
    for (int k = 0; k < arcs.size(); k++){
        if (arcs[k][3] >= 0) x[k] = model.addVar(0, inst.items[arcs[k][3]][2], 0, GRB_INTEGER);
		else x[k] = model.addVar(0, GRB_INFINITY, 0, GRB_INTEGER);
 	}
	
    // create linear expressions
    vector<GRBLinExpr> assigned(inst.n, 0); 
    vector<GRBLinExpr> cIn(inst.W+1, 0);  
	vector<GRBLinExpr> cOut(inst.W+1, 0); 	
    for (int k = 0; k < arcs.size(); k++){
		cIn[arcs[k][1]] += arcs[k][2] * x[k];
        cOut[arcs[k][0]] += arcs[k][2] * x[k];
		if(arcs[k][3] >= 0) assigned[arcs[k][3]] += x[k];
    }
    
    model.update();

    // create assignment constraints
    for (int j = 0; j < inst.n; j++){
		if(inst.items[j][2] == 0) continue; 		
		model.addConstr(assigned[j] == inst.items[j][2]); 					
	}
	
    // create flow conservation constraints 
    for (int i = 1; i < inst.W; i++){
		if (isWActive[i])		
			model.addConstr(cIn[i] == cOut[i]);
	}
	model.addConstr(cOut[0] == cIn[inst.W]);
	
	// set the objective: minimize z
    model.setObjective(cOut[0], GRB_MINIMIZE);

    // change some settings
    model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
    model.getEnv().set(GRB_IntParam_Threads, 1);
    model.getEnv().set(GRB_DoubleParam_TimeLimit, 3600);

    // find the optimal solution	
    model.optimize();

    // store the results in a Solution object
    sol.Nvar = model.get(GRB_IntAttr_NumVars);       
    sol.Nconstr = model.get(GRB_IntAttr_NumConstrs); 
    sol.Ncoeff = model.get(GRB_IntAttr_NumNZs);      	
	sol.opt = 0;
	sol.LB += ceil(model.get(GRB_DoubleAttr_ObjBound) - EPSILON);
		
	// if a solution has been found
    if (model.get(GRB_IntAttr_SolCount) >= 1) { 		
        sol.UB += ceil(model.get(GRB_DoubleAttr_ObjVal) - EPSILON);
		if(sol.LB == sol.UB) sol.opt = 1;
		
        // get bin for each item
        for (int k = 0; k < arcs.size(); k++){    
			if(arcs[k][3] >= 0){
				for (int l = 0; l < ceil(x[k].get(GRB_DoubleAttr_X) - EPSILON); l++)  
					sol.assignedBin[arcs[k][3]].push_back(arcs[k][0]);
            }
        }
    }
	/*GRBModel model2 = model.relax();
	model2.optimize();
	sol.contR += model2.get(GRB_DoubleAttr_ObjVal);*/
}

