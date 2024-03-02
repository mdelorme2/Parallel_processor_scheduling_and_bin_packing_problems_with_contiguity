#include "model.h"

int FLOW(const Instance& inst, Solution& sol, const int& DB, const int& TL) {
	
	// local variable
	vector<bool> isWActive(inst.W,false);
	
    // create a model
    GRBEnv env = GRBEnv();              	
    GRBModel model = GRBModel(env);         

    // declaration of the variables for the model
    vector<GRBVar> x (inst.arcs.size());
	
    // initizalization of the variables for the model
    for (int k = 0; k < inst.arcs.size(); k++){
		if(inst.RCs[k] + inst.contR - EPSILON > DB) continue; 
		isWActive[inst.arcs[k][0]] = true;
		isWActive[inst.arcs[k][1]] = true;
        if (inst.arcs[k][3] >= 0) x[k] = model.addVar(0, inst.items[inst.arcs[k][3]][2], 0, GRB_INTEGER);
		else x[k] = model.addVar(0, GRB_INFINITY, 0, GRB_INTEGER);
 	}
	
    // create linear expressions
    vector<GRBLinExpr> assigned(inst.n, 0); 
    vector<GRBLinExpr> cIn(inst.W+1, 0);  
	vector<GRBLinExpr> cOut(inst.W+1, 0); 	
    for (int k = 0; k < inst.arcs.size(); k++){
		if(inst.RCs[k] + inst.contR - EPSILON > DB) continue; 
		cIn[inst.arcs[k][1]] += inst.arcs[k][2] * x[k];
        cOut[inst.arcs[k][0]] += inst.arcs[k][2] * x[k];
		if(inst.arcs[k][3] >= 0) assigned[inst.arcs[k][3]] += x[k];
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
	model.addConstr(cOut[0] <= DB);
	
    // change some settings
    model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
    model.getEnv().set(GRB_IntParam_Threads, 1);
    model.getEnv().set(GRB_DoubleParam_TimeLimit, TL);

    // find the optimal solution	
    model.optimize();

    // store the results in a Solution object
    sol.Nvar = model.get(GRB_IntAttr_NumVars);       
    sol.Nconstr = model.get(GRB_IntAttr_NumConstrs); 
    sol.Ncoeff = model.get(GRB_IntAttr_NumNZs);      	
	sol.opt = 0;
		
	// if a solution has been found
    if (model.get(GRB_IntAttr_SolCount) >= 1) { 		
		sol.LB += DB;
        sol.UB += DB;
		sol.opt = 1;
		
        // get bin for each item
        for (int k = 0; k < inst.arcs.size(); k++){   
			if(inst.RCs[k] + inst.contR - EPSILON > DB) continue; 		
			if(inst.arcs[k][3] >= 0){
				for (int l = 0; l < ceil(x[k].get(GRB_DoubleAttr_X) - EPSILON); l++)  
					sol.assignedBin[inst.arcs[k][3]].push_back(inst.arcs[k][0]);
            }
        }
		return 0;
    }
	else{
		if(model.get(GRB_IntAttr_Status) == 9){
			sol.LB += DB;
			return -2;
		}
		return -1;
    }
}

void FLOWLP(Instance& inst) {
	
	vector<bool> isWActive(inst.W,false);
	
    // create a model
    GRBEnv env = GRBEnv();              	
    GRBModel model = GRBModel(env);         

    // declaration of the variables for the model
    vector<GRBVar> x (inst.arcs.size());
	
    // initizalization of the variables for the model
    for (int k = 0; k < inst.arcs.size(); k++){
		isWActive[inst.arcs[k][0]] = true;
		isWActive[inst.arcs[k][1]] = true;
        if (inst.arcs[k][3] >= 0) x[k] = model.addVar(0, inst.items[inst.arcs[k][3]][2], 0, GRB_CONTINUOUS);
		else x[k] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
 	}
	
    // create linear expressions
    vector<GRBLinExpr> assigned(inst.n, 0); 
    vector<GRBLinExpr> cIn(inst.W+1, 0);  
	vector<GRBLinExpr> cOut(inst.W+1, 0); 	
    for (int k = 0; k < inst.arcs.size(); k++){
		cIn[inst.arcs[k][1]] += inst.arcs[k][2] * x[k];
        cOut[inst.arcs[k][0]] += inst.arcs[k][2] * x[k];
		if(inst.arcs[k][3] >= 0) assigned[inst.arcs[k][3]] += x[k];
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

	// get reduced costs
	for (int k = 0; k < inst.arcs.size(); k++){    
		inst.RCs[k] = 0;  
	}
		
	inst.contR = model.get(GRB_DoubleAttr_ObjVal);
}


