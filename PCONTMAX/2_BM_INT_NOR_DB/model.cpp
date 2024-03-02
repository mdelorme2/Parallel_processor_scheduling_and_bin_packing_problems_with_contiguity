#include "model.h"

int BM(const Instance& inst, Solution& sol, const int& DB, const int& TL) {
	cout << "Test " << DB << "(" << inst.contR << ") with " << TL << "s left" << endl;
	
	// local variable
	vector<bool> isWActive(inst.W,false);
	
    // create a model
    GRBEnv env = GRBEnv();              	
    GRBModel model = GRBModel(env);         

    // declaration of the variables for the model
    vector<vector<GRBVar> > x;			
    x.resize(inst.W, vector<GRBVar>(inst.n)); 
	
    // initizalization of the variables for the model
    for (int j = 0; j < inst.n; j++) { 
		if(inst.items[j][2] == 0) continue; 
		for (int i = 0; i <= inst.W - inst.items[j][0]; i++) {       	     
			if(inst.NPs[j][i] && inst.RCs[j][i] + inst.contR - EPSILON <= DB){ 
				isWActive[i] = true;
				x[i][j] = model.addVar(0, inst.items[j][2], 0, GRB_INTEGER);
			}
        }
    }
    model.update();

    // create linear expressions
    vector<GRBLinExpr> assigned(inst.n, 0);  
    vector<GRBLinExpr> height(inst.W, 0);   
    for (int j = 0; j < inst.n; j++) { 	
		if(inst.items[j][2] == 0) continue; 
        for (int i = 0; i <= inst.W - inst.items[j][0]; i++) { 
			if(inst.NPs[j][i] && inst.RCs[j][i] + inst.contR - EPSILON <= DB){
				assigned[j] += x[i][j];	
				for (int ip = i; ip < i+inst.items[j][0]; ip++) { 	
					height[ip] 	+= inst.items[j][1] * x[i][j];
				}
			}
        }
    }
    model.update();

    // create assignment constraints
    for (int j = 0; j < inst.n; j++){     
		if(inst.items[j][2] == 0) continue; 	
		model.addConstr(assigned[j] == inst.items[j][2]); 					
	}

    // create height constraints 
    for (int i = 0; i < inst.W; i++){
		if (isWActive[i])		
			model.addConstr(height[i] <= DB);
	}
	
    // change some settings
    model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
    model.getEnv().set(GRB_IntParam_Threads, 1);
	model.getEnv().set(GRB_IntParam_MinRelNodes, 0);
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
        for (int j = 0; j < inst.n; j++) {    
			if(inst.items[j][2] == 0) continue; 		
			for (int i = 0; i <= inst.W - inst.items[j][0]; i++) { 
				if(inst.NPs[j][i] && inst.RCs[j][i] + inst.contR - EPSILON <= DB){  
					for (int k = 0; k < ceil(x[i][j].get(GRB_DoubleAttr_X) - EPSILON); k++) 
						sol.assignedBin[j].push_back(i);
				}
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

void BMLP(Instance& inst) {
	
	// local variable
	vector<bool> isWActive(inst.W,false);

    // create a model
    GRBEnv env = GRBEnv();              	
    GRBModel model = GRBModel(env);         

    // declaration of the variables for the model
    vector<vector<GRBVar> > x;			
    x.resize(inst.W, vector<GRBVar>(inst.n)); 
	GRBVar z = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
	GRBLinExpr obj = z;
	
    // initizalization of the variables for the model
    for (int j = 0; j < inst.n; j++) { 
		if(inst.items[j][2] == 0) continue; 
		for (int i = 0; i <= inst.W - inst.items[j][0]; i++) {       	     
			if(inst.NPs[j][i]){ 
				isWActive[i] = true;
				x[i][j] = model.addVar(0, inst.items[j][2], 0, GRB_CONTINUOUS);
			}
        }
    }
    model.update();

    // create linear expressions
    vector<GRBLinExpr> assigned(inst.n, 0);  
    vector<GRBLinExpr> height(inst.W, 0);   
    for (int j = 0; j < inst.n; j++) { 	
		if(inst.items[j][2] == 0) continue; 
        for (int i = 0; i <= inst.W - inst.items[j][0]; i++) { 
			if(inst.NPs[j][i]){
				assigned[j] += x[i][j];	
				for (int ip = i; ip < i+inst.items[j][0]; ip++) { 	
					height[ip] 	+= inst.items[j][1] * x[i][j];
				}
			}
        }
    }
    model.update();

    // create assignment constraints
    for (int j = 0; j < inst.n; j++){     
		if(inst.items[j][2] == 0) continue; 	
		model.addConstr(assigned[j] == inst.items[j][2]); 					
	}

    // create height constraints 
    for (int i = 0; i < inst.W; i++){
		if (isWActive[i])		
			model.addConstr(height[i] <= z);
	}
	
    // set the objective: minimize z
    model.setObjective(obj, GRB_MINIMIZE);

    // change some settings
    model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
    model.getEnv().set(GRB_IntParam_Threads, 1);
    model.getEnv().set(GRB_DoubleParam_TimeLimit, 3600);

    // find the optimal solution	
    model.optimize();

	// get reduced costs
	for (int j = 0; j < inst.n; j++) {    
		if(inst.items[j][2] == 0) continue; 		
			for (int i = 0; i <= inst.W - inst.items[j][0]; i++) { 
				if(inst.NPs[j][i]){  
					inst.RCs[j][i] = 0;
				}
            }
	}
		
	inst.contR = model.get(GRB_DoubleAttr_ObjVal);
}



