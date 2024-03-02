#include "model.h"

void BM(const Instance& inst, Solution& sol) {
	
	// local variable
	vector<bool> isWActive(inst.W,false);

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
	
    // create a model
    GRBEnv env = GRBEnv();              	
    GRBModel model = GRBModel(env);         

    // declaration of the variables for the model
    vector<vector<vector<GRBVar> > > x(inst.W); for(int i = 0; i < inst.W;i++) x[i].resize(inst.n);
	GRBVar z = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
	GRBLinExpr obj = z;
	
    // initizalization of the variables for the model
    for (int j = 0; j < inst.n; j++) { 
		if(inst.items[j][2] == 0) continue; 
		for (int i = 0; i <= inst.W - inst.items[j][0]; i++) {   
			if(NPs[j][i]){			
				isWActive[i] = true;
				x[i][j].resize(ceil(log2(inst.items[j][2]))+1);
				for (int k = 0; k < ceil(log2(inst.items[j][2]))+1; k++)
					x[i][j][k] = model.addVar(0, 1, 0, GRB_BINARY);
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
			if(NPs[j][i]){
				for (int k = 0; k < ceil(log2(inst.items[j][2]))+1; k++){
					assigned[j] += pow(2,k) * x[i][j][k];	
					for (int ip = i; ip < i+inst.items[j][0]; ip++) { 
						height[ip] 	+= inst.items[j][1] * x[i][j][k] * pow(2,k);
					}
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
        for (int j = 0; j < inst.n; j++) { 
			if(inst.items[j][2] == 0) continue; 		
			for (int i = 0; i <= inst.W - inst.items[j][0]; i++) { 
				if(NPs[j][i]){
					for (int k = 0; k < ceil(log2(inst.items[j][2]))+1; k++){
						if(ceil(x[i][j][k].get(GRB_DoubleAttr_X) - EPSILON) == 1){ 
							for(int l = 0; l < pow(2,k);l++){
								sol.assignedBin[j].push_back(i);
							}
						}
					}
				}
            }
        }
    }
	/*GRBModel model2 = model.relax();
	model2.optimize();
	sol.contR += model2.get(GRB_DoubleAttr_ObjVal);*/
}

