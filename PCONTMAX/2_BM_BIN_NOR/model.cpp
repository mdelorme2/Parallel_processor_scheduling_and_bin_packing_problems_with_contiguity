#include "model.h"

void BM(const Instance& inst, Solution& sol) {
	
	// from CSP to BPP
	vector<vector<int> > items2;
	int n2 = 0;
    for (int j = 0; j < inst.n; j++) { 
		for (int k = 0; k < inst.items[j][2]; k++){
			items2.push_back({inst.items[j][0],inst.items[j][1],j});
			n2++;
		}			
	}

	// local variable
	vector<bool> isWActive(inst.W,false);

	// compute normal patterms
	vector<vector<bool> > NPs (n2,vector<bool>(inst.W,false));
	for (int j = 0; j < n2; j++) {
		NPs[j][0] = true;
		for (int jp = 0; jp<n2;jp++){
			if(jp != j){
				for (int i = inst.W-items2[jp][0] - 1;i >= 0; i--){
					if (NPs[j][i]) NPs[j][i+items2[jp][0]] = true;
				}
			}
		}
	}
	
    // create a model
    GRBEnv env = GRBEnv();              	
    GRBModel model = GRBModel(env);         

    // declaration of the variables for the model
    vector<vector<GRBVar> > x;			
    x.resize(inst.W, vector<GRBVar>(n2)); 
	GRBVar z = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
	GRBLinExpr obj = z;
	
    // initizalization of the variables for the model
    for (int j = 0; j < n2; j++) { 
		for (int i = 0; i <= inst.W - items2[j][0]; i++) {
			if(NPs[j][i]){ 
				isWActive[i] = true;
				x[i][j] = model.addVar(0, 1, 0, GRB_BINARY);
			}
		}
    }
    model.update();

    // create linear expressions
    vector<GRBLinExpr> assigned(n2, 0);  
    vector<GRBLinExpr> height(inst.W, 0);   
    for (int j = 0; j < n2; j++) { 	
        for (int i = 0; i <= inst.W - items2[j][0]; i++){
			if(NPs[j][i]){ 
				assigned[j] += x[i][j];	
				for (int ip = i; ip < i+items2[j][0]; ip++){ 	
					height[ip] 	+= items2[j][1] * x[i][j];
				}
			}
		}
    }   
    model.update();

    // create assignment constraints
    for (int j = 0; j < n2; j++){      										
		model.addConstr(assigned[j] == 1); 						
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
        for (int j = 0; j < n2; j++) {                                        
			for (int i = 0; i <= inst.W- items2[j][0]; i++) {
				if(NPs[j][i] && ceil(x[i][j].get(GRB_DoubleAttr_X) - EPSILON) == 1) 
					sol.assignedBin[items2[j][2]].push_back(i);
            }
        }
    }
	/*GRBModel model2 = model.relax();
	model2.optimize();
	sol.contR += model2.get(GRB_DoubleAttr_ObjVal);*/
}

