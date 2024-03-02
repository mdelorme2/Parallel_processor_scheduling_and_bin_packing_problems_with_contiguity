#include "LB.h"

void LB(const Instance& inst, Solution& sol) {
	
	vector<vector<int> > items2 = inst.items;
	
	// scale the instance
	for (int j = 0; j < inst.n; j++)
		items2[j][0] *= 2;
	
	// compute normal patterms
	vector<vector<bool> > NPs (inst.n,vector<bool>(inst.W,false));
	for (int j = 0; j < inst.n; j++) {
		NPs[j][0] = true;
		for (int jp = 0; jp<inst.n;jp++){
			int lim = items2[jp][2];
			if (jp == j) lim--; 
			for(int l = 0; l<lim;l++){
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
    x.resize(inst.W, vector<GRBVar>(inst.n)); 
	GRBVar z = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
	GRBLinExpr obj = z;	
	
	// initizalization of the variables for the model
    for (int j = 0; j < inst.n; j++) { 
		for (int i = 0; i < inst.W; i++) {       	     
			if(NPs[j][i] && i <= inst.W - abs(inst.W - (i + items2[j][0]))){ 
				x[i][j] = model.addVar(0, items2[j][2], 0, GRB_CONTINUOUS);
			}
        }
    }
    model.update();

	// create linear expressions
    vector<GRBLinExpr> assigned(inst.n, 0); 
    vector<GRBLinExpr> height(inst.W, 0);   
    for (int j = 0; j < inst.n; j++) { 	
        for (int i = 0; i < inst.W; i++) { 
			if(NPs[j][i] && i <= inst.W - abs(inst.W - (i + items2[j][0]))){ 
				assigned[j] += x[i][j];	
				for (int ip = i; ip < min(inst.W,i+items2[j][0]); ip++) { 	
					height[ip] 	+= items2[j][1] * x[i][j];
				}
				for (int ip = inst.W-1; ip >= 2 * inst.W - (i + items2[j][0]); ip--) { 	
					height[ip] 	+= items2[j][1] * x[i][j];
				}
			}
        }
    }
    model.update();

    // create assignment constraints
    for (int j = 0; j < inst.n; j++)      										
		model.addConstr(assigned[j] == items2[j][2]); 						

    // create height constraints 
    for (int i = 0; i < inst.W; i++){		
		model.addConstr(height[i] <= 2*z);
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
	if (model.get(GRB_IntAttr_SolCount) >= 1){ 		
		sol.UB += ceil(model.get(GRB_DoubleAttr_ObjVal) - EPSILON);
		if(sol.LB == sol.UB) sol.opt = 1;
	}
}

