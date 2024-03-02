#include "LB.h"

void LB(const Instance& inst, Solution& sol) {
	
	// merge items with the same height
	vector<vector<int> > items2;
	vector<int> countItem (inst.W + 1, 0);
    for (int j = 0; j < inst.n; j++)
		countItem[inst.items[j][0]] += inst.items[j][1] * inst.items[j][2];
	for (int j = inst.W; j >= 0; j--){
		if(countItem[j] > 0)
			items2.push_back({j,countItem[j]});
	}

	// create item arcs
	vector<bool> isA (inst.W + 1, false); isA[0] = true;
	vector<vector<int> > arcs;
	for (int j = 0; j < items2.size(); j++) { 
		vector<bool> isS (inst.W  + 1, false);
		for (int i = inst.W ; i >= 0; i--){
			if(isA[i]){
				for (int k = 1; k <=  items2[j][1]; k++){
					if(i+items2[j][0]*k <= inst.W && !isS[i+items2[j][0]*(k-1)]){
						arcs.push_back({i+items2[j][0]*(k-1),i+items2[j][0]*k,j});
						isA[i+items2[j][0]*k] = true;	
						isS[i+items2[j][0]*(k-1)] = true;						
					}
					else
						break;
				}
			}
		}
	}
			
	GRBEnv env = GRBEnv();

	// create a model 	
	GRBModel model = GRBModel(env);         
	GRBLinExpr obj = 0;

	// declaration of the variables for the model
	vector<GRBVar> x (arcs.size());		

	// initizalization of the variables for the model
	for (int k = 0; k < arcs.size(); k++)
		x[k] = model.addVar(0, items2[arcs[k][2]][1], 0, GRB_CONTINUOUS);
	
	model.update();

	// create linear expressions
	vector<GRBLinExpr> assigned(items2.size(), 0); 
	vector<GRBLinExpr> cIn(inst.W + 1, 0);  
	vector<GRBLinExpr> cOut(inst.W + 1, 0); 	
	for (int k = 0; k < arcs.size(); k++){
		assigned[arcs[k][2]] += x[k];
		cIn[arcs[k][1]] += x[k];
		cOut[arcs[k][0]] += x[k];
	}

	model.update();
	
	// create assignment constraints
	for (int j = 0; j < items2.size(); j++)										
		model.addConstr(assigned[j] == items2[j][1]); 					

	// create flow conservation constraints 
	for (int i = 1; i < inst.W; i++)
		model.addConstr(cIn[i] >= cOut[i]);

	// set the objective: minimize z
	model.setObjective(cOut[0], GRB_MINIMIZE);

	// change some settings
	model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
	model.getEnv().set(GRB_IntParam_Threads, 1);
	model.getEnv().set(GRB_IntParam_Method, 2);
	model.getEnv().set(GRB_IntParam_MIPFocus, 1);
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

