#include "LB.h"

void LB(const Instance& inst, Solution& sol) {
	
	// create node matrix
	vector<vector<int> > matToId(inst.items.size()+1,vector<int>(inst.W+1,-1));
	vector<vector<int> > idToMat; 
	int nodeCount = 0;
	
	// create item arcs
	vector<bool> isA (inst.W+1, false); isA[0] = true;
	vector<vector<int> > arcs;
	int tail,head;
	for (int j = 0; j < inst.items.size(); j++) { 
		for (int i = inst.W ; i >= 0; i--){
			if(isA[i]){
				if(matToId[j][i] == -1){
					matToId[j][i] = nodeCount;
					idToMat.push_back({j,i});
					nodeCount++;					
				}
				tail = matToId[j][i];
				for (int k = 0; k <= inst.items[j][2]; k++){
					if(i+inst.items[j][0]*k <= inst.W){
						if(matToId[j+1][i+inst.items[j][0]*k] == -1){
							matToId[j+1][i+inst.items[j][0]*k] = nodeCount;
							idToMat.push_back({j+1,i+inst.items[j][0]*k});
							nodeCount++;					
						}
						head = matToId[j+1][i+inst.items[j][0]*k];
						arcs.push_back({tail,head,j,k});
						isA[i+inst.items[j][0]*k] = true;							
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
		x[k] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
	
	model.update();

	// create linear expressions
	vector<GRBLinExpr> assigned(inst.items.size(), 0); 
	vector<GRBLinExpr> cIn(idToMat.size(), 0);  
	vector<GRBLinExpr> cOut(idToMat.size(), 0); if(cOut.size() == 0) cOut.resize(1,0);	  		
	for (int k = 0; k < arcs.size(); k++){
		assigned[arcs[k][2]] += x[k] * arcs[k][3];
		cIn[arcs[k][1]] += x[k];
		cOut[arcs[k][0]] += x[k];
	}
    
	model.update();

	// create assignment constraints
	for (int j = 0; j < inst.items.size(); j++)										
		model.addConstr(assigned[j] == inst.items[j][2] * inst.items[j][1]); 					

	// create flow conservation constraints 
	for (int i = 1; i < idToMat.size(); i++)
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

