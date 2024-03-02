#include "LB.h"

void LB(const Instance& inst, Solution& sol) {
	
	// compute bounds
	int sumArea = 0;
	int LB = 0;
	int UB = 0;
	for (int j = 0; j < inst.n; j++){
		sumArea += inst.items[j][0] * inst.items[j][1] * inst.items[j][2];	
		UB += inst.items[j][1] * inst.items[j][2];
	}
	LB = sumArea/inst.W;

	// create node matrix
	vector<vector<int> > matToId(inst.items.size()+1,vector<int>(UB+1,-1));
	vector<vector<int> > idToMat; 
	int nodeCount = 0;
	
	// create item arcs
	vector<bool> isA (UB  + 1, false); isA[0] = true;
	vector<vector<int> > arcs;
	int tail,head;
	for (int j = 0; j < inst.items.size(); j++) { 
		for (int i = UB ; i >= 0; i--){
			if(isA[i]){
				if(matToId[j][i] == -1){
					matToId[j][i] = nodeCount;
					idToMat.push_back({j,i});
					nodeCount++;					
				}
				tail = matToId[j][i];
				for (int k = 0; k <= inst.items[j][2]; k++){
					if(i+inst.items[j][1]*k <= UB){
						if(matToId[j+1][i+inst.items[j][1]*k] == -1){
							matToId[j+1][i+inst.items[j][1]*k] = nodeCount;
							idToMat.push_back({j+1,i+inst.items[j][1]*k});
							nodeCount++;					
						}
						head = matToId[j+1][i+inst.items[j][1]*k];
						arcs.push_back({tail,head,j,k});
						isA[i+inst.items[j][1]*k] = true;							
					}
					else
						break;
				}
			}
		}
	}
			
	GRBEnv env = GRBEnv();
	for(;;){
		cout << "Try H = " << LB << endl;
		
		// create a model           
		GRBModel model = GRBModel(env);         
		GRBLinExpr obj = 0;
	
		// declaration of the variables for the model
		vector<GRBVar> x (arcs.size());		
	
		// initizalization of the variables for the model
		for (int k = 0; k < arcs.size(); k++){
			if(idToMat[arcs[k][1]][1] <= LB) 
				x[k] = model.addVar(0, inst.W, 0, GRB_INTEGER);
		}
		model.update();

		// create linear expressions
		vector<GRBLinExpr> assigned(inst.items.size(), 0); 
		vector<GRBLinExpr> cIn(idToMat.size(), 0);  
		vector<GRBLinExpr> cOut(idToMat.size(), 0); 		
		for (int k = 0; k < arcs.size(); k++){
			if(idToMat[arcs[k][1]][1] <= LB){
				assigned[arcs[k][2]] += x[k] * arcs[k][3];
				cIn[arcs[k][1]] += x[k];
				cOut[arcs[k][0]] += x[k];
			}
		}
    
		model.update();

		// create assignment constraints
		for (int j = 0; j < inst.items.size(); j++)										
			model.addConstr(assigned[j] == inst.items[j][2] * inst.items[j][0] ); 					

		// create flow conservation constraints 
		for (int i = 1; i < idToMat.size(); i++){
			if(idToMat[i][1] <= LB)
				model.addConstr(cIn[i] >= cOut[i]);
		}

		// set the objective: minimize z
		model.setObjective(cOut[0], GRB_MINIMIZE);

		// change some settings
		model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
		model.getEnv().set(GRB_IntParam_Threads, 1);
		model.getEnv().set(GRB_IntParam_Method, 2);
		model.getEnv().set(GRB_IntParam_MIPFocus, 1);
		model.getEnv().set(GRB_DoubleParam_TimeLimit, 3600);
		model.getEnv().set(GRB_DoubleParam_Cutoff, inst.W + 0.0001);

		// find the optimal solution		
		model.optimize();

		// store the results in a Solution object
		sol.Nvar = model.get(GRB_IntAttr_NumVars);       
		sol.Nconstr = model.get(GRB_IntAttr_NumConstrs); 
		sol.Ncoeff = model.get(GRB_IntAttr_NumNZs);      	
		sol.opt = 0;
		
		// if a solution has been found
		if (model.get(GRB_IntAttr_SolCount) >= 1){ 		
			/*for (int k = 0; k < arcs.size(); k++){
				if(idToMat[arcs[k][1]][1] <= LB){
					for (int j = 0; j < ceil(x[k].get(GRB_DoubleAttr_X) - EPSILON); j++){
						cout << idToMat[arcs[k][0]][0] << " " << idToMat[arcs[k][0]][1] << " " << idToMat[arcs[k][1]][0] << " " << idToMat[arcs[k][1]][1] << " " << arcs[k][2] << " " << arcs[k][3] << endl;
					}
				}
			}*/		
			UB = LB;
			sol.opt = 1;
			break;
		}
		else{
			if (model.get(GRB_IntAttr_Status) == 9)
				break;
			else 
				LB++;
		}
	}
	
	sol.UB += UB;
	sol.LB += LB;
}

