#include "LB.h"

void LB(const Instance& inst, Solution& sol) {
	
	// compute bounds
	int maxH = 0;
	int sumArea = 0;
	int LB = 0;
	int UB = 0;
	for (int j = 0; j < inst.n; j++){
		sumArea += inst.items[j][0] * inst.items[j][1] * inst.items[j][2];	
		maxH = max(maxH,inst.items[j][1]);
		UB += inst.items[j][1] * inst.items[j][2];
	}
	LB = sumArea/inst.W;

	// merge items with the same height
	vector<vector<int> > items2;
	vector<int> countItem (maxH + 1, 0);
    for (int j = 0; j < inst.n; j++)
		countItem[inst.items[j][1]] += inst.items[j][0] * inst.items[j][2];
	for (int j = maxH; j >= 0; j--){
		if(countItem[j] > 0)
			items2.push_back({j,countItem[j]});
	}

	// create item arcs
	vector<bool> isA (UB  + 1, false); isA[0] = true;
	vector<vector<int> > arcs;
	for (int j = 0; j < items2.size(); j++) { 
		vector<bool> isS (UB  + 1, false);
		for (int i = UB ; i >= 0; i--){
			if(isA[i]){
				for (int k = 1; k <=  items2[j][1]; k++){
					if(i+items2[j][0]*k <= UB && !isS[i+items2[j][0]*(k-1)]){
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
	for(;;){
		cout << "Try H = " << LB << endl;
		
		// create a model           
		GRBModel model = GRBModel(env);         
		GRBLinExpr obj = 0;
	
		// declaration of the variables for the model
		vector<GRBVar> x (arcs.size());		
	
		// initizalization of the variables for the model
		for (int k = 0; k < arcs.size(); k++)
			if(arcs[k][1] <= LB) x[k] = model.addVar(0, items2[arcs[k][2]][1], 0, GRB_CONTINUOUS);
		
		model.update();

		// create linear expressions
		vector<GRBLinExpr> assigned(items2.size(), 0); 
		vector<GRBLinExpr> cIn(LB+1, 0);  
		vector<GRBLinExpr> cOut(LB+1, 0); 	
		for (int k = 0; k < arcs.size(); k++){
			if(arcs[k][1] <= LB){
				assigned[arcs[k][2]] += x[k];
				cIn[arcs[k][1]] += x[k];
				cOut[arcs[k][0]] += x[k];
			}
		}
    
		model.update();

		// create assignment constraints
		for (int j = 0; j < items2.size(); j++)										
			model.addConstr(assigned[j] == items2[j][1]); 					

		// create flow conservation constraints 
		for (int i = 1; i < LB; i++)
			model.addConstr(cIn[i] >= cOut[i]);

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

