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

	vector<vector<int> > items2 = inst.items;
	
	// scale the instance
	for (int j = 0; j < inst.n; j++)
		items2[j][1] *= 2;

	// compute normal patterms
	vector<vector<bool> > NPs (inst.n,vector<bool>(UB + 1,false));
	for (int j = 0; j < inst.n; j++) {
		NPs[j][0] = true;
		for (int jp = 0; jp<inst.n;jp++){
			int lim = items2[jp][2];
			if (jp == j) lim--; 
			for(int l = 0; l<lim;l++){
				for (int i = UB-items2[jp][1] - 1;i >= 0; i--){
					if (NPs[j][i]) NPs[j][i+items2[jp][1]] = true;
				}
			}
		}
	}

	GRBEnv env = GRBEnv();
	for(;;){
		cout << "Try W = " << LB << endl;
		// create a model
		GRBEnv env = GRBEnv();              	
		GRBModel model = GRBModel(env);  
	
		// declaration of the variables for the model
		vector<vector<GRBVar> > x;				
		x.resize(LB, vector<GRBVar>(inst.n)); 
	
		// initizalization of the variables for the model
		for (int j = 0; j < inst.n; j++) { 
			for (int i = 0; i < LB; i++) {       	     
				if(NPs[j][i] && i <= LB - abs(LB - (i + items2[j][1]))){ 
					x[i][j] = model.addVar(0, items2[j][2], 0, GRB_CONTINUOUS);
				}
			}
		}
		model.update();

		// create linear expressions
		vector<GRBLinExpr> assigned(inst.n, 0); 
		vector<GRBLinExpr> height(LB, 0);   
		for (int j = 0; j < inst.n; j++) { 	
			for (int i = 0; i < LB; i++) { 
				if(NPs[j][i] && i <= LB - abs(LB - (i + items2[j][1]))){ 
					assigned[j] += x[i][j];	
					for (int ip = i; ip < min(LB,i+items2[j][1]); ip++) { 	
						height[ip] 	+= items2[j][0] * x[i][j];
					}
					for (int ip = LB-1; ip >= 2 * LB - (i + items2[j][1]); ip--) { 	
						height[ip] 	+= items2[j][0] * x[i][j];
					}
				}
			}
		}
		model.update();

		// create assignment constraints
		for (int j = 0; j < inst.n; j++)      										
			model.addConstr(assigned[j] == items2[j][2]); 						

		// create height constraints 
		for (int i = 0; i < LB; i++){		
			model.addConstr(height[i] <= 2*inst.W);
		}

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

