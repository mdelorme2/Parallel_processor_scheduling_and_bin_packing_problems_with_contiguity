#include "model.h"

void CPM(const Instance& inst, Solution& sol) {
	
	// from CSP to BPP
	vector<vector<int> > items2;
	int n2 = 0;
    for (int j = 0; j < inst.n; j++) { 
		for (int k = 0; k < inst.items[j][2]; k++){
			items2.push_back({inst.items[j][0],inst.items[j][1],j});
			n2++;
		}			
	}

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
	IloEnv env;     
	IloModel model(env);

	// declaration of the variables for the model
	IloIntervalVarArray x(env, n2);
	IloCumulFunctionExpr z (env);
    // initizalization of the variables for the model
    for (int j = 0; j < n2; j++) {	
		IloNumToNumStepFunction IF(env, 0, inst.W - items2[j][0] + 1, 100);
		for (int i = 0; i <= inst.W - items2[j][0]; i++){       	     
			if(NPs[j][i] == false) IF.setValue(i, i+1, 0);
		}
		x[j] = IloIntervalVar(env); 
		x[j].setStartMin(0);
		x[j].setStartMax(inst.W - items2[j][0]);
		//if(j > 0 && items2[j][2] == items2[j-1][2])  model.add(IloStartBeforeStart(env, x[j-1], x[j]));
		x[j].setSizeMin(items2[j][0]);
		x[j].setSizeMax(items2[j][0]);
		z += IloPulse(x[j], items2[j][1]);
		model.add(IloForbidStart(env, x[j], IF));
    }

    // set the objective: minimize z
	IloIntVar obj (env, 0, IloInfinity);
	model.add(z <= obj);
	IloObjective objective = IloMinimize(env, obj);
    model.add(objective);

    // change some settings
    IloCP cp(model);
    cp.setParameter(IloCP::TimeLimit, 3600);
    cp.setParameter(IloCP::Workers, 1);
	cp.setParameter(IloCP::LogPeriod, 500000); 

	// store the results in a Solution object
	sol.opt = 0;
	
    // find the optimal solution    
	if(cp.solve()){
		sol.UB += cp.getObjValue();
		sol.LB += cp.getObjBound();
		if(cp.getStatus() == 2){	
			sol.opt = 1;
			sol.LB = sol.UB;
		}
		
        // get bin for each item
        for (int j = 0; j < n2; j++){                                 
			cout << "Task "<< j << " starts at time " << cp.getStart(x[j]) << endl;
			sol.assignedBin[items2[j][2]].push_back(cp.getStart(x[j]));
		}
	}
}

