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

	// create a model
	IloEnv env;     
	IloModel model(env);

	// declaration of the variables for the model
	IloIntervalVarArray x(env, n2);
	IloCumulFunctionExpr cumul (env);
	IloIntExprArray ends(env);
	
    // initizalization of the variables for the model
    for (int j = 0; j < n2; j++) { 
		x[j] = IloIntervalVar(env); 
		x[j].setStartMin(0);
	//	if(j > 0 && items2[j][2] == items2[j-1][2])  model.add(IloStartBeforeStart(env, x[j-1], x[j]));
	    x[j].setSizeMin(items2[j][1]);
		x[j].setSizeMax(items2[j][1]);
		cumul += IloPulse(x[j], items2[j][0]);
		ends.add(IloEndOf(x[j]));
    }

	// capacity constraints
	model.add(cumul <= inst.W);	
	
    // set the objective: minimize the makespan
	IloIntVar z (env, 0, IloInfinity); 
	model.add(z >= IloMax(ends));
	IloObjective objective = IloMinimize(env, z);
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

