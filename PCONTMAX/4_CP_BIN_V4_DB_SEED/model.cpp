#include "model.h"

int CPM(const Instance& inst, Solution& sol, const int& DB, const int& TL, const int& seed) {
	
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
	IloCumulFunctionExpr z (env);
	
    // initizalization of the variables for the model
    for (int j = 0; j < n2; j++) { 
		x[j] = IloIntervalVar(env); 
		x[j].setStartMin(0);
		x[j].setStartMax(inst.W - items2[j][0]);
		if(j > 0 && items2[j][2] == items2[j-1][2])  model.add(IloStartBeforeStart(env, x[j-1], x[j]));
	    x[j].setSizeMin(items2[j][0]);
		x[j].setSizeMax(items2[j][0]);
		z += IloPulse(x[j], items2[j][1]);
    }

	
    // set the objective: minimize z
	model.add(z <= DB);
	
    // change some settings
    IloCP cp(model);
    cp.setParameter(IloCP::TimeLimit, 3600);
    cp.setParameter(IloCP::Workers, 1);
	cp.setParameter(IloCP::LogPeriod, 500000); 
	cp.setParameter(IloCP::RandomSeed, seed);
	
	// store the results in a Solution object
	sol.opt = 0;
	
    // find the optimal solution    
	if(cp.solve()){
		sol.UB += DB;
		sol.LB += DB;
		sol.opt = 1;
		
        // get bin for each item
        for (int j = 0; j < n2; j++){                                 
			cout << "Task "<< j << " starts at time " << cp.getStart(x[j]) << endl;
			sol.assignedBin[items2[j][2]].push_back(cp.getStart(x[j]));
		}
		return 0;
	}
	else{
		if (cp.getStatus() != 3)
			sol.LB += DB;
		return -1;
	}		
}

void BMLP(Instance& inst) {
	
	// local variable
	vector<bool> isWActive(inst.W,false);

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
		if(inst.items[j][2] == 0) continue; 
		for (int i = 0; i <= inst.W - inst.items[j][0]; i++) {       	     
			isWActive[i] = true;
			x[i][j] = model.addVar(0, inst.items[j][2], 0, GRB_CONTINUOUS);
        }
    }
    model.update();

    // create linear expressions
    vector<GRBLinExpr> assigned(inst.n, 0);  
    vector<GRBLinExpr> height(inst.W, 0);   
    for (int j = 0; j < inst.n; j++) { 	
		if(inst.items[j][2] == 0) continue; 
        for (int i = 0; i <= inst.W - inst.items[j][0]; i++) { 
			assigned[j] += x[i][j];	
			for (int ip = i; ip < i+inst.items[j][0]; ip++) { 	
				height[ip] 	+= inst.items[j][1] * x[i][j];
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
	inst.contR = model.get(GRB_DoubleAttr_ObjVal);
}



