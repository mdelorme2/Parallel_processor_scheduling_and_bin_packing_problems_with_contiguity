#include "model.h"

void CPM(const Instance& inst, Solution& sol, const int& seed) {
	
	// from CSP to BPP
	vector<vector<int> > items2;
	int n2 = 0;
    for (int j = 0; j < inst.n; j++) { 
		for (int k = 0; k < inst.items[j][2]; k++){
			items2.push_back({inst.items[j][0],inst.items[j][1],j});
			n2++;
		}			
	}

	int LUB = CPMUB(inst,seed);
	
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
		if(j == 0)
			x[j].setStartMax((LUB - items2[j][1])/2);
		else
			x[j].setStartMax(LUB - items2[j][1]);
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
	model.add(z >= LB(inst));
	IloObjective objective = IloMinimize(env, z);
    model.add(objective);

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

int CPMUB(const Instance& inst, const int& seed) {
	
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
		if(j > 0 && items2[j][2] == items2[j-1][2])  model.add(IloStartBeforeStart(env, x[j-1], x[j]));
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
    cp.setParameter(IloCP::SolutionLimit, 1);
    cp.setParameter(IloCP::Workers, 1);
	cp.setParameter(IloCP::LogPeriod, 500000); 
	cp.setParameter(IloCP::RandomSeed, seed);
	cp.solve();
	
    return cp.getObjValue();
}

int LB(const Instance& inst) {
	
	// merge items with the same height
	vector<vector<int> > items2;
	vector<int> countItem (inst.W + 1, 0);
	vector<int> countItem2 (inst.W + 1, 0);
    for (int j = 0; j < inst.n; j++){
		countItem[inst.items[j][0]] += inst.items[j][1] * inst.items[j][2];
		countItem2[inst.items[j][0]] += inst.items[j][2];
	}
	for (int j = inst.W; j >= 0; j--){
		if(countItem[j] > 0)
			items2.push_back({j,countItem[j],countItem2[j]});
	}

	// create node matrix
	vector<vector<int> > matToId(items2.size()+1,vector<int>(inst.W+1,-1));
	vector<vector<int> > idToMat; 
	int nodeCount = 0;
	
	// create item arcs
	vector<bool> isA (inst.W+1, false); isA[0] = true;
	vector<vector<int> > arcs;
	int tail,head;
	for (int j = 0; j < items2.size(); j++) { 
		for (int i = inst.W ; i >= 0; i--){
			if(isA[i]){
				if(matToId[j][i] == -1){
					matToId[j][i] = nodeCount;
					idToMat.push_back({j,i});
					nodeCount++;					
				}
				tail = matToId[j][i];
				for (int k = 0; k <=  items2[j][2]; k++){
					if(i+items2[j][0]*k <= inst.W){
						if(matToId[j+1][i+items2[j][0]*k] == -1){
							matToId[j+1][i+items2[j][0]*k] = nodeCount;
							idToMat.push_back({j+1,i+items2[j][0]*k});
							nodeCount++;					
						}
						head = matToId[j+1][i+items2[j][0]*k];
						arcs.push_back({tail,head,j,k});
						isA[i+items2[j][0]*k] = true;	
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
	vector<GRBLinExpr> assigned(items2.size(), 0); 
	vector<GRBLinExpr> cIn(idToMat.size(), 0);  
	vector<GRBLinExpr> cOut(idToMat.size(), 0); if(cOut.size() == 0) cOut.resize(1,0);	
	for (int k = 0; k < arcs.size(); k++){
		assigned[arcs[k][2]] += x[k] * arcs[k][3];
		cIn[arcs[k][1]] += x[k];
		cOut[arcs[k][0]] += x[k];
	}

	model.update();

	// create assignment constraints
	for (int j = 0; j < items2.size(); j++){										
		model.addConstr(assigned[j] == items2[j][1]); 					
	}
	// create flow conservation constraints 
	for (int i = 1; i < idToMat.size(); i++)
		model.addConstr(cIn[i] >= cOut[i]);

	// set the objective: minimize z
	model.setObjective(cOut[0], GRB_MINIMIZE);

	// change some settings
	model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
	model.getEnv().set(GRB_IntParam_Threads, 1);
	model.getEnv().set(GRB_IntParam_Method, 2);
	model.getEnv().set(GRB_DoubleParam_TimeLimit, 3600);

	// find the optimal solution		
	model.optimize();

	return ceil(model.get(GRB_DoubleAttr_ObjVal) - EPSILON);
}
