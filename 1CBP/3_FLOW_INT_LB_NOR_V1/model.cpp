#include "model.h"

void FLOW(const Instance& inst, Solution& sol) {
	
	// local variable
	int LUB = CPM(inst);
	int LLB = LB(inst);
	vector<bool> isHActive(LUB,false); 
	vector<int> isHActiveV; 

	// compute normal patterms
	vector<vector<bool> > NPs (inst.n,vector<bool>(LUB,false));
	for (int j = 0; j < inst.n; j++) {
		if(inst.items[j][2] == 0) continue; 
		NPs[j][0] = true;
		for (int jp = 0; jp<inst.n;jp++){
			int lim = inst.items[jp][2];
			if (jp == j) lim--; 
			for(int l = 0; l<lim;l++){
				for (int i = LUB-inst.items[jp][1] - 1;i >= 0; i--){
					if (NPs[j][i]) NPs[j][i+inst.items[jp][1]] = true;
				}
			}
		}
	}
	
	// create item arcs
	vector<vector<int> > arcs;
	for (int j = 0; j < inst.n; j++) {
		if(inst.items[j][2] == 0) continue; 
		for (int i = 0; i <= LUB - inst.items[j][1]; i++){
			if(NPs[j][i]){
				arcs.push_back({i,i+inst.items[j][1],inst.items[j][0],j});
				isHActive[i] = true;
				isHActive[i+inst.items[j][1]] = true;
			}
		}
	}
	cout << "There are " << arcs.size() << " arcs, loss arcs excluded" << endl;

	// create loss arcs
	for(int i = 0; i<LUB;i++)
		if(isHActive[i]) isHActiveV.push_back(i);
	isHActiveV.push_back(LUB);	
	for (int i = 0; i < isHActiveV.size()-1; i++)
		arcs.push_back({isHActiveV[i],isHActiveV[i+1],1,-1,1});
	cout << "There are " << arcs.size() << " arcs, loss arcs included" << endl;
	
    // create a model
    GRBEnv env = GRBEnv();              	
    GRBModel model = GRBModel(env);         

    // declaration of the variables for the model
	vector<GRBVar> x (arcs.size());		
    vector<GRBVar> zvec(LUB+1);
	GRBLinExpr obj = 0;
	
    // initizalization of the variables for the model
    for (int k = 0; k < arcs.size(); k++){
        if (arcs[k][3] >= 0) x[k] = model.addVar(0, inst.items[arcs[k][3]][2], 0, GRB_INTEGER);
		else x[k] = model.addVar(0, GRB_INFINITY, 0, GRB_INTEGER);
 	}
	for (int i = 0; i < LLB; i++) 
		zvec[i] =  model.addVar(0, 0, 0, GRB_BINARY);
	for (int i = LLB; i <= LUB; i++){ 
		if(i == LUB || isHActive[i])
			zvec[i] =  model.addVar(0, 1, 0, GRB_BINARY);
		else 
			zvec[i] =  model.addVar(0, 0, 0, GRB_BINARY);
	}
    model.update();

    // create linear expressions
    vector<GRBLinExpr> assigned(inst.n, 0); 
    vector<GRBLinExpr> cIn(LUB+1, 0);  
	vector<GRBLinExpr> cOut(LUB+1, 0); 	
	vector<vector<GRBLinExpr> > startsHere(LUB+1, vector<GRBLinExpr>(inst.n,0)); 
	vector<GRBLinExpr> endsHere(LUB+1, 0); 
	GRBLinExpr startsHere0 = 0;
    for (int k = 0; k < arcs.size(); k++){
		cIn[arcs[k][1]] += arcs[k][2] * x[k];
        cOut[arcs[k][0]] += arcs[k][2] * x[k];
		if(arcs[k][3] >= 0){
			assigned[arcs[k][3]] += x[k];
			startsHere[arcs[k][0]][arcs[k][3]] += x[k];
			endsHere[arcs[k][1]] += x[k];
		}
    }
    
    model.update();

    // create assignment constraints
	int nbTot = 0; 
    for (int j = 0; j < inst.n; j++){
		if(inst.items[j][2] == 0) continue; 		
		model.addConstr(assigned[j] == inst.items[j][2]);
		nbTot += inst.items[j][2];		
		startsHere0 += startsHere[0][j];
	}
	
    // create flow conservation constraints 
    for (int i = 1; i <= LUB; i++){
			obj += i * zvec[i];
			model.addConstr(cIn[i] == cOut[i] + inst.W * zvec[i]);
	}
	model.addConstr(cOut[0] == inst.W - zvec[0]);

	// symmetry breaking constraints
	if(nbTot > 0) model.addConstr(startsHere0 >= 1);
    for (int j = 0; j < inst.n; j++) { 	
		if(inst.items[j][2] == 0) continue; 
        for (int i = 1; i <= LUB - inst.items[j][1]; i++) { 
			if(NPs[j][i]){
				model.addConstr(startsHere[i][j] <= inst.items[j][2] * endsHere[i]);
			}
        }
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
    if (model.get(GRB_IntAttr_SolCount) >= 1) { 		
        sol.UB += ceil(model.get(GRB_DoubleAttr_ObjVal) - EPSILON);
		if(sol.LB == sol.UB) sol.opt = 1;
		
        // get bin for each item
        for (int k = 0; k < arcs.size(); k++){    
			if(arcs[k][3] >= 0){
				for (int l = 0; l < ceil(x[k].get(GRB_DoubleAttr_X) - EPSILON); l++)  
					sol.assignedBin[arcs[k][3]].push_back(arcs[k][0]);
            }
        }
    }
	/*GRBModel model2 = model.relax();
	model2.optimize();
	sol.contR += model2.get(GRB_DoubleAttr_ObjVal);*/
}

int CPM(const Instance& inst) {
	
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