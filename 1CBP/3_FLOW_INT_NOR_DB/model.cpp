#include "model.h"

int FLOW(const Instance& inst, Solution& sol, const int& DB, const int& TL) {
	cout << "Test " << DB << "(" << inst.contR << ") with " << TL << "s left" << endl;	
	
	vector<bool> isHActive(DB,false); 
	
    // create a model
    GRBEnv env = GRBEnv();              	
    GRBModel model = GRBModel(env);         

    // declaration of the variables for the model
	vector<GRBVar> x (inst.arcs.size());		
	
    // initizalization of the variables for the model
    for (int k = 0; k < inst.arcs.size(); k++){
		if(inst.arcs[k][1] > DB || inst.RCs[k] + inst.contR - EPSILON > DB) continue;
		if (inst.arcs[k][3] >= 0) x[k] = model.addVar(0, inst.items[inst.arcs[k][3]][2], 0, GRB_INTEGER);
		else x[k] = model.addVar(0, GRB_INFINITY, 0, GRB_INTEGER);
 	}
    model.update();

    // create linear expressions
    vector<GRBLinExpr> assigned(inst.n, 0); 
    vector<GRBLinExpr> cIn(DB+1, 0);  
	vector<GRBLinExpr> cOut(DB+1, 0); 	
    for (int k = 0; k < inst.arcs.size(); k++){
		if(inst.arcs[k][1] > DB || inst.RCs[k] + inst.contR - EPSILON > DB) continue;
		cIn[inst.arcs[k][1]] += inst.arcs[k][2] * x[k];
        cOut[inst.arcs[k][0]] += inst.arcs[k][2] * x[k];
		if(inst.arcs[k][3] >= 0) assigned[inst.arcs[k][3]] += x[k];
    }
    
    model.update();

    // create assignment constraints
    for (int j = 0; j < inst.n; j++){
		if(inst.items[j][2] == 0) continue; 		
		model.addConstr(assigned[j] == inst.items[j][2]);
	}
	
    // create flow conservation constraints 
    for (int i = 1; i < DB; i++){
		model.addConstr(cIn[i] == cOut[i]);
	}
	model.addConstr(cOut[0] == cIn[DB]);
	model.addConstr(cOut[0] <= inst.W);
	
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
    if (model.get(GRB_IntAttr_SolCount) >= 1) { 		
		sol.LB += DB;
        sol.UB += DB;
		sol.opt = 1;
		
       // get bin for each item
        for (int k = 0; k < inst.arcs.size(); k++){  
			if(inst.arcs[k][1] > DB || inst.RCs[k] + inst.contR - EPSILON > DB) continue;
			if(inst.arcs[k][3] >= 0){
				for (int l = 0; l < ceil(x[k].get(GRB_DoubleAttr_X) - EPSILON); l++)  
					sol.assignedBin[inst.arcs[k][3]].push_back(inst.arcs[k][0]);
            }
        }
		return 0;
    }
	else{
		if(model.get(GRB_IntAttr_Status) == 9){
			sol.LB += DB;
			return -2;
		}
		return -1;
	}
}

void FLOWLP(Instance& inst) {
	
    // create a model
	vector<bool> isHActive(inst.LUB,false);

    GRBEnv env = GRBEnv();              	
    GRBModel model = GRBModel(env);         

    // declaration of the variables for the model
	vector<GRBVar> x (inst.arcs.size());		
    vector<GRBVar> zvec(inst.LUB+1);
	GRBLinExpr obj = 0;
	
    // initizalization of the variables for the model
    for (int k = 0; k < inst.arcs.size(); k++){
        if (inst.arcs[k][3] >= 0) x[k] = model.addVar(0, inst.items[inst.arcs[k][3]][2], 0, GRB_CONTINUOUS);
		else x[k] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
 	}
	for (int i = 0; i <= inst.LUB; i++){ 
		zvec[i] =  model.addVar(0, 1, 0, GRB_CONTINUOUS);
	}
    model.update();

    // create linear expressions
    vector<GRBLinExpr> assigned(inst.n, 0); 
    vector<GRBLinExpr> cIn(inst.LUB+1, 0);  
	vector<GRBLinExpr> cOut(inst.LUB+1, 0); 	
    for (int k = 0; k < inst.arcs.size(); k++){
		cIn[inst.arcs[k][1]] += inst.arcs[k][2] * x[k];
        cOut[inst.arcs[k][0]] += inst.arcs[k][2] * x[k];
		if(inst.arcs[k][3] >= 0) assigned[inst.arcs[k][3]] += x[k];
    }
    
    model.update();

    // create assignment constraints
    for (int j = 0; j < inst.n; j++){
		if(inst.items[j][2] == 0) continue; 		
		model.addConstr(assigned[j] == inst.items[j][2]);
	}
	
    // create flow conservation constraints 
    for (int i = 1; i <= inst.LUB; i++){
			obj += i * zvec[i];
			model.addConstr(cIn[i] == cOut[i] + inst.W * zvec[i]);
	}
	model.addConstr(cOut[0] == inst.W - zvec[0]);
	
	// set the objective: minimize z
    model.setObjective(obj, GRB_MINIMIZE);

    // change some settings
    model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
    model.getEnv().set(GRB_IntParam_Threads, 1);
    model.getEnv().set(GRB_DoubleParam_TimeLimit, 3600);

    // find the optimal solution	
    model.optimize();

	// get reduced costs
	for (int k = 0; k < inst.arcs.size(); k++){    
		inst.RCs[k] = 0;  
	}
		
	inst.contR = model.get(GRB_DoubleAttr_ObjVal);
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
