#include "model.h"

int APTP(const Instance& inst, Solution& sol, const int& DB, const int& TL) {
	cout << "Test " << DB << "(" << inst.contR << ") with " << TL << "s left" << endl;
	
	vector<bool> isHActive(DB,false);
	
    GRBEnv env = GRBEnv();              	
    GRBModel model = GRBModel(env);         

    // declaration of the variables for the model
    vector<vector<GRBVar> > x(DB, vector<GRBVar>(inst.n));			
	
    // initizalization of the variables for the model
    for (int j = 0; j < inst.n; j++) {  
		if(inst.items[j][2] == 0) continue; 
		for (int i = 0; i <= DB - inst.items[j][1]; i++) {
			if(inst.NPs[j][i] && inst.RCs[j][i] + inst.contR - EPSILON <= DB){
				isHActive[i] = true;
				x[i][j] = model.addVar(0, inst.items[j][2], 0, GRB_INTEGER);
			}
        }
    }

    // create linear expressions
	vector<GRBLinExpr> load(DB, 0); 
	vector<GRBLinExpr> assigned(inst.n, 0);
	
    for (int j = 0; j < inst.n; j++) { 	
		if(inst.items[j][2] == 0) continue; 
        for (int i = 0; i <= DB - inst.items[j][1]; i++){
			if(inst.NPs[j][i] && inst.RCs[j][i] + inst.contR - EPSILON <= DB){
				for (int ip = i; ip < i+inst.items[j][1]; ip++){ 	
					load[ip] += inst.items[j][0] * x[i][j];
				}
				assigned[j] += x[i][j];
			}
		}
    }   
    model.update();

    // create assignment constraints
    for (int j = 0; j < inst.n; j++){      	
		if(inst.items[j][2] == 0) continue; 
		model.addConstr(assigned[j] == inst.items[j][2]); 						
	}

    // create load constraints 
    for (int i = 0; i < DB; i++){
		model.addConstr(load[i] <= inst.W);
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
    if (model.get(GRB_IntAttr_SolCount) >= 1) { 		
		sol.LB += DB;
        sol.UB += DB;
		sol.opt = 1;
		
        // get bin for each item
        for (int j = 0; j < inst.n; j++) {
			if(inst.items[j][2] == 0) continue; 			
			for (int i = 0; i <= DB - inst.items[j][1]; i++){
				if(inst.NPs[j][i] && inst.RCs[j][i] + inst.contR - EPSILON <= DB){
					for (int k = 0; k < ceil(x[i][j].get(GRB_DoubleAttr_X) - EPSILON); k++) 		
						sol.assignedBin[j].push_back(i);
				}
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

void APTPLP(Instance& inst) {
	
    // create a model
	vector<bool> isHActive(inst.LUB,false);
	
    GRBEnv env = GRBEnv();              	
    GRBModel model = GRBModel(env);         

    // declaration of the variables for the model
    vector<vector<GRBVar> > x(inst.LUB, vector<GRBVar>(inst.n));			
    vector<GRBVar> zvec(inst.LUB);
	GRBLinExpr obj = 0;
	
    // initizalization of the variables for the model
    for (int j = 0; j < inst.n; j++) {  
		if(inst.items[j][2] == 0) continue; 
		for (int i = 0; i <= inst.LUB - inst.items[j][1]; i++) {
			if(inst.NPs[j][i]){
				isHActive[i] = true;
				x[i][j] = model.addVar(0, inst.items[j][2], 0, GRB_CONTINUOUS);
			}
        }
    }

	for (int i = 0; i < inst.LUB; i++) 
		zvec[i] =  model.addVar(0, 1, 0, GRB_CONTINUOUS);
    model.update();

    // create linear expressions
	vector<GRBLinExpr> load(inst.LUB, 0); 
	vector<GRBLinExpr> assigned(inst.n, 0);
	
    for (int j = 0; j < inst.n; j++) { 	
		if(inst.items[j][2] == 0) continue; 
        for (int i = 0; i <= inst.LUB - inst.items[j][1]; i++){
			if(inst.NPs[j][i]){
				for (int ip = i; ip < i+inst.items[j][1]; ip++){ 	
					load[ip] += inst.items[j][0] * x[i][j];
				}
				assigned[j] += x[i][j];
			}
		}
    }   
    model.update();

	for (int i = 0; i < inst.LUB; i++)
		obj += zvec[i];
	model.update();

    // create assignment constraints
    for (int j = 0; j < inst.n; j++){      	
		if(inst.items[j][2] == 0) continue; 
		model.addConstr(assigned[j] == inst.items[j][2]); 						
	}

    // create load constraints 
    for (int i = 0; i < inst.LUB; i++){
		model.addConstr(load[i] <= inst.W * zvec[i]);
	}
    
	// create symmetry breaking constraints
    for (int i = 1; i < inst.LUB; i++){
        model.addConstr(zvec[i] <= zvec[i-1]);
	}
	
    // set the objective: minimize z
    model.setObjective(obj, GRB_MINIMIZE);

    // change some settings
    model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
    model.getEnv().set(GRB_IntParam_Threads, 1);
    model.getEnv().set(GRB_DoubleParam_TimeLimit, 3600);

    // find the optimal solution	
    model.optimize();

	// get reduced costs
	for (int j = 0; j < inst.n; j++) {    
		if(inst.items[j][2] == 0) continue; 		
        for (int i = 0; i <= inst.LUB - inst.items[j][1]; i++){
				if(inst.NPs[j][i]){  
					inst.RCs[j][i] = 0;
				}
            }
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