This repository contains the code for all algorithms discussed in the paper "Solving the parallel processor scheduling and bin packing problems with contiguity constraints: mathematical models and computational studies" by Fatih Burak Akcay and Maxence Delorme. 

Our algorithms are coded in C++ and use the commercial solver Gurobi for the ILP models and Cplex for the CP models. 
The code is divided over 2 main folders (one for the P|cont|Cmax and one for the 1CBP), each containing the subfolders for each method discussed. The different folders correspond to the following methods in our paper:
Note: Preprocessing is included in all models unless explicitly stated.

For P|cont|Cmax main folder:
- 1_BKRS_BIN					      | Model BKRS with binary decision variables
- 1_BKRS_BIN_NP					    | Model BKRS with binary decision variables without preprocessing
- 2_BM_BIN						      | Model BM with binary decision variables
- 2_BM_BIN_MIM 				    	| Model BM with binary decision variables and MIM patterns
- 2_BM_BIN_NOR				    	| Model BM with binary decision variables and normal patterns
- 2_BM_BIN_NP				        | Model BM with binary expansion without preprocessing
- 2_BM_BINEXP_NOR				    | Model BM with binary expansion and normal patterns
- 2_BM_INT_NOR					    | Model BM with integer decision variables and normal patterns
- 2_BM_INT_NOR_DB				    | Model BM with integer decision variables and normal patterns + destructive bounds
- 2_BM_INT_NOR_DBRCVF	    	| Model BM with integer decision variables and normal patterns + destructive bounds and reduced-cost variable fixing
- 2_BM_INT_NOR_SEED			    | Model BM with integer decision variables and normal patterns (different seed values)
- 2_BM_INT_NOR_V1			    	| Model BM with integer decision variables and normal patterns + symmetry breaking constraints of type (59)+(60)
- 2_BM_INT_NOR_V2			    	| Model BM with integer decision variables and normal patterns + symmetry breaking constraints of type (59)+(61)
- 2_BM_INT_NOR_V3			    	| Model BM with integer decision variables and normal patterns + symmetry breaking constraints of type (63)
- 2_BM_INT_NOR_V4			    	| Model BM with integer decision variables and normal patterns + symmetry breaking constraints of type (64) 
- 3_FLOW_BIN					      | Model FLOW-PCC with binary decision variables
- 3_FLOW_BIN_MIM				    | Model FLOW-PCC with binary decision variables and MIM patterns
- 3_FLOW_BIN_NOR				    | Model FLOW-PCC with binary decision variables and normal patterns
- 3_FLOW_BIN_NP					    | Model FLOW-PCC with binary decision variables without preprocessing
- 3_FLOW_BINEXP_NOR			    | Model FLOW-PCC with binary expansion and normal patterns
- 3_FLOW_INT_NOR				    | Model FLOW-PCC with integer decision variables and normal patterns
- 3_FLOW_INT_NOR_DB			    | Model FLOW-PCC with integer decision variables and normal patterns + destructive bounds
- 3_FLOW_INT_NOR_DBRCVF			| Model FLOW-PCC with integer decision variables and normal patterns + destructive bounds and reduced-cost variable fixing
- 3_FLOW_INT_NOR_V1			  	| Model FLOW-PCC with integer decision variables and normal patterns + symmetry breaking constraints of type (59)+(60)
- 3_FLOW_INT_NOR_V2			  	| Model FLOW-PCC with integer decision variables and normal patterns + symmetry breaking constraints of type (59)+(61)
- 3_FLOW_INT_NOR_V3			  	| Model FLOW-PCC with integer decision variables and normal patterns + symmetry breaking constraints of type (63)
- 3_FLOW_INT_NOR_V4			  	| Model FLOW-PCC with integer decision variables and normal patterns + symmetry breaking constraints of type (64) 
- 4_CP_BIN						      | Model CP-PCC
- 4_CP_BIN_MIM				    	| Model CP-PCC with MIM patterns 
- 4_CP_BIN_NOR					    | Model CP-PCC with normal patterns
- 4_CP_BIN_NP					      | Model CP-PCC without preprocessing
- 4_CP_BIN_V1				      	| Model CP-PCC + symmetry breaking constraints of type (59)+(60)
- 4_CP_BIN_V2				      	| Model CP-PCC + symmetry breaking constraints of type (59)+(61)
- 4_CP_BIN_V3					      | Model CP-PCC + symmetry breaking constraints of type (64)
- 4_CP_BIN_V4				      	| Model CP-PCC + symmetry breaking constraints of type (65)
- 4_CP_BIN_V4_DB			    	| Model CP-PCC + symmetry breaking constraints of type (65) and destructive bounds
- 4_CP_BIN_V4_DB_SEED		  	| Model CP-PCC + symmetry breaking constraints of type (65) and destructive bounds (different seed values)
- 5_LB_CO					        	| Model LB_3 (P|conf|Cmax)
- 5_LB_CO_LP				      	| Model LB_3 (P|conf|Cmax) LP relaxation
- 5_LB_CSP					      	| Model LB_2 (P||Cmax)
- 5_LB_CSP_LP					      | Model LB_2 (P||Cmax) LP relaxation
- 5_LB_RE					        	| Model LB_4 (REFLECT)
- 5_LB_RE_LP					      | Model LB_4 (REFLECT) LP relaxation

For 1CBP main folder:
- 1_APT_BIN					      	| Model APT with binary decision variables
- 1_APT_BIN_LB					    | Model APT with binary decision variables (+ a valid lower bound)
- 2_APTP_BIN				      	| Model APTP with binary decision variables
- 2_APTP_BIN_LB				     	| Model APTP with binary decision variables (+ a valid lower bound)
- 2_APTP_BIN_LB_NOR	    		| Model APTP with binary decision variables and normal patterns (+ a valid lower bound)
- 2_APTP_BINEXP_LB_NOR			| Model APTP with binary expansion and normal patterns (+ a valid lower bound)
- 2_APTP_INT_LB_NOR				  | Model APTP with integer decision variables and normal patterns (+ a valid lower bound)
- 2_APTP_INT_LB_NOR_V1			| Model APTP with integer decision variables and normal patterns + symmetry breaking constraints of type (59)+(60) (+ a valid lower bound)
- 2_APTP_INT_LB_NOR_V2			| Model APTP with integer decision variables and normal patterns + symmetry breaking constraints of type (59)+(61) (+ a valid lower bound)
- 2_APTP_INT_LB_NOR_V3			| Model APTP with integer decision variables and normal patterns + symmetry breaking constraints of type (63) (+ a valid lower bound)
- 2_APTP_INT_LB_NOR_V4			| Model APTP with integer decision variables and normal patterns + symmetry breaking constraints of type (64) (+ a valid lower bound)
- 2_APTP_INT_NOR_DB			  	| Model APTP with integer decision variables and normal patterns + destructive bounds (+ a valid lower bound)
- 2_APTP_INT_NOR_DBRCVF			| Model APTP with integer decision variables and normal patterns + destructive bounds and reduced-cost variable fixing
- 3_FLOW_BIN					      | Model FLOW-CBP with binary decision varibles
- 3_FLOW_BIN_LB				     	| Model FLOW-CBP with binary decision variables (+ a valid lower bound)
- 3_FLOW_BIN_LB_NOR				  | Model FLOW-CBP with binary decision variables and normal patterns (+ a valid lower bound)
- 3_FLOW_BINEXP_LB_NOR			| Model FLOW-CBP with binary expansion and normal patterns (+ a valid lower bound)
- 3_FLOW_INT_LB_NOR				  | Model FLOW-CBP with integer decision variables and normal patterns (+ a valid lower bound)
- 3_FLOW_INT_LB_NOR_V1			| Model FLOW-CBP with integer decision variables and normal patterns + symmetry breaking constraints of type (59)+(60) (+ a valid lower bound)
- 3_FLOW_INT_LB_NOR_V2			| Model FLOW-CBP with integer decision variables and normal patterns + symmetry breaking constraints of type (59)+(61) (+ a valid lower bound)
- 3_FLOW_INT_LB_NOR_V3			| Model FLOW-CBP with integer decision variables and normal patterns + symmetry breaking constraints of type (63) (+ a valid lower bound)
- 3_FLOW_INT_LB_NOR_V4			| Model FLOW-CBP with integer decision variables and normal patterns + symmetry breaking constraints of type (64) (+ a valid lower bound)
- 3_FLOW_INT_NOR_DB				  | Model FLOW-CBP with integer decision variables and normal patterns + destructive bounds 
- 3_FLOW_INT_NOR_DBRCVF			| Model FLOW-CBP with integer decision variables and normal patterns + destructive bounds and reduced-cost variable fixing
- 3_FLOW_INT_NOR_DBRCVF_SEED| Model FLOW-CBP with integer decision variables and normal patterns + destructive bounds and reduced-cost variable fixing (different seed values)
- 4_CP_CBP						      | Model CP-CBP
- 4_CP_CBP_LB					      | Model CP-CBP (+ a valid lower bound)
- 4_CP_CBP_LB_NOR			    	| Model CP-CBP with normal patterns (+ a valid lower bound)
- 4_CP_CBP_LB_V1				    | Model CP-PCC + symmetry breaking constraints of type (59)+(60) (+ a valid lower bound)
- 4_CP_CBP_LB_V2				    | Model CP-PCC + symmetry breaking constraints of type (59)+(61) (+ a valid lower bound)
- 4_CP_CBP_LB_V3				    | Model CP-PCC + symmetry breaking constraints of type (64) (+ a valid lower bound)
- 4_CP_CBP_LB_V3_SEED			  | Model CP-PCC + symmetry breaking constraints of type (64) (different seed values) (+ a valid lower bound)
- 4_CP_CBP_LB_V4				    | Model CP-PCC + symmetry breaking constraints of type (65) (+ a valid lower bound)
- 4_CP_CBP_V4_DB			    	| Model CP-PCC + symmetry breaking constraints of type (65) with destructive bounds 
- 5_LB_CO						        | Model LB_3 (arcflow formulation of the one-dimensional cutting stock problem with conflict constraints)
- 5_LB_CO_LP					      | Model LB_3 (arcflow formulation of the one-dimensional cutting stock problem with conflict constraints) LP relaxation
- 5_LB_CSP					      	| Model LB_2 (arcflow formulation of the one-dimensional cutting stock problem)
- 5_LB_CSP_LP					      | Model LB_2 (arcflow formulation of the one-dimensional cutting stock problem) LP relaxation

Each folder contains the same substructure. For example, 1_BKRS_BIN	contains the following files:
- helper_functions.cpp		  | Contains a number of secondary functions (this file is usually the same for each subfolder)
- helper_functions.h		  	| The header file corresponding to helper_functions.cpp (this file is usually the same for each subfolder)
- main.cpp						      | The front-end code for using the method  
- main.h						        | The header file corresponding to main.cpp 
- model.cpp						      | Contains the model of the tested method
- model.h						        | The header file corresponding to model.cpp 
- makefile						      | Used for compiling under linux (it needs to be updated by the user)
Note: for lower bound models, main.cpp (as well as main.h) is named LB.cpp (LB.h)

********
Once compiled, the following command can be used to run the algorithm:
	./PROGRAM "./PATH_INSTANCE" "NAME_INSTANCE" "./PATH_AND_NAME_OUTPUT_GENERAL" 
where
- PROGRAM is the name of the compiled software 
- ./PATH_INSTANCE is the relative path of the folder where the instance to solve is located
- NAME_INSTANCE is the name of the instance to solve
- ./PATH_AND_NAME_OUTPUT_GENERAL is the name of the file (together with its relative path) where performance metrics (such as the optimality status, the CPU time required, or the number of variables) are stored after solving an instance
********

Moreover, "_INPUT.rar" contains a txt-file for each of our test instances. There are 7 main folders, each corresponding to a different dataset:
- BENG		| A set of ten instances proposed by Bengtsson)
- BKW		  | A set of thirteen SCP instances proposed by Burke et al.)
- CGCUT   | A set of three instances proposed by Christofides and Whitlock)
- CLASS		| Ten sets of fifty instances proposed by Berkey and Wang \cite{BW87} (first six classes) and Martello and Vigo \cite{MV98} (last four classes))
- GCUT		| A set of thirteen instances proposed by Beasley)
- HT		  | A set of nine instances proposed by Hopper and Turton)
- NGCUT		| A set of twelve instances proposed by Beasley) 

Each txt-file is structured as follows:
- the first line contains the number of item types (m)
- the second line contains the width of the strip (W)
- the remaining (m) lines all contain, for each item type:
    	- the item dimensions (width and height) and its demand (w_j h_j d_j) 

