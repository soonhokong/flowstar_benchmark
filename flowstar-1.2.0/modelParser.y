%{
	/*---
  	Flow*: A Taylor Model Based Flowpipe analyzer.
  	Authors: Xin Chen, Erika Abraham and Sriram Sankaranarayanan.
  	Email: Xin Chen <xin.chen@cs.rwth-aachen.de> if you have questions or comments.
  
  	The code is released as is under the GNU General Public License (GPL). Please consult the file LICENSE.txt for
  	further information.
  	---*/


	#include "modelParser.h"

	extern int yyerror(const char *);
	extern int yyerror(string);
	extern int yylex();
	extern int yyparse();
	bool err;
%}

%union
{
	double dblVal;
	string *identifier;
	vector<Interval> *intVec;
	vector<int> *iVec;
	vector<double> *dVec;
	vector<Monomial> *monoVec;
	vector<Polynomial> *polyVec;
	Monomial *mono;
	Polynomial *poly;
	TaylorModelVec *tmVec;
	Matrix *mat;
	vector<vector<double> > *dVecVec;
	vector<PolynomialConstraint> *vecConstraints;
	ResetMap *resetMap;
	Flowpipe *pFlowpipe;
	TaylorModel *ptm;
	Interval *pint;
	vector<string> *strVec;
	TreeNode *pNode;
}


%token<dblVal> NUM
%token<identifier> IDENT
%token STATEVAR TMVAR TM EQ GEQ LEQ ASSIGN END
%token MODE INIT BELONGSTO
%token POLYODE1 POLYODE2
%token VISUALIZE PARAAGGREG INTAGGREG TMAGGREG
%token OUTPUT
%token CONTINUOUS HYBRID
%token SETTING
%token FIXEDST FIXEDORD ADAPTIVEST ADAPTIVEORD
%token MIN MAX
%token REMEST
%token INTERVAL OCTAGON GRID
%token QRPRECOND IDPRECOND
%token TIME
%token MODES JUMPS INV GUARD RESET START MAXJMPS
%token PRINTON PRINTOFF UNSAFESET
%token CONTINUOUSFLOW HYBRIDFLOW
%token TAYLOR_PICARD TAYLOR_REMAINDER TAYLOR_POLYNOMIAL
%token EXP SIN COS LOG SQRT
%token NPODE_TAYLOR CUTOFF PRECISION
%token GNUPLOT MATLAB COMPUTATIONPATHS


%type <poly> polynomial
%type <poly> ODEpolynomial
%type <poly> interval_polynomial
%type <tmVec> ode
%type <iVec> orders
%type <pFlowpipe> init
%type <resetMap> reset
%type <dVec> real_valued_vector
%type <dVecVec> real_valued_vectors
%type <dVec> vector_components
%type <vecConstraints> polynomial_constraints
%type <tmVec> taylor_model
%type <tmVec> interval_taylor_model
%type <intVec> taylor_model_domain
%type <intVec> intervals
%type <intVec> remainders
%type <ptm> non_polynomial_rhs_picard
%type <pint> non_polynomial_rhs_remainder
%type <poly> non_polynomial_rhs_no_remainder
%type <identifier> non_polynomial_rhs_string
%type <strVec> npode
%type <pNode> computation_path


%left GEQ LEQ EQ 
%left '+' '-'
%left '*' '/'
%nonassoc uminus
%right '^'

%start model

%%

model: CONTINUOUS '{' continuous '}'
{
	int mkres = mkdir(outputDir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for output files.\n");
		exit(1);
	}

	mkres = mkdir(imageDir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for images.\n");
		exit(1);
	}

	clock_t begin, end;
	begin = clock();
	continuousProblem.run();
	end = clock();
	printf("%ld flowpipes computed.\n", continuousProblem.numOfFlowpipes());
	printf("time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);
	
	continuousProblem.bSafetyChecking = false;

	printf("Preparing for plotting and dumping...\n");
	continuousProblem.composition();
	printf("Done.\n");

	continuousProblem.plot_2D();

	char filename[NAME_SIZE+10];
	sprintf(filename, "%s%s.flow", outputDir, continuousProblem.outputFileName);
	FILE *fpDumping = fopen(filename, "w");

	if(fpDumping == NULL)
	{
		printf("Can not create the dumping file.\n");
		exit(1);
	}
	
	printf("Dumping the Taylor model flowpipes...\n");
	continuousProblem.dump(fpDumping);
	printf("Done.\n");

	fclose(fpDumping);
}
|
CONTINUOUS '{' continuous '}' unsafe_continuous
{
	int mkres = mkdir(outputDir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for output files.\n");
		exit(1);
	}

	mkres = mkdir(imageDir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for images.\n");
		exit(1);
	}

	clock_t begin, end;
	begin = clock();
	continuousProblem.run();
	end = clock();
	printf("%ld flowpipes computed.\n", continuousProblem.numOfFlowpipes());
	printf("time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);

	continuousProblem.bSafetyChecking = true;

	printf("Preparing for plotting and dumping...\n");
	continuousProblem.composition();
	printf("Done.\n");

	continuousProblem.plot_2D();

	char filename[NAME_SIZE+10];
	sprintf(filename, "%s%s.flow", outputDir, continuousProblem.outputFileName);
	FILE *fpDumping = fopen(filename, "w");

	if(fpDumping == NULL)
	{
		printf("Can not create the dumping file.\n");
		exit(1);
	}
	
	printf("Dumping the Taylor model flowpipes...\n");
	continuousProblem.dump(fpDumping);
	printf("Done.\n");

	fclose(fpDumping);

	int checkingResult;
	printf("Safety checking ...\n");
	begin = clock();
	checkingResult = continuousProblem.safetyChecking();
	end = clock();
	printf("Done.\n");
	printf("time cost for safety checking: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);
	printf("Result: ");

	switch(checkingResult)
	{
	case UNSAFE:
		printf("UNSAFE\n");
		break;
	case SAFE:
		printf("SAFE\n");
		break;
	case UNKNOWN:
		printf("UNKNOWN\n");
		break;
	}
}
|
HYBRID '{' hybrid '}'
{
	int mkres = mkdir(outputDir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for output files.\n");
		exit(1);
	}

	mkres = mkdir(imageDir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for images.\n");
		exit(1);
	}

	clock_t begin, end;
	begin = clock();
	hybridProblem.run();
	end = clock();
	printf("%ld flowpipes computed.\n", hybridProblem.numOfFlowpipes());
	printf("time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);

	hybridProblem.bSafetyChecking = false;

	printf("Preparing for plotting and dumping...\n");
	printf("Done.\n");

	hybridProblem.plot_2D();

	char filename[NAME_SIZE+10];
	sprintf(filename, "%s%s.flow", outputDir, hybridProblem.outputFileName);
	FILE *fpDumping = fopen(filename, "w");

	if(fpDumping == NULL)
	{
		printf("Can not create the dumping file.\n");
		exit(1);
	}

	printf("Dumping the Taylor model flowpipes...\n");
	hybridProblem.dump(fpDumping);
	printf("Done.\n");

	fclose(fpDumping);
}
|
HYBRID '{' hybrid '}' unsafe_hybrid
{
	int mkres = mkdir(outputDir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for output files.\n");
		exit(1);
	}

	mkres = mkdir(imageDir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for images.\n");
		exit(1);
	}

	clock_t begin, end;
	begin = clock();
	hybridProblem.run();
	end = clock();
	printf("%ld flowpipes computed.\n", hybridProblem.numOfFlowpipes());
	printf("time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);

	hybridProblem.bSafetyChecking = true;

	printf("Preparing for plotting and dumping...\n");
	printf("Done.\n");

	hybridProblem.plot_2D();

	char filename[NAME_SIZE+10];
	sprintf(filename, "%s%s.flow", outputDir, hybridProblem.outputFileName);
	FILE *fpDumping = fopen(filename, "w");

	if(fpDumping == NULL)
	{
		printf("Can not create the dumping file.\n");
		exit(1);
	}
	
	printf("Dumping the Taylor model flowpipes...\n");
	hybridProblem.dump(fpDumping);
	printf("Done.\n");

	fclose(fpDumping);

	int checkingResult;
	printf("Safety checking ...\n");
	begin = clock();
	checkingResult = hybridProblem.safetyChecking();
	end = clock();
	printf("Done.\n");
	printf("time cost for safety checking: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);
	printf("Result: ");

	switch(checkingResult)
	{
	case UNSAFE:
		printf("UNSAFE\n");
		break;
	case SAFE:
		printf("SAFE\n");
		break;
	case UNKNOWN:
		printf("UNKNOWN\n");
		break;
	}
}
|
stateVarDecls plotting OUTPUT IDENT unsafe_continuous CONTINUOUSFLOW '{' tmVarDecls continuous_flowpipes '}'
{
	clock_t begin, end;
	strcpy(continuousProblem.outputFileName, $4->c_str());

	continuousProblem.plot_2D();

	int checkingResult;
	printf("Safety checking ...\n");
	begin = clock();
	checkingResult = continuousProblem.safetyChecking();
	end = clock();
	printf("Done.\n");
	printf("time cost for safety checking: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);
	printf("Result: ");

	switch(checkingResult)
	{
	case UNSAFE:
		printf("UNSAFE\n");
		break;
	case SAFE:
		printf("SAFE\n");
		break;
	case UNKNOWN:
		printf("UNKNOWN\n");
		break;
	}
}
|
stateVarDecls plotting OUTPUT IDENT CONTINUOUSFLOW '{' tmVarDecls continuous_flowpipes '}'
{
	strcpy(continuousProblem.outputFileName, $4->c_str());
	continuousProblem.plot_2D();
}
|
stateVarDecls modeDecls COMPUTATIONPATHS '{' computation_paths '}' plotting OUTPUT IDENT unsafe_hybrid HYBRIDFLOW '{' hybrid_flowpipes '}'
{
	clock_t begin, end;
	strcpy(hybridProblem.outputFileName, $9->c_str());
	generateNodeSeq(hybridProblem.traceNodes, hybridProblem.traceTree);

	hybridProblem.plot_2D();

	int checkingResult;
	printf("Safety checking ...\n");
	begin = clock();
	checkingResult = hybridProblem.safetyChecking();
	end = clock();
	printf("Done.\n");
	printf("time cost for safety checking: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);
	printf("Result: ");

	switch(checkingResult)
	{
	case UNSAFE:
		printf("UNSAFE\n");
		break;
	case SAFE:
		printf("SAFE\n");
		break;
	case UNKNOWN:
		printf("UNKNOWN\n");
		break;
	}

	delete $9;
}
|
stateVarDecls modeDecls COMPUTATIONPATHS '{' computation_paths '}' plotting OUTPUT IDENT HYBRIDFLOW '{' hybrid_flowpipes '}'
{
	strcpy(hybridProblem.outputFileName, $9->c_str());
	generateNodeSeq(hybridProblem.traceNodes, hybridProblem.traceTree);
	hybridProblem.plot_2D();

	delete $9;
}
|
TAYLOR_PICARD '{' non_polynomial_rhs_picard '}'
{
	$3->getExpansion(parseResult.expansion);
	parseResult.remainder = $3->getRemainder();
	delete $3;
}
|
TAYLOR_REMAINDER '{' non_polynomial_rhs_remainder '}'
{
	parseResult.remainder = (*$3);
	delete $3;
}
|
TAYLOR_POLYNOMIAL '{' non_polynomial_rhs_no_remainder '}'
{
	parseResult.expansion = (*$3);
	delete $3;
}
;

continuous_flowpipes: continuous_flowpipes '{' interval_taylor_model taylor_model_domain '}'
{
	continuousProblem.flowpipesCompo.push_back(*$3);
	continuousProblem.domains.push_back(*$4);

	delete $3;
	delete $4;
}
|
'{' interval_taylor_model taylor_model_domain '}'
{
	continuousProblem.flowpipesCompo.push_back(*$2);
	continuousProblem.domains.push_back(*$3);

	delete $2;
	delete $3;
}
;

modeDecls: modeDecls IDENT '{' polynomial_constraints '}'
{
	TaylorModelVec tmvDummy;
	vector<Interval> intVecDummy;
	hybridProblem.declareMode(*$2, tmvDummy, intVecDummy, *$4, 0);

	delete $2;
	delete $4;
}
|
IDENT '{' polynomial_constraints '}'
{
	TaylorModelVec tmvDummy;
	vector<Interval> intVecDummy;
	hybridProblem.declareMode(*$1, tmvDummy, intVecDummy, *$3, 0);

	delete $1;
	delete $3;
}
;

hybrid_flowpipes: hybrid_flowpipes IDENT '{' tmVarDecls continuous_flowpipes '}'
{
	int id = hybridProblem.getIDForMode(*$2);
	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	hybridProblem.modeIDs.push_back(id);
	hybridProblem.flowpipesCompo.push_back(continuousProblem.flowpipesCompo);
	hybridProblem.domains.push_back(continuousProblem.domains);

	continuousProblem.flowpipesCompo.clear();
	continuousProblem.domains.clear();
	continuousProblem.tmVarTab.clear();
	continuousProblem.tmVarNames.clear();
}
|
IDENT '{' tmVarDecls continuous_flowpipes '}'
{
	int id = hybridProblem.getIDForMode(*$1);
	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	hybridProblem.modeIDs.push_back(id);
	hybridProblem.flowpipesCompo.push_back(continuousProblem.flowpipesCompo);
	hybridProblem.domains.push_back(continuousProblem.domains);

	continuousProblem.flowpipesCompo.clear();
	continuousProblem.domains.clear();
	continuousProblem.tmVarTab.clear();
	continuousProblem.tmVarNames.clear();
}
;

computation_paths: computation_paths computation_path ';'
{
}
|
computation_path ';'
{
}
;

computation_path: computation_path '(' NUM ',' '[' NUM ',' NUM ']' ')' '-' '>' IDENT
{
	int id = hybridProblem.getIDForMode(*$13);
	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*$13).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	list<TreeNode *>::iterator iter = $$->children.begin();
	bool found = false;
	for(; iter!=$$->children.end(); ++iter)
	{
		if((*iter)->jumpID == $3 && (*iter)->modeID == id)
		{
			$$ = *iter;
			found = true;
			break;
		}
	}

	if(!found)
	{
		Interval I($6, $8);
		TreeNode *tmp = new TreeNode((int)$3, id, I);
		tmp->parent = $$;
		$$->children.push_back(tmp);
		$$ = tmp;
	}

	delete $13;
}
|
IDENT
{
	int id = hybridProblem.getIDForMode(*$1);
	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if(hybridProblem.traceTree == NULL)
	{
		Interval intZero;
		hybridProblem.traceTree = new TreeNode(0, id, intZero);
		$$ = hybridProblem.traceTree;
	}
	else
	{
		if(hybridProblem.traceTree->modeID == id)
		{
			$$ = hybridProblem.traceTree;
		}
		else
		{
			parseError("Invalid computation path.", lineNum);
			exit(1);
		}
	}

	delete $1;
}
;

print: PRINTON
{
	continuousProblem.bPrint = true;
	hybridProblem.bPrint = true;
}
|
PRINTOFF
{
	continuousProblem.bPrint = false;
	hybridProblem.bPrint = false;
}
;

unsafe_continuous: UNSAFESET '{' polynomial_constraints '}'
{
	continuousProblem.unsafeSet = *$3;
	delete $3;
}
;

unsafe_hybrid: UNSAFESET '{' hybrid_constraints '}'
{
}
;

hybrid_constraints: hybrid_constraints IDENT '{' polynomial_constraints '}'
{
	int id = hybridProblem.getIDForMode(*$2);
	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	hybridProblem.unsafeSet[id] = *$4;
	hybridProblem.bVecUnderCheck[id] = true;
	delete $4;
}
|
{
	vector<PolynomialConstraint> vecEmpty;
	for(int i=0; i<hybridProblem.modeNames.size(); ++i)
	{
		hybridProblem.unsafeSet.push_back(vecEmpty);
		hybridProblem.bVecUnderCheck.push_back(false);
	}
}
;

polynomial_constraints: polynomial_constraints ODEpolynomial LEQ NUM
{
	if($2->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	$$ = $1;
	Interval B($4);
	PolynomialConstraint pc(*$2, B);
	$$->push_back(pc);

	delete $2;
}
|
polynomial_constraints ODEpolynomial GEQ NUM
{
	if($2->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	$$ = $1;
	Interval I(-1);
	$2->mul_assign(I);

	Interval B(-$4);
	PolynomialConstraint pc(*$2, B);
	$$->push_back(pc);

	delete $2;
}
|
polynomial_constraints ODEpolynomial EQ NUM
{
	if($2->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	$$ = $1;
	Interval B($4);
	PolynomialConstraint pc1(*$2, B);
	$$->push_back(pc1);

	Interval I(-1);
	$2->mul_assign(I);
	Interval mB(-$4);
	PolynomialConstraint pc2(*$2, mB);
	$$->push_back(pc2);

	delete $2;
}
|
polynomial_constraints ODEpolynomial BELONGSTO '[' NUM ',' NUM ']'
{
	if($2->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	PolynomialConstraint pc1(*$2, $7);
	$$->push_back(pc1);

	Interval I(-1);
	$2->mul_assign(I);
	PolynomialConstraint pc2(*$2, -$5);
	$$->push_back(pc2);

	delete $2;
}
|
{
	$$ = new vector<PolynomialConstraint>(0);
}
;

continuous: stateVarDecls SETTING '{' settings print '}' POLYODE1 '{' ode '}' INIT '{' init '}'
{
	ContinuousSystem system(*$9, gUncertainties, *$13);
	continuousProblem.system = system;
	continuousProblem.integrationScheme = LOW_DEGREE;

	delete $9;
	delete $13;
}
|
stateVarDecls SETTING '{' settings print '}' POLYODE2 '{' ode '}' INIT '{' init '}'
{
	ContinuousSystem system(*$9, gUncertainties, *$13);
	continuousProblem.system = system;
	continuousProblem.integrationScheme = HIGH_DEGREE;

	delete $9;
	delete $13;
}
|
stateVarDecls SETTING '{' settings print '}' NPODE_TAYLOR '{' npode '}' INIT '{' init '}'
{
	ContinuousSystem system(*$9, gUncertainties, *$13);
	continuousProblem.system = system;
	continuousProblem.integrationScheme = NONPOLY_TAYLOR;

	delete $9;
	delete $13;
}
;

hybrid: stateVarDecls SETTING '{' settings MAXJMPS NUM print '}' MODES '{' modes '}' JUMPS '{' jumps '}' INIT '{' hybrid_init '}'
{
	if($6 < 0)
	{
		parseError("The maximum jump depth should be a nonnegative integer.", lineNum);
		exit(1);
	}

	hybridProblem.maxJumps = (int)$6;
}
;

hybrid_init: IDENT '{' intervals '}'
{
	int id = hybridProblem.getIDForMode(*$1);
	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	Interval intZero;
	Flowpipe initialSet(*$3, intZero);
	hybridProblem.initialConfig(id, initialSet);

	int numVars = hybridProblem.stateVarNames.size();

	string tVar("local_t");
	hybridProblem.declareTMVar(tVar);
	continuousProblem.declareTMVar(tVar);

	char name[NAME_SIZE];

	for(int i=0; i<numVars; ++i)
	{
		sprintf(name, "%s%d", local_var_name, i+1);
		string tmVarName(name);
		hybridProblem.declareTMVar(tmVarName);
		continuousProblem.declareTMVar(tmVarName);
	}

	delete $1;
	delete $3;
}
|
IDENT '{' tmVarDecls taylor_model taylor_model_domain '}'
{
	int id = hybridProblem.getIDForMode(*$1);
	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	Flowpipe initialSet(*$4, *$5);
	hybridProblem.initialConfig(id, initialSet);

	delete $4;
	delete $5;
}
;

modes: modes IDENT '{' POLYODE1 '{' ode '}' INV '{' polynomial_constraints '}' '}'
{
	hybridProblem.declareMode(*$2, *$6, gUncertainties, *$10, LOW_DEGREE);

	delete $2;
	delete $6;
	delete $10;
}
|
IDENT '{' POLYODE1 '{' ode '}' INV '{' polynomial_constraints '}' '}'
{
	hybridProblem.declareMode(*$1, *$5, gUncertainties, *$9, LOW_DEGREE);

	delete $1;
	delete $5;
	delete $9;
}
|
modes IDENT '{' POLYODE2 '{' ode '}' INV '{' polynomial_constraints '}' '}'
{
	hybridProblem.declareMode(*$2, *$6, gUncertainties, *$10, HIGH_DEGREE);

	delete $2;
	delete $6;
	delete $10;
}
|
IDENT '{' POLYODE2 '{' ode '}' INV '{' polynomial_constraints '}' '}'
{
	hybridProblem.declareMode(*$1, *$5, gUncertainties, *$9, HIGH_DEGREE);

	delete $1;
	delete $5;
	delete $9;
}
|
modes IDENT '{' NPODE_TAYLOR '{' npode '}' INV '{' polynomial_constraints '}' '}'
{
	hybridProblem.declareMode(*$2, *$6, gUncertainties, *$10, NONPOLY_TAYLOR);

	delete $2;
	delete $6;
	delete $10;
}
|
IDENT '{' NPODE_TAYLOR '{' npode '}' INV '{' polynomial_constraints '}' '}'
{
	hybridProblem.declareMode(*$1, *$5, gUncertainties, *$9, NONPOLY_TAYLOR);

	delete $1;
	delete $5;
	delete $9;
}
;

jumps: jumps IDENT '-' '>' IDENT GUARD '{' polynomial_constraints '}' RESET '{' reset '}' PARAAGGREG '{' real_valued_vectors '}'
{
	int startID = hybridProblem.getIDForMode(*$2);
	if(startID < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	int endID = hybridProblem.getIDForMode(*$5);
	if(endID < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*$5).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if($16->size() > 0)
	{
		hybridProblem.declareTrans(startID, endID, *$8, *$12, PARA_AGGREG, *$16);
	}
	else
	{
		vector<vector<double> > emptyVec;
		hybridProblem.declareTrans(startID, endID, *$8, *$12, PARA_AGGREG, emptyVec);
	}
}
|
jumps IDENT '-' '>' IDENT GUARD '{' polynomial_constraints '}' RESET '{' reset '}' INTAGGREG
{
	int startID = hybridProblem.getIDForMode(*$2);
	if(startID < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	int endID = hybridProblem.getIDForMode(*$5);
	if(endID < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*$5).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	vector<vector<double> > empty;
	hybridProblem.declareTrans(startID, endID, *$8, *$12, INTERVAL_AGGREG, empty);
}
|
{
	hybridProblem.declareTrans();
}
;

reset: reset IDENT '\'' ASSIGN ODEpolynomial '+' '[' NUM ',' NUM ']'
{
	$$ = $1;

	int id = hybridProblem.getIDForStateVar(*$2);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if($8 > $10)
	{
		parseError("Invalid remainder interval.", lineNum);
		exit(1);
	}

	Interval I($8, $10);
	TaylorModel tmTemp(*$5, I);
	$$->tmvReset.tms[id] = tmTemp;

	delete $5;
}
|
reset IDENT '\'' ASSIGN ODEpolynomial
{
	$$ = $1;

	int id = hybridProblem.getIDForStateVar(*$2);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	Interval intZero;
	TaylorModel tmTemp(*$5, intZero);
	$$->tmvReset.tms[id] = tmTemp;

	delete $5;
}
|
{
	int numVars = hybridProblem.stateVarNames.size();

	Matrix coefficients_identity_reset(numVars, numVars+1);

	for(int i=0; i<numVars; ++i)
	{
		coefficients_identity_reset.set(1, i, i+1);
	}

	TaylorModelVec tmvReset(coefficients_identity_reset);

	$$ = new ResetMap(tmvReset);
}
;

real_valued_vectors: real_valued_vectors real_valued_vector
{
	$$->push_back(*$2);
	delete $2;
}
|
{
	$$ = new vector<vector<double> >(0);
}
;

real_valued_vector: '[' vector_components ']'
{
	int rangeDim = $2->size();

	if(rangeDim != hybridProblem.stateVarNames.size())
	{
		parseError("The vector dimension should be equivalent to the system dimension.", lineNum);
		exit(1);
	}

	$$ = new vector<double>(0);

	for(int i=0; i<rangeDim; ++i)
	{
		$$->push_back(0);
	}

	bool bZero = true;
	for(int i=0; i<rangeDim; ++i)
	{
		if((*$2)[i] < -THRESHOLD_LOW || (*$2)[i] > THRESHOLD_LOW)
		{
			if(bZero)
			{
				bZero = false;
			}
		}

		(*$$)[i] = (*$2)[i];
	}

	if(bZero)
	{
		parseError("A template vector should not be zero.", lineNum);
		exit(1);
	}

	delete $2;
}
;

vector_components: vector_components ',' IDENT ':' NUM
{
	int id = hybridProblem.getIDForStateVar(*$3);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	$$ = $1;
	(*$$)[id] = $5;
	delete $3;
}
|
IDENT ':' NUM
{
	int num = hybridProblem.stateVarNames.size();
	$$ = new vector<double>(num);

	int id = hybridProblem.getIDForStateVar(*$1);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	(*$$)[id] = $3;
	delete $1;
}
;

stateVarDecls: STATEVAR stateIdDeclList
{
}
;

stateIdDeclList: stateIdDeclList ',' IDENT
{
	if(!continuousProblem.declareStateVar(*$3))
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s has already been declared.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	hybridProblem.declareStateVar(*$3);
	delete $3;
}
|
IDENT
{
	if(!continuousProblem.declareStateVar(*$1))
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s has already been declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	hybridProblem.declareStateVar(*$1);
	delete $1;
}
;

settings: FIXEDST NUM TIME NUM remainder_estimation precondition plotting FIXEDORD NUM CUTOFF NUM PRECISION NUM OUTPUT IDENT
{
	int order = (int)$9;

	if(order <= 0)
	{
		parseError("Orders should be larger than zero.", lineNum);
		exit(1);
	}

	continuousProblem.bAdaptiveSteps = false;
	continuousProblem.step = $2;
	continuousProblem.time = $4;
	continuousProblem.bAdaptiveOrders = false;
	continuousProblem.orderType = UNIFORM;
	continuousProblem.orders.push_back(order);
	continuousProblem.globalMaxOrder = order;

	hybridProblem.bAdaptiveSteps = false;
	hybridProblem.step = $2;
	hybridProblem.time = $4;
	hybridProblem.bAdaptiveOrders = false;
	hybridProblem.orderType = UNIFORM;
	hybridProblem.orders.push_back(order);
	hybridProblem.globalMaxOrder = order;

	if($11 <= 0)
	{
		parseError("The cutoff threshold should be a positive number.", lineNum);
		exit(1);
	}

	cutoff_threshold = $11;
	intervalNumPrecision = (int)$13;

	strcpy(continuousProblem.outputFileName, (*$15).c_str());
	strcpy(hybridProblem.outputFileName, (*$15).c_str());

	delete $15;
}
|
FIXEDST NUM TIME NUM remainder_estimation precondition plotting ADAPTIVEORD '{' MIN NUM ',' MAX NUM '}' CUTOFF NUM PRECISION NUM OUTPUT IDENT
{
	int minOrder = (int)$11;
	int maxOrder = (int)$14;

	if(minOrder <= 0 || maxOrder <= 0)
	{
		parseError("Orders should be larger than zero.", lineNum);
		exit(1);
	}

	if(minOrder > maxOrder)
	{
		parseError("MAX order should be no smaller than MIN order.", lineNum);
		exit(1);
	}

	continuousProblem.bAdaptiveSteps = false;
	continuousProblem.step = $2;
	continuousProblem.time = $4;
	continuousProblem.bAdaptiveOrders = true;
	continuousProblem.orderType = UNIFORM;
	continuousProblem.orders.push_back(minOrder);
	continuousProblem.maxOrders.push_back(maxOrder);
	continuousProblem.globalMaxOrder = maxOrder;

	hybridProblem.bAdaptiveSteps = false;
	hybridProblem.step = $2;
	hybridProblem.time = $4;
	hybridProblem.bAdaptiveOrders = true;
	hybridProblem.orderType = UNIFORM;
	hybridProblem.orders.push_back(minOrder);
	hybridProblem.maxOrders.push_back(maxOrder);
	hybridProblem.globalMaxOrder = maxOrder;

	if($17 <= 0)
	{
		parseError("The cutoff threshold should be a positive number.", lineNum);
		exit(1);
	}

	cutoff_threshold = $17;
	intervalNumPrecision = (int)$19;

	strcpy(continuousProblem.outputFileName, (*$21).c_str());
	strcpy(hybridProblem.outputFileName, (*$21).c_str());

	delete $21;
}
|
FIXEDST NUM TIME NUM remainder_estimation precondition plotting FIXEDORD '{' orders '}' CUTOFF NUM PRECISION NUM OUTPUT IDENT
{
	continuousProblem.bAdaptiveSteps = false;
	continuousProblem.step = $2;
	continuousProblem.time = $4;
	continuousProblem.bAdaptiveOrders = false;
	continuousProblem.orderType = MULTI;
	continuousProblem.orders = *$10;

	hybridProblem.bAdaptiveSteps = false;
	hybridProblem.step = $2;
	hybridProblem.time = $4;
	hybridProblem.bAdaptiveOrders = false;
	hybridProblem.orderType = MULTI;
	hybridProblem.orders = *$10;

	for(int i=0; i<$10->size(); ++i)
	{
		if((*$10)[i] <= 0)
		{
			parseError("Orders should be larger than zero.", lineNum);
			exit(1);
		}
	}

	int maxOrder = (*$10)[0];
	for(int i=1; i<$10->size(); ++i)
	{
		if(maxOrder < (*$10)[i])
		{
			maxOrder = (*$10)[i];
		}
	}

	continuousProblem.globalMaxOrder = maxOrder;
	hybridProblem.globalMaxOrder = maxOrder;

	if($13 <= 0)
	{
		parseError("The cutoff threshold should be a positive number.", lineNum);
		exit(1);
	}

	cutoff_threshold = $13;
	intervalNumPrecision = (int)$15;

	strcpy(continuousProblem.outputFileName, (*$17).c_str());
	strcpy(hybridProblem.outputFileName, (*$17).c_str());

	delete $10;
	delete $17;
}
|
FIXEDST NUM TIME NUM remainder_estimation precondition plotting ADAPTIVEORD '{' MIN '{' orders '}' ',' MAX '{' orders '}' '}' CUTOFF NUM PRECISION NUM OUTPUT IDENT
{
	continuousProblem.bAdaptiveSteps = false;
	continuousProblem.step = $2;
	continuousProblem.time = $4;
	continuousProblem.bAdaptiveOrders = true;
	continuousProblem.orderType = MULTI;
	continuousProblem.orders = *$12;
	continuousProblem.maxOrders = *$17;

	hybridProblem.bAdaptiveSteps = false;
	hybridProblem.step = $2;
	hybridProblem.time = $4;
	hybridProblem.bAdaptiveOrders = true;
	hybridProblem.orderType = MULTI;
	hybridProblem.orders = *$12;
	hybridProblem.maxOrders = *$17;

	if($12->size() != $17->size())
	{
		parseError("Orders are not properly specified.", lineNum);
		exit(1);
	}

	for(int i=0; i<$17->size(); ++i)
	{
		if((*$12)[i] <= 0 || (*$17)[i] <= 0)
		{
			parseError("Orders should be larger than zero.", lineNum);
			exit(1);
		}

		if((*$12)[i] > (*$17)[i])
		{
			parseError("MAX order should be no smaller than MIN order.", lineNum);
			exit(1);
		}
	}

	int maxOrder = (*$17)[0];
	for(int i=1; i<$17->size(); ++i)
	{
		if(maxOrder < (*$17)[i])
		{
			maxOrder = (*$17)[i];
		}
	}

	continuousProblem.globalMaxOrder = maxOrder;
	hybridProblem.globalMaxOrder = maxOrder;

	if($21 <= 0)
	{
		parseError("The cutoff threshold should be a positive number.", lineNum);
		exit(1);
	}

	cutoff_threshold = $21;
	intervalNumPrecision = (int)$23;

	strcpy(continuousProblem.outputFileName, (*$25).c_str());
	strcpy(hybridProblem.outputFileName, (*$25).c_str());

	delete $12;
	delete $17;
	delete $25;
}
|
ADAPTIVEST '{' MIN NUM ',' MAX NUM '}' TIME NUM remainder_estimation precondition plotting FIXEDORD NUM CUTOFF NUM PRECISION NUM OUTPUT IDENT
{
	if($4 > $7)
	{
		parseError("MIN step should be no larger than MAX step.", lineNum);
		exit(1);
	}

	int order = (int)$15;

	if(order <= 0)
	{
		parseError("Orders should be larger than zero.", lineNum);
		exit(1);
	}

	continuousProblem.bAdaptiveSteps = true;
	continuousProblem.step = $7;
	continuousProblem.miniStep = $4;
	continuousProblem.time = $10;
	continuousProblem.bAdaptiveOrders = false;
	continuousProblem.orderType = UNIFORM;
	continuousProblem.orders.push_back(order);
	continuousProblem.globalMaxOrder = order;

	hybridProblem.bAdaptiveSteps = true;
	hybridProblem.step = $7;
	hybridProblem.miniStep = $4;
	hybridProblem.time = $10;
	hybridProblem.bAdaptiveOrders = false;
	hybridProblem.orderType = UNIFORM;
	hybridProblem.orders.push_back(order);
	hybridProblem.globalMaxOrder = order;

	if($17 <= 0)
	{
		parseError("The cutoff threshold should be a positive number.", lineNum);
		exit(1);
	}

	cutoff_threshold = $17;
	intervalNumPrecision = (int)$19;

	strcpy(continuousProblem.outputFileName, (*$21).c_str());
	strcpy(hybridProblem.outputFileName, (*$21).c_str());

	delete $21;
}
|
ADAPTIVEST '{' MIN NUM ',' MAX NUM '}' TIME NUM remainder_estimation precondition plotting FIXEDORD '{' orders '}' CUTOFF NUM PRECISION NUM OUTPUT IDENT
{
	if($4 > $7)
	{
		parseError("MIN step should be no larger than MAX step.", lineNum);
		exit(1);
	}

	for(int i=0; i<$16->size(); ++i)
	{
		if((*$16)[i] <= 0)
		{
			parseError("Orders should be larger than zero.", lineNum);
			exit(1);
		}
	}

	continuousProblem.bAdaptiveSteps = true;
	continuousProblem.step = $7;
	continuousProblem.miniStep = $4;
	continuousProblem.time = $10;
	continuousProblem.bAdaptiveOrders = false;
	continuousProblem.orderType = MULTI;
	continuousProblem.orders = *$16;

	hybridProblem.bAdaptiveSteps = true;
	hybridProblem.step = $7;
	hybridProblem.miniStep = $4;
	hybridProblem.time = $10;
	hybridProblem.bAdaptiveOrders = false;
	hybridProblem.orderType = MULTI;
	hybridProblem.orders = *$16;

	int maxOrder = (*$16)[0];
	for(int i=1; i<$16->size(); ++i)
	{
		if(maxOrder < (*$16)[i])
		{
			maxOrder = (*$16)[i];
		}
	}
	
	continuousProblem.globalMaxOrder = maxOrder;
	hybridProblem.globalMaxOrder = maxOrder;

	if($19 <= 0)
	{
		parseError("The cutoff threshold should be a positive number.", lineNum);
		exit(1);
	}

	cutoff_threshold = $19;
	intervalNumPrecision = (int)$21;

	strcpy(continuousProblem.outputFileName, (*$23).c_str());
	strcpy(hybridProblem.outputFileName, (*$23).c_str());

	delete $16;
	delete $23;
}
;

remainder_estimation: REMEST NUM
{
	if($2 <= 0)
	{
		parseError("Remainder estimation should be a positive number.", lineNum);
		exit(1);
	}

	Interval I(-$2, $2);

	for(int i=0; i<continuousProblem.stateVarNames.size(); ++i)
	{
		continuousProblem.estimation.push_back(I);
		hybridProblem.estimation.push_back(I);
	}
}
|
REMEST '{' remainders '}'
{
	for(int i=0; i<$3->size(); ++i)
	{
		if((*$3)[i].inf() >= (*$3)[i].sup() - THRESHOLD_LOW)
		{
			parseError("Invalid remainder estimation.", lineNum);
			exit(1);
		}
	}

	continuousProblem.estimation = *$3;
	hybridProblem.estimation = *$3;
	delete $3;
}
;

remainders: remainders ',' IDENT ':' '[' NUM ',' NUM ']'
{
	$$ = $1;
	int id = continuousProblem.getIDForStateVar(*$3);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if($6 >= $8)
	{
		parseError("Invalid remainder estimation.", lineNum);
		exit(1);
	}

	Interval I($6,$8);
	(*$$)[id] = I;
	delete $3;
}
|
IDENT ':' '[' NUM ',' NUM ']'
{
	int numVars = continuousProblem.stateVarNames.size();
	$$ = new vector<Interval>(numVars);

	int id = continuousProblem.getIDForStateVar(*$1);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if($4 >= $6)
	{
		parseError("Invalid remainder estimation.", lineNum);
		exit(1);
	}

	Interval I($4,$6);
	(*$$)[id] = I;
	delete $1;
}
;

orders: orders ',' IDENT ':' NUM
{
	$$ = $1;
	int id = continuousProblem.getIDForStateVar(*$3);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	(*$$)[id] = (int)$5;
	delete $3;
}
|
IDENT ':' NUM
{
	int numVars = continuousProblem.stateVarNames.size();
	$$ = new vector<int>(numVars);
	for(int i=0; i<numVars; ++i)
	{
		(*$$)[i] = 0;
	}

	int id = continuousProblem.getIDForStateVar(*$1);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	(*$$)[id] = (int)$3;
	delete $1;
}
;

precondition: QRPRECOND
{
	continuousProblem.precondition = QR_PRE;
	hybridProblem.precondition = QR_PRE;
}
|
IDPRECOND
{
	continuousProblem.precondition = ID_PRE;
	hybridProblem.precondition = ID_PRE;
}
;

plotting: GNUPLOT INTERVAL IDENT ',' IDENT
{
	int x = continuousProblem.getIDForStateVar(*$3);
	int y = continuousProblem.getIDForStateVar(*$5);

	if(x < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if(y < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$5).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	continuousProblem.outputAxes.push_back(x);
	continuousProblem.outputAxes.push_back(y);
	continuousProblem.plotSetting = PLOT_INTERVAL;
	continuousProblem.plotFormat = PLOT_GNUPLOT;

	hybridProblem.outputAxes.push_back(x);
	hybridProblem.outputAxes.push_back(y);
	hybridProblem.plotSetting = PLOT_INTERVAL;
	hybridProblem.plotFormat = PLOT_GNUPLOT;

	delete $3;
	delete $5;
}
|
GNUPLOT OCTAGON IDENT ',' IDENT
{
	int x = continuousProblem.getIDForStateVar(*$3);
	int y = continuousProblem.getIDForStateVar(*$5);

	if(x < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State Variable %s is not declared.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if(y < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$5).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	continuousProblem.outputAxes.push_back(x);
	continuousProblem.outputAxes.push_back(y);
	continuousProblem.plotSetting = PLOT_OCTAGON;
	continuousProblem.plotFormat = PLOT_GNUPLOT;

	hybridProblem.outputAxes.push_back(x);
	hybridProblem.outputAxes.push_back(y);
	hybridProblem.plotSetting = PLOT_OCTAGON;
	hybridProblem.plotFormat = PLOT_GNUPLOT;

	delete $3;
	delete $5;
}
|
GNUPLOT GRID NUM IDENT ',' IDENT
{
	int x = continuousProblem.getIDForStateVar(*$4);
	int y = continuousProblem.getIDForStateVar(*$6);

	if(x < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$4).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if(y < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$6).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	continuousProblem.outputAxes.push_back(x);
	continuousProblem.outputAxes.push_back(y);
	continuousProblem.plotSetting = PLOT_GRID;
	continuousProblem.numSections = (int)$3;
	continuousProblem.plotFormat = PLOT_GNUPLOT;

	hybridProblem.outputAxes.push_back(x);
	hybridProblem.outputAxes.push_back(y);
	hybridProblem.plotSetting = PLOT_GRID;
	hybridProblem.numSections = (int)$3;
	hybridProblem.plotFormat = PLOT_GNUPLOT;

	delete $4;
	delete $6;
}
|
MATLAB INTERVAL IDENT ',' IDENT
{
	int x = continuousProblem.getIDForStateVar(*$3);
	int y = continuousProblem.getIDForStateVar(*$5);

	if(x < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if(y < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$5).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	continuousProblem.outputAxes.push_back(x);
	continuousProblem.outputAxes.push_back(y);
	continuousProblem.plotSetting = PLOT_INTERVAL;
	continuousProblem.plotFormat = PLOT_MATLAB;

	hybridProblem.outputAxes.push_back(x);
	hybridProblem.outputAxes.push_back(y);
	hybridProblem.plotSetting = PLOT_INTERVAL;
	hybridProblem.plotFormat = PLOT_MATLAB;

	delete $3;
	delete $5;
}
|
MATLAB OCTAGON IDENT ',' IDENT
{
	int x = continuousProblem.getIDForStateVar(*$3);
	int y = continuousProblem.getIDForStateVar(*$5);

	if(x < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State Variable %s is not declared.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if(y < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$5).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	continuousProblem.outputAxes.push_back(x);
	continuousProblem.outputAxes.push_back(y);
	continuousProblem.plotSetting = PLOT_OCTAGON;
	continuousProblem.plotFormat = PLOT_MATLAB;

	hybridProblem.outputAxes.push_back(x);
	hybridProblem.outputAxes.push_back(y);
	hybridProblem.plotSetting = PLOT_OCTAGON;
	hybridProblem.plotFormat = PLOT_MATLAB;

	delete $3;
	delete $5;
}
|
MATLAB GRID NUM IDENT ',' IDENT
{
	int x = continuousProblem.getIDForStateVar(*$4);
	int y = continuousProblem.getIDForStateVar(*$6);

	if(x < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$4).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if(y < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$6).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	continuousProblem.outputAxes.push_back(x);
	continuousProblem.outputAxes.push_back(y);
	continuousProblem.plotSetting = PLOT_GRID;
	continuousProblem.numSections = (int)$3;
	continuousProblem.plotFormat = PLOT_MATLAB;

	hybridProblem.outputAxes.push_back(x);
	hybridProblem.outputAxes.push_back(y);
	hybridProblem.plotSetting = PLOT_GRID;
	hybridProblem.numSections = (int)$3;
	hybridProblem.plotFormat = PLOT_MATLAB;

	delete $4;
	delete $6;
}
;

init: tmVarDecls taylor_model taylor_model_domain
{
	$$ = new Flowpipe(*$2, *$3);

	delete $2;
	delete $3;
}
|
intervals
{
	Interval intZero;
	$$ = new Flowpipe(*$1, intZero);

	delete $1;
}

tmVarDecls: TMVAR tmIdDeclList
{
}
;

tmIdDeclList: tmIdDeclList ',' IDENT
{
	if(!continuousProblem.declareTMVar(*$3))
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "TM variable %s has already been declared.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	hybridProblem.declareTMVar(*$3);
	delete $3;
}
|
IDENT
{
	string tVar("local_t");
	continuousProblem.declareTMVar(tVar);
	hybridProblem.declareTMVar(tVar);

	if(!continuousProblem.declareTMVar(*$1))
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "TM variable %s has already been declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	hybridProblem.declareTMVar(*$1);
	delete $1;
}
;

taylor_model: taylor_model IDENT EQ polynomial '+' '[' NUM ',' NUM ']'
{
	int id = continuousProblem.getIDForStateVar(*$2);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if($7 > $9)
	{
		parseError("Invalid interval remainder.", lineNum);
		exit(1);
	}

	Interval I($7,$9);
	TaylorModel tmTemp(*$4, I);
	$$ = $1;
	$$->tms[id] = tmTemp;

	delete $2;
	delete $4;
}
|
{
	TaylorModel tmEmpty;
	$$ = new TaylorModelVec;

	for(int i=0; i<continuousProblem.stateVarNames.size(); ++i)
	{
		$$->tms.push_back(tmEmpty);
	}
}
;

taylor_model_domain: taylor_model_domain IDENT BELONGSTO '[' NUM ',' NUM ']'
{
	int id = continuousProblem.getIDForTMVar(*$2);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "TM variable %s is not declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if($5 > $7)
	{
		parseError("Invalid interval.", lineNum);
		exit(1);
	}

	Interval I($5,$7);
	$$ = $1;
	(*$$)[id] = I;

	delete $2;
}
|
IDENT BELONGSTO '[' NUM ',' NUM ']'
{
	$$ = new vector<Interval>( continuousProblem.tmVarNames.size() );

	Interval intZero;
	(*$$)[0] = intZero;

	int id = continuousProblem.getIDForTMVar(*$1);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "TM variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if($4 > $6)
	{
		parseError("Invalid interval.", lineNum);
		exit(1);
	}

	Interval I($4,$6);
	(*$$)[id] = I;

	delete $1;
}
;

intervals: intervals IDENT BELONGSTO '[' NUM ',' NUM ']'
{
	int id = continuousProblem.getIDForStateVar(*$2);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if($5 > $7)
	{
		parseError("Invalid interval.", lineNum);
		exit(1);
	}

	Interval I($5,$7);
	$$ = $1;
	(*$$)[id] = I;

	delete $2;
}
|
{
	int numVars = continuousProblem.stateVarNames.size();
	$$ = new vector<Interval>(numVars);

	string tVar("local_t");
	continuousProblem.declareTMVar(tVar);

	char name[NAME_SIZE];

	for(int i=0; i<numVars; ++i)
	{
		sprintf(name, "%s%d", local_var_name, i+1);
		string tmVarName(name);
		continuousProblem.declareTMVar(tmVarName);
	}
}
;

ode: ode IDENT '\'' EQ ODEpolynomial
{
	$$ = $1;

	int id = continuousProblem.getIDForStateVar(*$2);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	Interval intZero;
	TaylorModel tmTemp(*$5, intZero);
	$$->tms[id] = tmTemp;

	delete $2;
	delete $5;
}
|
ode IDENT '\'' EQ ODEpolynomial '+' '[' NUM ',' NUM ']'
{
	$$ = $1;

	int id = continuousProblem.getIDForStateVar(*$2);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if($8 > $10)
	{
		parseError("Invalid interval", lineNum);
		exit(1);
	}

	Interval uncertainty($8, $10);
	gUncertainties[id] = uncertainty;

	Interval intZero;
	TaylorModel tmTemp(*$5, intZero);
	$$->tms[id] = tmTemp;

	delete $2;
	delete $5;
}
|
{
	int numVars = continuousProblem.stateVarNames.size();

	$$ = new TaylorModelVec;
	TaylorModel tmTemp;
	Interval intZero;

	gUncertainties.clear();

	for(int i=0; i<numVars; ++i)
	{
		$$->tms.push_back(tmTemp);
		gUncertainties.push_back(intZero);
	}
}
;

npode: npode IDENT '\'' EQ non_polynomial_rhs_string
{
	$$ = $1;

	int id = continuousProblem.getIDForStateVar(*$2);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	(*$$)[id] = (*$5);

	delete $2;
	delete $5;
}
|
npode IDENT '\'' EQ non_polynomial_rhs_string '+' '[' NUM ',' NUM ']'
{
	$$ = $1;

	int id = continuousProblem.getIDForStateVar(*$2);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	(*$$)[id] = (*$5);

	if($8 > $10)
	{
		parseError("Invalid interval", lineNum);
		exit(1);
	}

	Interval uncertainty($8, $10);
	gUncertainties[id] = uncertainty;

	delete $2;
	delete $5;
}
|
{
	int numVars = continuousProblem.stateVarNames.size();
	$$ = new vector<string>;

	string empty;
	Interval intZero;

	gUncertainties.clear();

	for(int i=0; i<numVars; ++i)
	{
		$$->push_back(empty);
		gUncertainties.push_back(intZero);
	}
}
;

polynomial: polynomial '+' polynomial
{
	$$ = $1;
	(*$$) += (*$3);

	delete $3;
}
|
polynomial '-' polynomial
{
	$$ = $1;
	(*$$) -= (*$3);

	delete $3;
}
|
'(' polynomial ')'
{
	$$ = $2; 
}
|
polynomial '*' polynomial
{
	$$ = $1;
	(*$$) *= (*$3);

	delete $3;
}
|
polynomial '^' NUM
{
	int exp = (int) $3;

	if(exp == 0)
	{
		Interval I(1);
		$$ = new Polynomial(I, continuousProblem.tmVarNames.size());
	}
	else
	{
		$$ = new Polynomial(*$1);

		for(int i=1; i<exp; ++i)
		{
			(*$$) *= (*$1);
		}
	}

	delete $1;
}
|
'-' polynomial %prec uminus
{
	Interval I(-1);
	$$ = $2;
	$$->mul_assign(I);
}
|
IDENT
{
	int id = continuousProblem.getIDForTMVar(*$1);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "TM variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	int numVars = continuousProblem.tmVarNames.size();
	Interval I(1);

	vector<int> degrees;
	for(int i=0; i<numVars; ++i)
	{
		degrees.push_back(0);
	}

	degrees[id] = 1;
	Monomial monomial(I, degrees);

	$$ = new Polynomial(monomial);
	delete $1;
}
|
NUM
{
	int numVars = continuousProblem.tmVarNames.size();
	Interval I($1);
	$$ = new Polynomial(I, numVars);
}
;

ODEpolynomial: ODEpolynomial '+' ODEpolynomial
{
	$$ = $1;
	(*$$) += (*$3);

	delete $3;
}
|
ODEpolynomial '-' ODEpolynomial
{
	$$ = $1;
	(*$$) -= (*$3);

	delete $3;
}
|
'(' ODEpolynomial ')'
{
	$$ = $2; 
}
|
ODEpolynomial '*' ODEpolynomial
{
	$$ = $1;
	(*$$) *= (*$3);

	delete $3;
}
|
ODEpolynomial '^' NUM
{
	int exp = (int) $3;

	if(exp == 0)
	{
		Interval I(1);
		$$ = new Polynomial(I, continuousProblem.stateVarNames.size()+1);
	}
	else
	{
		$$ = new Polynomial(*$1);

		for(int i=1; i<exp; ++i)
		{
			(*$$) *= (*$1);
		}
	}

	delete $1;
}
|
'-' ODEpolynomial %prec uminus
{
	Interval I(-1);
	$$ = $2;
	$$->mul_assign(I);
}
|
IDENT
{
	int id = continuousProblem.getIDForStateVar(*$1);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	int numVars = continuousProblem.stateVarNames.size()+1;
	Interval I(1);

	vector<int> degrees;
	for(int i=0; i<numVars; ++i)
	{
		degrees.push_back(0);
	}

	degrees[id+1] = 1;
	Monomial monomial(I, degrees);

	$$ = new Polynomial(monomial);
	delete $1;
}
|
NUM
{
	int numVars = continuousProblem.stateVarNames.size()+1;
	Interval I($1);
	$$ = new Polynomial(I, numVars);
}
;

interval_taylor_model: interval_taylor_model IDENT EQ interval_polynomial '+' '[' NUM ',' NUM ']'
{
	int id = continuousProblem.getIDForStateVar(*$2);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if($7 > $9)
	{
		parseError("Invalid interval remainder.", lineNum);
		exit(1);
	}

	Interval I($7,$9);
	TaylorModel tmTemp(*$4, I);
	$$ = $1;
	$$->tms[id] = tmTemp;

	delete $2;
	delete $4;
}
|
IDENT EQ interval_polynomial '+' '[' NUM ',' NUM ']'
{
	TaylorModel tmEmpty;
	$$ = new TaylorModelVec;

	for(int i=0; i<continuousProblem.stateVarNames.size(); ++i)
	{
		$$->tms.push_back(tmEmpty);
	}

	int id = continuousProblem.getIDForStateVar(*$1);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if($6 > $8)
	{
		parseError("Invalid interval remainder.", lineNum);
		exit(1);
	}

	Interval I($6,$8);
	TaylorModel tmTemp(*$3, I);

	$$->tms[id] = tmTemp;

	delete $1;
	delete $3;
}
;

interval_polynomial: interval_polynomial '+' interval_polynomial
{
	$$ = $1;
	(*$$) += (*$3);

	delete $3;
}
|
interval_polynomial '-' interval_polynomial
{
	$$ = $1;
	(*$$) -= (*$3);

	delete $3;
}
|
'(' interval_polynomial ')'
{
	$$ = $2; 
}
|
interval_polynomial '*' interval_polynomial
{
	$$ = $1;
	(*$$) *= (*$3);

	delete $3;
}
|
interval_polynomial '^' NUM
{
	int exp = (int) $3;

	if(exp == 0)
	{
		Interval I(1);
		$$ = new Polynomial(I, continuousProblem.tmVarNames.size());
	}
	else
	{
		$$ = new Polynomial(*$1);

		for(int i=1; i<exp; ++i)
		{
			(*$$) *= (*$1);
		}
	}

	delete $1;
}
|
'-' interval_polynomial %prec uminus
{
	Interval I(-1);
	$$ = $2;
	$$->mul_assign(I);
}
|
IDENT
{
	int id = continuousProblem.getIDForTMVar(*$1);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "TM variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	int numVars = continuousProblem.tmVarNames.size();
	Interval I(1);

	vector<int> degrees;
	for(int i=0; i<numVars; ++i)
	{
		degrees.push_back(0);
	}

	degrees[id] = 1;
	Monomial monomial(I, degrees);

	$$ = new Polynomial(monomial);
	delete $1;
}
|
'[' NUM ',' NUM ']'
{
	int numVars = continuousProblem.tmVarNames.size();
	Interval I($2, $4);
	$$ = new Polynomial(I, numVars);
}
;




















non_polynomial_rhs_picard: non_polynomial_rhs_picard '+' non_polynomial_rhs_picard
{
	$$ = $1;
	$1->add_assign(*$3);
	delete $3;
}
|
non_polynomial_rhs_picard '-' non_polynomial_rhs_picard
{
	$$ = $1;
	$1->sub_assign(*$3);
	delete $3;
}
|
non_polynomial_rhs_picard '*' non_polynomial_rhs_picard
{
	$$ = $1;

	Interval intPoly1, intPoly2, intTrunc;

	$3->polyRangeNormal(intPoly2, parseSetting.step_exp_table);
	$1->mul_insert_ctrunc_normal_assign(intPoly1, intTrunc, *$3, intPoly2, parseSetting.step_exp_table, parseSetting.order);

	parseSetting.ranges.push_back(intPoly1);
	parseSetting.ranges.push_back(intPoly2);
	parseSetting.ranges.push_back(intTrunc);

	delete $3;
}
|
'(' non_polynomial_rhs_picard ')'
{
	$$ = $2;
}
|
non_polynomial_rhs_picard '/' non_polynomial_rhs_picard
{
	TaylorModel tmTemp;
	$3->rec_taylor(tmTemp, parseSetting.ranges, parseSetting.step_exp_table, continuousProblem.stateVarNames.size()+1, parseSetting.order);

	$$ = $1;

	Interval intPoly1, intPoly2, intTrunc;

	tmTemp.polyRangeNormal(intPoly2, parseSetting.step_exp_table);
	$1->mul_insert_ctrunc_normal_assign(intPoly1, intTrunc, tmTemp, intPoly2, parseSetting.step_exp_table, parseSetting.order);

	parseSetting.ranges.push_back(intPoly1);
	parseSetting.ranges.push_back(intPoly2);
	parseSetting.ranges.push_back(intTrunc);

	delete $3;
}
|
EXP '(' non_polynomial_rhs_picard ')'
{
	TaylorModel tmTemp;
	$3->exp_taylor(tmTemp, parseSetting.ranges, parseSetting.step_exp_table, continuousProblem.stateVarNames.size()+1, parseSetting.order);

	*$3 = tmTemp;
	$$ = $3;
}
|
SIN '(' non_polynomial_rhs_picard ')'
{
	TaylorModel tmTemp;
	$3->sin_taylor(tmTemp, parseSetting.ranges, parseSetting.step_exp_table, continuousProblem.stateVarNames.size()+1, parseSetting.order);

	*$3 = tmTemp;
	$$ = $3;
}
|
COS '(' non_polynomial_rhs_picard ')'
{
	TaylorModel tmTemp;
	$3->cos_taylor(tmTemp, parseSetting.ranges, parseSetting.step_exp_table, continuousProblem.stateVarNames.size()+1, parseSetting.order);

	*$3 = tmTemp;
	$$ = $3;
}
|
LOG '(' non_polynomial_rhs_picard ')'
{
	TaylorModel tmTemp;
	$3->log_taylor(tmTemp, parseSetting.ranges, parseSetting.step_exp_table, continuousProblem.stateVarNames.size()+1, parseSetting.order);

	*$3 = tmTemp;
	$$ = $3;
}
|
SQRT '(' non_polynomial_rhs_picard ')'
{
	TaylorModel tmTemp;
	$3->sqrt_taylor(tmTemp, parseSetting.ranges, parseSetting.step_exp_table, continuousProblem.stateVarNames.size()+1, parseSetting.order);

	*$3 = tmTemp;
	$$ = $3;
}
|
non_polynomial_rhs_picard '^' NUM
{
	int exp = (int)$3;

	if(exp == 0)
	{
		Interval I(1);
		$$ = new TaylorModel(I, continuousProblem.stateVarNames.size()+1);
	}
	else
	{
		$$ = new TaylorModel(*$1);

		Interval intPoly1, intPoly2, intTrunc;
		$1->polyRangeNormal(intPoly2, parseSetting.step_exp_table);

		for(int i=2; i<=exp; ++i)
		{
			$$->mul_insert_ctrunc_normal_assign(intPoly1, intTrunc, *$1, intPoly2, parseSetting.step_exp_table, parseSetting.order);

			parseSetting.ranges.push_back(intPoly1);
			parseSetting.ranges.push_back(intPoly2);
			parseSetting.ranges.push_back(intTrunc);
		}
	}

	delete $1;
}
|
'-' non_polynomial_rhs_picard %prec uminus
{
	Interval I(-1);
	$$ = $2;
	$$->mul_assign(I);
}
|
IDENT
{
	int id = continuousProblem.getIDForStateVar(*$1);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	$$ = new TaylorModel;
	(*$$) = parseSetting.flowpipe.tms[id];

	delete $1;
}
|
'[' NUM ',' NUM ']'
{
	int numVars = continuousProblem.stateVarNames.size()+1;
	Interval I($2, $4);
	$$ = new TaylorModel(I, numVars);
}
;
















non_polynomial_rhs_remainder: non_polynomial_rhs_remainder '+' non_polynomial_rhs_remainder
{
	$$ = $1;
	(*$$) += (*$3);
	delete $3;
}
|
non_polynomial_rhs_remainder '-' non_polynomial_rhs_remainder
{
	$$ = $1;
	(*$$) -= (*$3);
	delete $3;
}
|
non_polynomial_rhs_remainder '*' non_polynomial_rhs_remainder
{
	$$ = new Interval;

	(*$$) = (*parseSetting.iterRange) * (*$3);
	++parseSetting.iterRange;
	(*$$) += (*parseSetting.iterRange) * (*$1);
	(*$$) += (*$1) * (*$3);
	++parseSetting.iterRange;
	(*$$) += (*parseSetting.iterRange);
	++parseSetting.iterRange;

	delete $1;
	delete $3;
}
|
'(' non_polynomial_rhs_remainder ')'
{
	$$ = $2;
}
|
non_polynomial_rhs_remainder '/' non_polynomial_rhs_remainder
{
	Interval intTemp;
	rec_taylor_only_remainder(intTemp, *$3, parseSetting.iterRange, parseSetting.order);

	$$ = new Interval;

	(*$$) = (*parseSetting.iterRange) * intTemp;
	++parseSetting.iterRange;
	(*$$) += (*parseSetting.iterRange) * (*$1);
	(*$$) += (*$1) * intTemp;
	++parseSetting.iterRange;
	(*$$) += (*parseSetting.iterRange);
	++parseSetting.iterRange;

	delete $1;
	delete $3;
}
|
EXP '(' non_polynomial_rhs_remainder ')'
{
	Interval intTemp;
	exp_taylor_only_remainder(intTemp, *$3, parseSetting.iterRange, parseSetting.order);

	(*$3) = intTemp;
	$$ = $3;
}
|
SIN '(' non_polynomial_rhs_remainder ')'
{
	Interval intTemp;
	sin_taylor_only_remainder(intTemp, *$3, parseSetting.iterRange, parseSetting.order);

	(*$3) = intTemp;
	$$ = $3;
}
|
COS '(' non_polynomial_rhs_remainder ')'
{
	Interval intTemp;
	cos_taylor_only_remainder(intTemp, *$3, parseSetting.iterRange, parseSetting.order);

	(*$3) = intTemp;
	$$ = $3;
}
|
LOG '(' non_polynomial_rhs_remainder ')'
{
	Interval intTemp;
	log_taylor_only_remainder(intTemp, *$3, parseSetting.iterRange, parseSetting.order);

	(*$3) = intTemp;
	$$ = $3;
}
|
SQRT '(' non_polynomial_rhs_remainder ')'
{
	Interval intTemp;
	sqrt_taylor_only_remainder(intTemp, *$3, parseSetting.iterRange, parseSetting.order);

	(*$3) = intTemp;
	$$ = $3;
}
|
non_polynomial_rhs_remainder '^' NUM
{
	int exp = (int)$3;

	if(exp == 0)
	{
		Interval intZero;
		(*$1) = intZero;
		$$ = $1;
	}
	else
	{
		$$ = new Interval(*$1);

		for(int i=2; i<=exp; ++i)
		{
			Interval intTemp;
			intTemp = (*parseSetting.iterRange) * (*$1);
			++parseSetting.iterRange;
			intTemp += (*parseSetting.iterRange) * (*$$);
			intTemp += (*$1) * (*$$);
			++parseSetting.iterRange;
			intTemp += (*parseSetting.iterRange);
			++parseSetting.iterRange;
			
			(*$$) = intTemp;
		}
	}

	delete $1;
}
|
'-' non_polynomial_rhs_remainder %prec uminus
{
	$$ = $2;
	$$->inv_assign();
}
|
IDENT
{
	int id = continuousProblem.getIDForStateVar(*$1);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	$$ = new Interval;
	(*$$) = parseSetting.flowpipe.tms[id].getRemainder();

	delete $1;
}
|
'[' NUM ',' NUM ']'
{
	$$ = new Interval;
}
;

















non_polynomial_rhs_no_remainder: non_polynomial_rhs_no_remainder '+' non_polynomial_rhs_no_remainder
{
	$$ = $1;
	(*$$) += (*$3);

	delete $3;
}
|
non_polynomial_rhs_no_remainder '-' non_polynomial_rhs_no_remainder
{
	$$ = $1;
	(*$$) -= (*$3);

	delete $3;
}
|
non_polynomial_rhs_no_remainder '*' non_polynomial_rhs_no_remainder
{
	$$ = $1;
	(*$$) *= (*$3);
	$$->nctrunc(parseSetting.order);
	$$->cutoff();

	delete $3;
}
|
'(' non_polynomial_rhs_no_remainder ')'
{
	$$ = $2;
}
|
non_polynomial_rhs_no_remainder '/' non_polynomial_rhs_no_remainder
{
	Polynomial polyTemp;
	$3->rec_taylor(polyTemp, continuousProblem.stateVarNames.size()+1, parseSetting.order);

	(*$1) *= polyTemp;
	$1->nctrunc(parseSetting.order);
	$$ = $1;
	$$->cutoff();

	delete $3;
}
|
EXP '(' non_polynomial_rhs_no_remainder ')'
{
	Polynomial polyTemp;
	$3->exp_taylor(polyTemp, continuousProblem.stateVarNames.size()+1, parseSetting.order);

	(*$3) = polyTemp;
	$$ = $3;
}
|
SIN '(' non_polynomial_rhs_no_remainder ')'
{
	Polynomial polyTemp;
	$3->sin_taylor(polyTemp, continuousProblem.stateVarNames.size()+1, parseSetting.order);

	(*$3) = polyTemp;
	$$ = $3;
}
|
COS '(' non_polynomial_rhs_no_remainder ')'
{
	Polynomial polyTemp;
	$3->cos_taylor(polyTemp, continuousProblem.stateVarNames.size()+1, parseSetting.order);

	(*$3) = polyTemp;
	$$ = $3;
}
|
LOG '(' non_polynomial_rhs_no_remainder ')'
{
	Polynomial polyTemp;
	$3->log_taylor(polyTemp, continuousProblem.stateVarNames.size()+1, parseSetting.order);

	(*$3) = polyTemp;
	$$ = $3;
}
|
SQRT '(' non_polynomial_rhs_no_remainder ')'
{
	Polynomial polyTemp;
	$3->sqrt_taylor(polyTemp, continuousProblem.stateVarNames.size()+1, parseSetting.order);

	(*$3) = polyTemp;
	$$ = $3;
}
|
non_polynomial_rhs_no_remainder '^' NUM
{
	Interval I(1);
	Polynomial polyTemp(I, continuousProblem.stateVarNames.size()+1);

	int degree = (int)$3;
	for(int i=0; i<degree; ++i)
	{
		polyTemp *= *$1;
		polyTemp.nctrunc(parseSetting.order);
	}

	$$ = $1;
	(*$$) = polyTemp;
	$$->cutoff();
}
|
'-' non_polynomial_rhs_no_remainder %prec uminus
{
	$$ = $2;
	$$->inv_assign();
}
|
IDENT
{
	int id = continuousProblem.getIDForStateVar(*$1);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	$$ = new Polynomial;
	parseSetting.flowpipe.tms[id].getExpansion(*$$);
}
|
'[' NUM ',' NUM ']'
{
	Interval I($2, $4);
	$$ = new Polynomial(I, continuousProblem.stateVarNames.size()+1);
}
;














non_polynomial_rhs_string: non_polynomial_rhs_string '+' non_polynomial_rhs_string
{
	(*$1) += ' ';
	(*$1) += '+';
	(*$1) += ' ';
	(*$1) += (*$3);

	$$ = $1;
	delete $3;
}
|
non_polynomial_rhs_string '-' non_polynomial_rhs_string
{
	(*$1) += ' ';
	(*$1) += '-';
	(*$1) += ' ';
	(*$1) += (*$3);

	$$ = $1;
	delete $3;
}
|
non_polynomial_rhs_string '*' non_polynomial_rhs_string
{
	(*$1) += ' ';
	(*$1) += '*';
	(*$1) += ' ';
	(*$1) += (*$3);

	$$ = $1;
	delete $3;
}
|
'(' non_polynomial_rhs_string ')'
{
	string str;
	str += '(';
	str += (*$2);
	str += ')';
	(*$2) = str;

	$$ = $2;
}
|
non_polynomial_rhs_string '/' non_polynomial_rhs_string
{
	(*$1) += ' ';
	(*$1) += '/';
	(*$1) += ' ';
	(*$1) += (*$3);

	$$ = $1;
	delete $3;
}
|
EXP '(' non_polynomial_rhs_string ')'
{
	string str("exp");
	str += '(';
	str += (*$3);
	str += ')';
	(*$3) = str;

	$$ = $3;
}
|
SIN '(' non_polynomial_rhs_string ')'
{
	string str("sin");
	str += '(';
	str += (*$3);
	str += ')';
	(*$3) = str;

	$$ = $3;
}
|
COS '(' non_polynomial_rhs_string ')'
{
	string str("cos");
	str += '(';
	str += (*$3);
	str += ')';
	(*$3) = str;

	$$ = $3;
}
|
LOG '(' non_polynomial_rhs_string ')'
{
	string str("log");
	str += '(';
	str += (*$3);
	str += ')';
	(*$3) = str;

	$$ = $3;
}
|
SQRT '(' non_polynomial_rhs_string ')'
{
	string str("sqrt");
	str += '(';
	str += (*$3);
	str += ')';
	(*$3) = str;

	$$ = $3;
}
|
non_polynomial_rhs_string '^' NUM
{
	(*$1) += '^';

	char strNum[NUM_LENGTH];
	sprintf(strNum, "%d", (int)$3);
	string num(strNum);
	(*$1) += num;

	$$ = $1;
}
|
'-' non_polynomial_rhs_string %prec uminus
{
	string str;
	str += '-';
	str += (*$2);
	(*$2) = str;

	$$ = $2;
}
|
IDENT
{
	int id = continuousProblem.getIDForStateVar(*$1);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	$$ = $1;
}
|
NUM
{
	$$ = new string;
	char strNum[NUM_LENGTH];
	sprintf(strNum, "%.20e", $1);
	string num(strNum);

	(*$$) += '[';
	(*$$) += num;
	(*$$) += ' ';
	(*$$) += ',';
	(*$$) += ' ';
	(*$$) += num;
	(*$$) += ']';
}
;
















%%

int yyerror(const char * what)
{
	fprintf(stderr, "Error line %d: %s\n", lineNum,what);
	err=true;
	return 1;
}

int yyerror(string what)
{
	cerr << "Error line "<<lineNum<<" "<<what<<endl;
	err=true;
	return 1;
}
