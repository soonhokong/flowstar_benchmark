/*---
  Flow*: A Taylor Model Based Flowpipe analyzer.
  Authors: Xin Chen, Erika Abraham and Sriram Sankaranarayanan.
  Email: Xin Chen <xin.chen@cs.rwth-aachen.de> if you have questions or comments.

  The code is released as is under the GNU General Public License (GPL). Please consult the file LICENSE.txt for
  further information.
---*/

#include "Hybrid.h"

// class ResetMap

ResetMap::ResetMap()
{
}

ResetMap::ResetMap(const TaylorModelVec & tmv)
{
	tmvReset = tmv;
}

ResetMap::ResetMap(const ResetMap & reset)
{
	tmvReset = reset.tmvReset;
}

ResetMap::~ResetMap()
{
}

void ResetMap::reset(TaylorModelVec & result, const TaylorModelVec & tmv, const vector<Interval> & domain, const int order) const
{
	vector<Interval> tmvPolyRange;
	tmv.polyRange(tmvPolyRange, domain);
	tmvReset.insert_ctrunc(result, tmv, tmvPolyRange, domain, order);
}

ResetMap & ResetMap::operator = (const ResetMap & reset)
{
	if(this == &reset)
		return *this;

	tmvReset = reset.tmvReset;
	return *this;
}























// class DiscTrans

DiscTrans::DiscTrans()
{
}

DiscTrans::DiscTrans(const int id, const int start, const int target, const vector<PolynomialConstraint> & lcs, const ResetMap & reset)
{
	jumpID = id;
	startID = start;
	targetID = target;
	guard = lcs;
	resetMap = reset;
}

DiscTrans::DiscTrans(const DiscTrans & trans)
{
	jumpID = trans.jumpID;
	startID = trans.startID;
	targetID = trans.targetID;
	guard = trans.guard;
	resetMap = trans.resetMap;
}

DiscTrans::~DiscTrans()
{
	guard.clear();
}

DiscTrans & DiscTrans::operator = (const DiscTrans & trans)
{
	if(this == &trans)
		return *this;

	jumpID = trans.jumpID;
	startID = trans.startID;
	targetID = trans.targetID;
	guard = trans.guard;
	resetMap = trans.resetMap;
	return *this;
}


























// computation tree

TreeNode::TreeNode(const int jump, const int mode, const Interval & t)
{
	jumpID = jump;
	modeID = mode;
	localTime = t;
	parent = NULL;
}

TreeNode::TreeNode(const TreeNode & node)
{
	jumpID = node.jumpID;
	modeID = node.modeID;
	localTime = node.localTime;
	parent = node.parent;
	children = node.children;
}

TreeNode::~TreeNode()
{
	list<TreeNode *>::iterator iter = children.begin();

	for(; iter!=children.end(); ++iter)
	{
		delete *iter;
	}

	children.clear();
}

void TreeNode::dump(FILE *fp, const string & prefix, const vector<string> & modeNames) const
{
	char buffer[NAME_SIZE];

	if(jumpID == 0)
	{
		sprintf(buffer, "%s", modeNames[modeID].c_str());
	}
	else
	{
		string strTime;
		localTime.toString(strTime);

		sprintf(buffer, " ( %d , %s ) -> %s", jumpID, strTime.c_str(), modeNames[modeID].c_str());
	}

	string strTemp(buffer);
	string strPath = prefix + strTemp;

	if(children.size() == 0)
	{
		fprintf(fp, "%s;\n\n", strPath.c_str());
	}
	else
	{
		list<TreeNode *>::const_iterator iter = children.begin();

		for(; iter!=children.end(); ++iter)
		{
			(*iter)->dump(fp, strPath, modeNames);
		}
	}
}

TreeNode & TreeNode::operator = (const TreeNode & node)
{
	if(this == &node)
		return *this;

	jumpID = node.jumpID;
	modeID = node.modeID;
	localTime = node.localTime;
	parent = node.parent;
	children = node.children;
	return *this;
}

































// class HybridSystem

HybridSystem::HybridSystem()
{
}

HybridSystem::HybridSystem(const vector<int> & modes_input, const vector<TaylorModelVec> & odes_input, const vector<vector<HornerForm> > & hfOdes_input,
		const vector<vector<string> > & strOdes_input, const vector<vector<PolynomialConstraint> > & invariants_input, const vector<vector<DiscTrans> > & transitions_input,
		const int initMode,	const vector<vector<Interval> > & uncertainties_input, const vector<vector<Interval> > & uncertainty_centers_input, const Flowpipe & initSet)
{
	modes				=	modes_input;
	odes				=	odes_input;
	hfOdes				=	hfOdes_input;
	strOdes				=	strOdes_input;
	invariants			=	invariants_input;
	transitions			=	transitions_input;
	initialMode			=	initMode;
	initialSet			=	initSet;
	uncertainties		=	uncertainties_input;
	uncertainty_centers	=	uncertainty_centers_input;
}

HybridSystem::HybridSystem(const HybridSystem & hybsys)
{
	modes				=	hybsys.modes;
	odes				=	hybsys.odes;
	hfOdes				=	hybsys.hfOdes;
	strOdes				=	hybsys.strOdes;
	invariants			=	hybsys.invariants;
	transitions			=	hybsys.transitions;
	initialMode			=	hybsys.initialMode;
	initialSet			=	hybsys.initialSet;
	uncertainties		=	hybsys.uncertainties;
	uncertainty_centers	=	hybsys.uncertainty_centers;
}

HybridSystem::~HybridSystem()
{
	modes.clear();
	odes.clear();
	hfOdes.clear();
	strOdes.clear();
	invariants.clear();
	transitions.clear();
	uncertainties.clear();
	uncertainty_centers.clear();
}

HybridSystem & HybridSystem::operator = (const HybridSystem & hybsys)
{
	if(this == &hybsys)
		return *this;

	modes				=	hybsys.modes;
	odes				=	hybsys.odes;
	hfOdes				=	hybsys.hfOdes;
	strOdes				=	hybsys.strOdes;
	invariants			=	hybsys.invariants;
	transitions			=	hybsys.transitions;
	initialMode			=	hybsys.initialMode;
	initialSet			=	hybsys.initialSet;
	uncertainties		=	hybsys.uncertainties;
	uncertainty_centers	=	hybsys.uncertainty_centers;

	return *this;
}

// for low-degree ODEs
// fixed step sizes and orders

bool HybridSystem::reach_continuous_low_degree(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp, const double step,
		const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
		vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;
	Interval intZero;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);

	resultsCompo.clear();

	TaylorModelVec tmvTemp;
	initFp.composition_normal(tmvTemp, step_exp_table);
	resultsCompo.push_back(tmvTemp);
	domains.push_back(initFp.domain);

	Flowpipe newFlowpipe, currentFlowpipe = initFp;

	vector<Polynomial> polyODE;
	for(int i=0; i<odes[mode].tms.size(); ++i)
	{
		polyODE.push_back(odes[mode].tms[i].expansion);
	}

	vector<HornerForm> taylorExpansion;
	computeTaylorExpansion(taylorExpansion, polyODE, order);

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_low_degree(newFlowpipe, hfOdes[mode], taylorExpansion, precondition, step_exp_table, step_end_exp_table, order, estimation, uncertainties[mode]);

		if(bvalid)
		{
			// over-approximate the flowpipe/invariant intersection
			TaylorModelVec tmvCompo;
			newFlowpipe.composition_normal(tmvCompo, step_exp_table);

			vector<Interval> contracted_domain = newFlowpipe.domain;
			vector<bool> local_boundary_intersected;
			int type = contract_interval_arithmetic(tmvCompo, contracted_domain, invariants[mode], local_boundary_intersected);

			if(type >= 0)
			{
				// collect the intersected invariant boundary
				if(invariant_boundary_intersected.size() != local_boundary_intersected.size())
				{
					invariant_boundary_intersected = local_boundary_intersected;
				}
				else
				{
					for(int i=0; i<local_boundary_intersected.size(); ++i)
					{
						if(local_boundary_intersected[i])
						{
							invariant_boundary_intersected[i] = true;
						}
					}
				}
			}

			switch(type)
			{
			case -1:	// invariant violated
				return true;
			case 0:		// domain is not contracted
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 1: 	// time interval is not contracted
				newFlowpipe.normalize();
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 2: 	// time interval is contracted
				if(contracted_domain[0] > intZero)
				{
					return true;
				}
				else
				{
					resultsCompo.push_back(tmvCompo);
					domains.push_back(contracted_domain);
					return true;
				}
			}

			t += step;

			if(bPrint)
			{
				printf("mode: %s,\t", modeNames[mode].c_str());
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("order = %d\n", order);
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			return false;
		}
	}

	return true;
}

bool HybridSystem::reach_continuous_low_degree(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
		const double step, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition,
		const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
		vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;
	Interval intZero;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	resultsCompo.clear();

	TaylorModelVec tmvTemp;
	initFp.composition_normal(tmvTemp, step_exp_table);
	resultsCompo.push_back(tmvTemp);
	domains.push_back(initFp.domain);

	Flowpipe newFlowpipe, currentFlowpipe = initFp;

	vector<Polynomial> polyODE;
	for(int i=0; i<odes[mode].tms.size(); ++i)
	{
		polyODE.push_back(odes[mode].tms[i].expansion);
	}

	vector<HornerForm> taylorExpansion;
	computeTaylorExpansion(taylorExpansion, polyODE, orders);

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_low_degree(newFlowpipe, hfOdes[mode], taylorExpansion, precondition, step_exp_table, step_end_exp_table, orders, estimation, uncertainties[mode]);

		if(bvalid)
		{
			// over-approximate the flowpipe/invariant intersection
			TaylorModelVec tmvCompo;
			newFlowpipe.composition_normal(tmvCompo, step_exp_table);

			vector<Interval> contracted_domain = newFlowpipe.domain;
			vector<bool> local_boundary_intersected;
			int type = contract_interval_arithmetic(tmvCompo, contracted_domain, invariants[mode], local_boundary_intersected);

			if(type >= 0)
			{
				// collect the intersected invariant boundary
				if(invariant_boundary_intersected.size() != local_boundary_intersected.size())
				{
					invariant_boundary_intersected = local_boundary_intersected;
				}
				else
				{
					for(int i=0; i<local_boundary_intersected.size(); ++i)
					{
						if(local_boundary_intersected[i])
						{
							invariant_boundary_intersected[i] = true;
						}
					}
				}
			}

			switch(type)
			{
			case -1:	// invariant violated
				return true;
			case 0:		// domain is not contracted
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 1: 	// time interval is not contracted
				newFlowpipe.normalize();
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 2: 	// time interval is contracted
				if(contracted_domain[0] > intZero)
				{
					return true;
				}
				else
				{
					resultsCompo.push_back(tmvCompo);
					domains.push_back(contracted_domain);
					return true;
				}
			}

			t += step;

			if(bPrint)
			{
				printf("mode: %s,\t", modeNames[mode].c_str());
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("orders:\t");
				int num = orders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), orders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), orders[num]);
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			return false;
		}
	}

	return true;
}

// adaptive step sizes and fixed orders

bool HybridSystem::reach_continuous_low_degree(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
		const double step, const double miniStep, const double time, const int order, const int precondition,
		const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
		vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;
	Interval intZero;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);

	resultsCompo.clear();

	TaylorModelVec tmvTemp;
	initFp.composition_normal(tmvTemp, step_exp_table);
	resultsCompo.push_back(tmvTemp);
	domains.push_back(initFp.domain);

	Flowpipe newFlowpipe, currentFlowpipe = initFp;

	double newStep = 0;

	vector<Polynomial> polyODE;
	for(int i=0; i<odes[mode].tms.size(); ++i)
	{
		polyODE.push_back(odes[mode].tms[i].expansion);
	}

	vector<HornerForm> taylorExpansion;
	computeTaylorExpansion(taylorExpansion, polyODE, order);

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_low_degree(newFlowpipe, hfOdes[mode], taylorExpansion, precondition, step_exp_table, step_end_exp_table, newStep, miniStep, order, estimation, uncertainties[mode]);

		if(bvalid)
		{
			// over-approximate the flowpipe/invariant intersection
			TaylorModelVec tmvCompo;
			newFlowpipe.composition_normal(tmvCompo, step_exp_table);

			vector<Interval> contracted_domain = newFlowpipe.domain;
			vector<bool> local_boundary_intersected;
			int type = contract_interval_arithmetic(tmvCompo, contracted_domain, invariants[mode], local_boundary_intersected);

			if(type >= 0)
			{
				// collect the intersected invariant boundary
				if(invariant_boundary_intersected.size() != local_boundary_intersected.size())
				{
					invariant_boundary_intersected = local_boundary_intersected;
				}
				else
				{
					for(int i=0; i<local_boundary_intersected.size(); ++i)
					{
						if(local_boundary_intersected[i])
						{
							invariant_boundary_intersected[i] = true;
						}
					}
				}
			}

			switch(type)
			{
			case -1:	// invariant violated
				return true;
			case 0:		// domain is not contracted
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 1: 	// time interval is not contracted
				newFlowpipe.normalize();
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 2: 	// time interval is contracted
				if(contracted_domain[0] > intZero)
				{
					return true;
				}
				else
				{
					resultsCompo.push_back(tmvCompo);
					domains.push_back(contracted_domain);
					return true;
				}
			}

			t += step_exp_table[1].sup();

			if(bPrint)
			{
				printf("mode: %s,\t", modeNames[mode].c_str());
				printf("time = %f,\t", t);
				printf("step = %f,\t", step_exp_table[1].sup());
				printf("order = %d\n", order);
			}

			newStep = step_exp_table[1].sup() * LAMBDA_UP;
			if(newStep > step - THRESHOLD_HIGH)
			{
				newStep = 0;
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			return false;
		}
	}

	return true;
}

bool HybridSystem::reach_continuous_low_degree(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
		const double step, const double miniStep, const double time, const vector<int> & orders, const int globalMaxOrder,
		const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
		vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;
	Interval intZero;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	resultsCompo.clear();

	TaylorModelVec tmvTemp;
	initFp.composition_normal(tmvTemp, step_exp_table);
	resultsCompo.push_back(tmvTemp);
	domains.push_back(initFp.domain);

	Flowpipe newFlowpipe, currentFlowpipe = initFp;

	double newStep = 0;

	vector<Polynomial> polyODE;
	for(int i=0; i<odes[mode].tms.size(); ++i)
	{
		polyODE.push_back(odes[mode].tms[i].expansion);
	}

	vector<HornerForm> taylorExpansion;
	computeTaylorExpansion(taylorExpansion, polyODE, orders);

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_low_degree(newFlowpipe, hfOdes[mode], taylorExpansion, precondition, step_exp_table, step_end_exp_table, newStep, miniStep, orders, globalMaxOrder, estimation, uncertainties[mode]);

		if(bvalid)
		{
			// over-approximate the flowpipe/invariant intersection
			TaylorModelVec tmvCompo;
			newFlowpipe.composition_normal(tmvCompo, step_exp_table);

			vector<Interval> contracted_domain = newFlowpipe.domain;
			vector<bool> local_boundary_intersected;
			int type = contract_interval_arithmetic(tmvCompo, contracted_domain, invariants[mode], local_boundary_intersected);

			if(type >= 0)
			{
				// collect the intersected invariant boundary
				if(invariant_boundary_intersected.size() != local_boundary_intersected.size())
				{
					invariant_boundary_intersected = local_boundary_intersected;
				}
				else
				{
					for(int i=0; i<local_boundary_intersected.size(); ++i)
					{
						if(local_boundary_intersected[i])
						{
							invariant_boundary_intersected[i] = true;
						}
					}
				}
			}

			switch(type)
			{
			case -1:	// invariant violated
				return true;
			case 0:		// domain is not contracted
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 1: 	// time interval is not contracted
				newFlowpipe.normalize();
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 2: 	// time interval is contracted
				if(contracted_domain[0] > intZero)
				{
					return true;
				}
				else
				{
					resultsCompo.push_back(tmvCompo);
					domains.push_back(contracted_domain);
					return true;
				}
			}

			t += step_exp_table[1].sup();

			if(bPrint)
			{
				printf("mode: %s,\t", modeNames[mode].c_str());
				printf("time = %f,\t", t);
				printf("step = %f,\t", step_exp_table[1].sup());
				printf("orders:\t");
				int num = orders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), orders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), orders[num]);
			}

			newStep = step_exp_table[1].sup() * LAMBDA_UP;
			if(newStep > step - THRESHOLD_HIGH)
			{
				newStep = 0;
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			return false;
		}
	}

	return true;
}

// adaptive orders and fixed step sizes

bool HybridSystem::reach_continuous_low_degree(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
		const double step, const double time, const int order, const int maxOrder, const int precondition, const vector<Interval> & estimation,
		const bool bPrint, const vector<string> & stateVarNames, vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;
	Interval intZero;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*maxOrder);

	resultsCompo.clear();

	TaylorModelVec tmvTemp;
	initFp.composition_normal(tmvTemp, step_exp_table);
	resultsCompo.push_back(tmvTemp);
	domains.push_back(initFp.domain);

	Flowpipe newFlowpipe, currentFlowpipe = initFp;

	int newOrder = order;
	int localMaxOrder = order;

	vector<Polynomial> polyODE;
	for(int i=0; i<odes[mode].tms.size(); ++i)
	{
		polyODE.push_back(odes[mode].tms[i].expansion);
	}

	vector<HornerForm> taylorExpansion;
	computeTaylorExpansion(taylorExpansion, polyODE, order);

	vector<vector<HornerForm> > expansions;
	expansions.push_back(taylorExpansion);

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_low_degree(newFlowpipe, hfOdes[mode], expansions[newOrder-order], precondition, step_exp_table, step_end_exp_table, newOrder, maxOrder, estimation, uncertainties[mode]);

		if(bvalid)
		{
			// over-approximate the flowpipe/invariant intersection
			TaylorModelVec tmvCompo;
			newFlowpipe.composition_normal(tmvCompo, step_exp_table);

			vector<Interval> contracted_domain = newFlowpipe.domain;
			vector<bool> local_boundary_intersected;
			int type = contract_interval_arithmetic(tmvCompo, contracted_domain, invariants[mode], local_boundary_intersected);

			if(type >= 0)
			{
				// collect the intersected invariant boundary
				if(invariant_boundary_intersected.size() != local_boundary_intersected.size())
				{
					invariant_boundary_intersected = local_boundary_intersected;
				}
				else
				{
					for(int i=0; i<local_boundary_intersected.size(); ++i)
					{
						if(local_boundary_intersected[i])
						{
							invariant_boundary_intersected[i] = true;
						}
					}
				}
			}

			switch(type)
			{
			case -1:	// invariant violated
				return true;
			case 0:		// domain is not contracted
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 1: 	// time interval is not contracted
				newFlowpipe.normalize();
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 2: 	// time interval is contracted
				if(contracted_domain[0] > intZero)
				{
					return true;
				}
				else
				{
					resultsCompo.push_back(tmvCompo);
					domains.push_back(contracted_domain);
					return true;
				}
			}

			t += step;

			if(bPrint)
			{
				printf("mode: %s,\t", modeNames[mode].c_str());
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("order = %d\n", newOrder);
			}

			if(newOrder > order)
			{
				if(newOrder > localMaxOrder)
				{
					for(int i=localMaxOrder+1; i<=newOrder; ++i)
					{
						vector<HornerForm> newTaylorExpansion;
						computeTaylorExpansion(newTaylorExpansion, polyODE, i);
						expansions.push_back(newTaylorExpansion);
					}

					localMaxOrder = newOrder;
				}

				--newOrder;
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			return false;
		}
	}

	return true;
}

bool HybridSystem::reach_continuous_low_degree(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
		const double step, const double time, const vector<int> & orders, const vector<int> & maxOrders, const int globalMaxOrder,
		const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
		vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;
	Interval intZero;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	resultsCompo.clear();

	TaylorModelVec tmvTemp;
	initFp.composition_normal(tmvTemp, step_exp_table);
	resultsCompo.push_back(tmvTemp);
	domains.push_back(initFp.domain);

	Flowpipe newFlowpipe, currentFlowpipe = initFp;

	vector<int> newOrders = orders;
	vector<int> localMaxOrders = orders;

	vector<Polynomial> polyODE;
	for(int i=0; i<odes[mode].tms.size(); ++i)
	{
		polyODE.push_back(odes[mode].tms[i].expansion);
	}

	vector<HornerForm> taylorExpansionHF;
	vector<Polynomial> taylorExpansionMF;
	vector<Polynomial> highestTerms;

	computeTaylorExpansion(taylorExpansionHF, taylorExpansionMF, highestTerms, polyODE, orders);

	vector<vector<HornerForm> > expansions;
	vector<HornerForm> emptySet;
	for(int i=0; i<taylorExpansionHF.size(); ++i)
	{
		expansions.push_back(emptySet);
		expansions[i].push_back(taylorExpansionHF[i]);
	}

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int localMaxOrder = newOrders[0];
		for(int i=1; i<newOrders.size(); ++i)
		{
			if(localMaxOrder < newOrders[i])
				localMaxOrder = newOrders[i];
		}

		bool bvalid = currentFlowpipe.advance_low_degree(newFlowpipe, hfOdes[mode], taylorExpansionHF, precondition, step_exp_table, step_end_exp_table, newOrders, maxOrders, estimation, uncertainties[mode]);

		if(bvalid)
		{
			// over-approximate the flowpipe/invariant intersection
			TaylorModelVec tmvCompo;
			newFlowpipe.composition_normal(tmvCompo, step_exp_table);

			vector<Interval> contracted_domain = newFlowpipe.domain;
			vector<bool> local_boundary_intersected;
			int type = contract_interval_arithmetic(tmvCompo, contracted_domain, invariants[mode], local_boundary_intersected);

			if(type >= 0)
			{
				// collect the intersected invariant boundary
				if(invariant_boundary_intersected.size() != local_boundary_intersected.size())
				{
					invariant_boundary_intersected = local_boundary_intersected;
				}
				else
				{
					for(int i=0; i<local_boundary_intersected.size(); ++i)
					{
						if(local_boundary_intersected[i])
						{
							invariant_boundary_intersected[i] = true;
						}
					}
				}
			}

			switch(type)
			{
			case -1:	// invariant violated
				return true;
			case 0:		// domain is not contracted
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 1: 	// time interval is not contracted
			{
				newFlowpipe.normalize();
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			}
			case 2: 	// time interval is contracted
				if(contracted_domain[0] > intZero)
				{
					return true;
				}
				else
				{
					resultsCompo.push_back(tmvCompo);
					domains.push_back(contracted_domain);
					return true;
				}
			}

			t += step;

			if(bPrint)
			{
				printf("mode: %s,\t", modeNames[mode].c_str());
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("orders:\t");
				int num = orders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), newOrders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), newOrders[num]);
			}

			for(int i=0; i<newOrders.size(); ++i)
			{
				if(newOrders[i] > orders[i])
				{
					--newOrders[i];

					if(newOrders[i] > localMaxOrders[i])
					{
						for(int j=localMaxOrders[i]; j<newOrders[i]; ++j)
						{
							HornerForm newTaylorExpansionHF;
							Polynomial newTaylorExpansionMF;

							increaseExpansionOrder(newTaylorExpansionHF, newTaylorExpansionMF, highestTerms[i], taylorExpansionMF[i], polyODE, j);

							expansions[i].push_back(newTaylorExpansionHF);
							taylorExpansionMF[i] = newTaylorExpansionMF;
						}
					}

					localMaxOrders[i] = newOrders[i];

					taylorExpansionHF[i] = expansions[i][newOrders[i]-orders[i]];
				}
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			return false;
		}
	}

	return true;
}

// for high-degree ODEs
// fixed step sizes and orders

bool HybridSystem::reach_continuous_high_degree(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp, const double step,
		const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
		vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;
	Interval intZero;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);

	resultsCompo.clear();

	TaylorModelVec tmvTemp;
	initFp.composition_normal(tmvTemp, step_exp_table);
	resultsCompo.push_back(tmvTemp);
	domains.push_back(initFp.domain);

	Flowpipe newFlowpipe, currentFlowpipe = initFp;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_high_degree(newFlowpipe, hfOdes[mode], precondition, step_exp_table, step_end_exp_table, order, estimation, uncertainties[mode]);

		if(bvalid)
		{
			// over-approximate the flowpipe/invariant intersection
			TaylorModelVec tmvCompo;
			newFlowpipe.composition_normal(tmvCompo, step_exp_table);

			vector<Interval> contracted_domain = newFlowpipe.domain;
			vector<bool> local_boundary_intersected;
			int type = contract_interval_arithmetic(tmvCompo, contracted_domain, invariants[mode], local_boundary_intersected);

			if(type >= 0)
			{
				// collect the intersected invariant boundary
				if(invariant_boundary_intersected.size() != local_boundary_intersected.size())
				{
					invariant_boundary_intersected = local_boundary_intersected;
				}
				else
				{
					for(int i=0; i<local_boundary_intersected.size(); ++i)
					{
						if(local_boundary_intersected[i])
						{
							invariant_boundary_intersected[i] = true;
						}
					}
				}
			}

			switch(type)
			{
			case -1:	// invariant violated
				return true;
			case 0:		// domain is not contracted
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 1: 	// time interval is not contracted
				newFlowpipe.normalize();
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 2: 	// time interval is contracted
				if(contracted_domain[0] > intZero)
				{
					return true;
				}
				else
				{
					resultsCompo.push_back(tmvCompo);
					domains.push_back(contracted_domain);
					return true;
				}
			}

			t += step;

			if(bPrint)
			{
				printf("mode: %s,\t", modeNames[mode].c_str());
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("order = %d\n", order);
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			return false;
		}
	}

	return true;
}

bool HybridSystem::reach_continuous_high_degree(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
		const double step, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition,
		const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
		vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;
	Interval intZero;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	resultsCompo.clear();

	TaylorModelVec tmvTemp;
	initFp.composition_normal(tmvTemp, step_exp_table);
	resultsCompo.push_back(tmvTemp);
	domains.push_back(initFp.domain);

	Flowpipe newFlowpipe, currentFlowpipe = initFp;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_high_degree(newFlowpipe, hfOdes[mode], precondition, step_exp_table, step_end_exp_table, orders, globalMaxOrder, estimation, uncertainties[mode]);

		if(bvalid)
		{
			// over-approximate the flowpipe/invariant intersection
			TaylorModelVec tmvCompo;
			newFlowpipe.composition_normal(tmvCompo, step_exp_table);

			vector<Interval> contracted_domain = newFlowpipe.domain;
			vector<bool> local_boundary_intersected;
			int type = contract_interval_arithmetic(tmvCompo, contracted_domain, invariants[mode], local_boundary_intersected);

			if(type >= 0)
			{
				// collect the intersected invariant boundary
				if(invariant_boundary_intersected.size() != local_boundary_intersected.size())
				{
					invariant_boundary_intersected = local_boundary_intersected;
				}
				else
				{
					for(int i=0; i<local_boundary_intersected.size(); ++i)
					{
						if(local_boundary_intersected[i])
						{
							invariant_boundary_intersected[i] = true;
						}
					}
				}
			}

			switch(type)
			{
			case -1:	// invariant violated
				return true;
			case 0:		// domain is not contracted
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 1: 	// time interval is not contracted
				newFlowpipe.normalize();
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 2: 	// time interval is contracted
				if(contracted_domain[0] > intZero)
				{
					return true;
				}
				else
				{
					resultsCompo.push_back(tmvCompo);
					domains.push_back(contracted_domain);
					return true;
				}
			}

			t += step;

			if(bPrint)
			{
				printf("mode: %s,\t", modeNames[mode].c_str());
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("orders:\t");
				int num = orders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), orders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), orders[num]);
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			return false;
		}
	}

	return true;
}

// adaptive step sizes and fixed orders

bool HybridSystem::reach_continuous_high_degree(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
		const double step, const double miniStep, const double time, const int order, const int precondition,
		const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
		vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;
	Interval intZero;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);

	resultsCompo.clear();

	TaylorModelVec tmvTemp;
	initFp.composition_normal(tmvTemp, step_exp_table);
	resultsCompo.push_back(tmvTemp);
	domains.push_back(initFp.domain);

	Flowpipe newFlowpipe, currentFlowpipe = initFp;

	double newStep = 0;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_high_degree(newFlowpipe, hfOdes[mode], precondition, step_exp_table, step_end_exp_table, newStep, miniStep, order, estimation, uncertainties[mode]);

		if(bvalid)
		{
			// over-approximate the flowpipe/invariant intersection
			TaylorModelVec tmvCompo;
			newFlowpipe.composition_normal(tmvCompo, step_exp_table);

			vector<Interval> contracted_domain = newFlowpipe.domain;
			vector<bool> local_boundary_intersected;
			int type = contract_interval_arithmetic(tmvCompo, contracted_domain, invariants[mode], local_boundary_intersected);

			if(type >= 0)
			{
				// collect the intersected invariant boundary
				if(invariant_boundary_intersected.size() != local_boundary_intersected.size())
				{
					invariant_boundary_intersected = local_boundary_intersected;
				}
				else
				{
					for(int i=0; i<local_boundary_intersected.size(); ++i)
					{
						if(local_boundary_intersected[i])
						{
							invariant_boundary_intersected[i] = true;
						}
					}
				}
			}

			switch(type)
			{
			case -1:	// invariant violated
				return true;
			case 0:		// domain is not contracted
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 1: 	// time interval is not contracted
				newFlowpipe.normalize();
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 2: 	// time interval is contracted
				if(contracted_domain[0] > intZero)
				{
					return true;
				}
				else
				{
					resultsCompo.push_back(tmvCompo);
					domains.push_back(contracted_domain);
					return true;
				}
			}

			t += step_exp_table[1].sup();

			if(bPrint)
			{
				printf("mode: %s,\t", modeNames[mode].c_str());
				printf("time = %f,\t", t);
				printf("step = %f,\t", step_exp_table[1].sup());
				printf("order = %d\n", order);
			}

			newStep = step_exp_table[1].sup() * LAMBDA_UP;
			if(newStep > step - THRESHOLD_HIGH)
			{
				newStep = 0;
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			return false;
		}
	}

	return true;
}

bool HybridSystem::reach_continuous_high_degree(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
		const double step, const double miniStep, const double time, const vector<int> & orders, const int globalMaxOrder,
		const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
		vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;
	Interval intZero;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	resultsCompo.clear();

	TaylorModelVec tmvTemp;
	initFp.composition_normal(tmvTemp, step_exp_table);
	resultsCompo.push_back(tmvTemp);
	domains.push_back(initFp.domain);

	Flowpipe newFlowpipe, currentFlowpipe = initFp;

	double newStep = 0;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_high_degree(newFlowpipe, hfOdes[mode], precondition, step_exp_table, step_end_exp_table, newStep, miniStep, orders, globalMaxOrder, estimation, uncertainties[mode]);

		if(bvalid)
		{
			// over-approximate the flowpipe/invariant intersection
			TaylorModelVec tmvCompo;
			newFlowpipe.composition_normal(tmvCompo, step_exp_table);

			vector<Interval> contracted_domain = newFlowpipe.domain;
			vector<bool> local_boundary_intersected;
			int type = contract_interval_arithmetic(tmvCompo, contracted_domain, invariants[mode], local_boundary_intersected);

			if(type >= 0)
			{
				// collect the intersected invariant boundary
				if(invariant_boundary_intersected.size() != local_boundary_intersected.size())
				{
					invariant_boundary_intersected = local_boundary_intersected;
				}
				else
				{
					for(int i=0; i<local_boundary_intersected.size(); ++i)
					{
						if(local_boundary_intersected[i])
						{
							invariant_boundary_intersected[i] = true;
						}
					}
				}
			}

			switch(type)
			{
			case -1:	// invariant violated
				return true;
			case 0:		// domain is not contracted
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 1: 	// time interval is not contracted
				newFlowpipe.normalize();
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 2: 	// time interval is contracted
				if(contracted_domain[0] > intZero)
				{
					return true;
				}
				else
				{
					resultsCompo.push_back(tmvCompo);
					domains.push_back(contracted_domain);
					return true;
				}
			}

			t += step_exp_table[1].sup();

			if(bPrint)
			{
				printf("mode: %s,\t", modeNames[mode].c_str());
				printf("time = %f,\t", t);
				printf("step = %f,\t", step_exp_table[1].sup());
				printf("orders:\t");
				int num = orders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), orders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), orders[num]);
			}

			newStep = step_exp_table[1].sup() * LAMBDA_UP;
			if(newStep > step - THRESHOLD_HIGH)
			{
				newStep = 0;
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			return false;
		}
	}

	return true;
}

// adaptive orders and fixed step sizes

bool HybridSystem::reach_continuous_high_degree(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
		const double step, const double time, const int order, const int maxOrder, const int precondition, const vector<Interval> & estimation,
		const bool bPrint, const vector<string> & stateVarNames, vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;
	Interval intZero;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*maxOrder);

	resultsCompo.clear();

	TaylorModelVec tmvTemp;
	initFp.composition_normal(tmvTemp, step_exp_table);
	resultsCompo.push_back(tmvTemp);
	domains.push_back(initFp.domain);

	Flowpipe newFlowpipe, currentFlowpipe = initFp;

	int newOrder = order;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_high_degree(newFlowpipe, hfOdes[mode], precondition, step_exp_table, step_end_exp_table, newOrder, maxOrder, estimation, uncertainties[mode]);

		if(bvalid)
		{
			// over-approximate the flowpipe/invariant intersection
			TaylorModelVec tmvCompo;
			newFlowpipe.composition_normal(tmvCompo, step_exp_table);

			vector<Interval> contracted_domain = newFlowpipe.domain;
			vector<bool> local_boundary_intersected;
			int type = contract_interval_arithmetic(tmvCompo, contracted_domain, invariants[mode], local_boundary_intersected);

			if(type >= 0)
			{
				// collect the intersected invariant boundary
				if(invariant_boundary_intersected.size() != local_boundary_intersected.size())
				{
					invariant_boundary_intersected = local_boundary_intersected;
				}
				else
				{
					for(int i=0; i<local_boundary_intersected.size(); ++i)
					{
						if(local_boundary_intersected[i])
						{
							invariant_boundary_intersected[i] = true;
						}
					}
				}
			}

			switch(type)
			{
			case -1:	// invariant violated
				return true;
			case 0:		// domain is not contracted
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 1: 	// time interval is not contracted
				newFlowpipe.normalize();
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 2: 	// time interval is contracted
				if(contracted_domain[0] > intZero)
				{
					return true;
				}
				else
				{
					resultsCompo.push_back(tmvCompo);
					domains.push_back(contracted_domain);
					return true;
				}
			}

			t += step;

			if(bPrint)
			{
				printf("mode: %s,\t", modeNames[mode].c_str());
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("order = %d\n", newOrder);
			}

			if(newOrder > order)
			{
				--newOrder;
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			return false;
		}
	}

	return true;
}

bool HybridSystem::reach_continuous_high_degree(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
		const double step, const double time, const vector<int> & orders, const vector<int> & maxOrders, const int globalMaxOrder,
		const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
		vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;
	Interval intZero;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	resultsCompo.clear();

	TaylorModelVec tmvTemp;
	initFp.composition_normal(tmvTemp, step_exp_table);
	resultsCompo.push_back(tmvTemp);
	domains.push_back(initFp.domain);

	Flowpipe newFlowpipe, currentFlowpipe = initFp;

	vector<int> newOrders = orders;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int localMaxOrder = newOrders[0];
		for(int i=1; i<newOrders.size(); ++i)
		{
			if(localMaxOrder < newOrders[i])
				localMaxOrder = newOrders[i];
		}

		bool bvalid = currentFlowpipe.advance_high_degree(newFlowpipe, hfOdes[mode], precondition, step_exp_table, step_end_exp_table, newOrders, localMaxOrder, maxOrders, estimation, uncertainties[mode]);

		if(bvalid)
		{
			// over-approximate the flowpipe/invariant intersection
			TaylorModelVec tmvCompo;
			newFlowpipe.composition_normal(tmvCompo, step_exp_table);

			vector<Interval> contracted_domain = newFlowpipe.domain;
			vector<bool> local_boundary_intersected;
			int type = contract_interval_arithmetic(tmvCompo, contracted_domain, invariants[mode], local_boundary_intersected);

			if(type >= 0)
			{
				// collect the intersected invariant boundary
				if(invariant_boundary_intersected.size() != local_boundary_intersected.size())
				{
					invariant_boundary_intersected = local_boundary_intersected;
				}
				else
				{
					for(int i=0; i<local_boundary_intersected.size(); ++i)
					{
						if(local_boundary_intersected[i])
						{
							invariant_boundary_intersected[i] = true;
						}
					}
				}
			}

			switch(type)
			{
			case -1:	// invariant violated
				return true;
			case 0:		// domain is not contracted
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 1: 	// time interval is not contracted
				newFlowpipe.normalize();
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 2: 	// time interval is contracted
				if(contracted_domain[0] > intZero)
				{
					return true;
				}
				else
				{
					resultsCompo.push_back(tmvCompo);
					domains.push_back(contracted_domain);
					return true;
				}
			}

			t += step;

			if(bPrint)
			{
				printf("mode: %s,\t", modeNames[mode].c_str());
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("orders:\t");
				int num = orders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), newOrders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), newOrders[num]);
			}

			for(int i=0; i<newOrders.size(); ++i)
			{
				if(newOrders[i] > orders[i])
					--newOrders[i];
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			return false;
		}
	}

	return true;
}


// for non-polynomial ODEs (using Taylor approximations)
// fixed step sizes and orders
bool HybridSystem::reach_continuous_non_polynomial_taylor(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
		const double step, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint,
		const vector<string> & stateVarNames, vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;
	Interval intZero;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);

	resultsCompo.clear();

	TaylorModelVec tmvTemp;
	initFp.composition_normal(tmvTemp, step_exp_table);
	resultsCompo.push_back(tmvTemp);
	domains.push_back(initFp.domain);

	Flowpipe newFlowpipe, currentFlowpipe = initFp;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_non_polynomial_taylor(newFlowpipe, strOdes[mode], precondition, step_exp_table, step_end_exp_table, order, estimation, uncertainties[mode], uncertainty_centers[mode]);

		if(bvalid)
		{
			// over-approximate the flowpipe/invariant intersection
			TaylorModelVec tmvCompo;
			newFlowpipe.composition_normal(tmvCompo, step_exp_table);

			vector<Interval> contracted_domain = newFlowpipe.domain;
			vector<bool> local_boundary_intersected;
			int type = contract_interval_arithmetic(tmvCompo, contracted_domain, invariants[mode], local_boundary_intersected);

			if(type >= 0)
			{
				// collect the intersected invariant boundary
				if(invariant_boundary_intersected.size() != local_boundary_intersected.size())
				{
					invariant_boundary_intersected = local_boundary_intersected;
				}
				else
				{
					for(int i=0; i<local_boundary_intersected.size(); ++i)
					{
						if(local_boundary_intersected[i])
						{
							invariant_boundary_intersected[i] = true;
						}
					}
				}
			}

			switch(type)
			{
			case -1:	// invariant violated
				return true;
			case 0:		// domain is not contracted
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 1: 	// time interval is not contracted
				newFlowpipe.normalize();
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 2: 	// time interval is contracted
				if(contracted_domain[0] > intZero)
				{
					return true;
				}
				else
				{
					resultsCompo.push_back(tmvCompo);
					domains.push_back(contracted_domain);
					return true;
				}
			}

			t += step;

			if(bPrint)
			{
				printf("mode: %s,\t", modeNames[mode].c_str());
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("order = %d\n", order);
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			return false;
		}
	}

	return true;
}

bool HybridSystem::reach_continuous_non_polynomial_taylor(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
		const double step, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition,
		const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
		vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;
	Interval intZero;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	resultsCompo.clear();

	TaylorModelVec tmvTemp;
	initFp.composition_normal(tmvTemp, step_exp_table);
	resultsCompo.push_back(tmvTemp);
	domains.push_back(initFp.domain);

	Flowpipe newFlowpipe, currentFlowpipe = initFp;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_non_polynomial_taylor(newFlowpipe, strOdes[mode], precondition, step_exp_table, step_end_exp_table, orders, globalMaxOrder, estimation, uncertainties[mode], uncertainty_centers[mode]);

		if(bvalid)
		{
			// over-approximate the flowpipe/invariant intersection
			TaylorModelVec tmvCompo;
			newFlowpipe.composition_normal(tmvCompo, step_exp_table);

			vector<Interval> contracted_domain = newFlowpipe.domain;
			vector<bool> local_boundary_intersected;
			int type = contract_interval_arithmetic(tmvCompo, contracted_domain, invariants[mode], local_boundary_intersected);

			if(type >= 0)
			{
				// collect the intersected invariant boundary
				if(invariant_boundary_intersected.size() != local_boundary_intersected.size())
				{
					invariant_boundary_intersected = local_boundary_intersected;
				}
				else
				{
					for(int i=0; i<local_boundary_intersected.size(); ++i)
					{
						if(local_boundary_intersected[i])
						{
							invariant_boundary_intersected[i] = true;
						}
					}
				}
			}

			switch(type)
			{
			case -1:	// invariant violated
				return true;
			case 0:		// domain is not contracted
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 1: 	// time interval is not contracted
				newFlowpipe.normalize();
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 2: 	// time interval is contracted
				if(contracted_domain[0] > intZero)
				{
					return true;
				}
				else
				{
					resultsCompo.push_back(tmvCompo);
					domains.push_back(contracted_domain);
					return true;
				}
			}

			t += step;

			if(bPrint)
			{
				printf("mode: %s,\t", modeNames[mode].c_str());
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("orders:\t");
				int num = orders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), orders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), orders[num]);
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			return false;
		}
	}

	return true;
}


// adaptive step sizes and fixed orders
bool HybridSystem::reach_continuous_non_polynomial_taylor(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
		const double step, const double miniStep, const double time, const int order, const int precondition,
		const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
		vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;
	Interval intZero;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);

	resultsCompo.clear();

	TaylorModelVec tmvTemp;
	initFp.composition_normal(tmvTemp, step_exp_table);
	resultsCompo.push_back(tmvTemp);
	domains.push_back(initFp.domain);

	Flowpipe newFlowpipe, currentFlowpipe = initFp;

	double newStep = 0;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_non_polynomial_taylor(newFlowpipe, strOdes[mode], precondition, step_exp_table, step_end_exp_table, newStep, miniStep, order, estimation, uncertainties[mode], uncertainty_centers[mode]);

		if(bvalid)
		{
			// over-approximate the flowpipe/invariant intersection
			TaylorModelVec tmvCompo;
			newFlowpipe.composition_normal(tmvCompo, step_exp_table);

			vector<Interval> contracted_domain = newFlowpipe.domain;
			vector<bool> local_boundary_intersected;
			int type = contract_interval_arithmetic(tmvCompo, contracted_domain, invariants[mode], local_boundary_intersected);

			if(type >= 0)
			{
				// collect the intersected invariant boundary
				if(invariant_boundary_intersected.size() != local_boundary_intersected.size())
				{
					invariant_boundary_intersected = local_boundary_intersected;
				}
				else
				{
					for(int i=0; i<local_boundary_intersected.size(); ++i)
					{
						if(local_boundary_intersected[i])
						{
							invariant_boundary_intersected[i] = true;
						}
					}
				}
			}

			switch(type)
			{
			case -1:	// invariant violated
				return true;
			case 0:		// domain is not contracted
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 1: 	// time interval is not contracted
				newFlowpipe.normalize();
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 2: 	// time interval is contracted
				if(contracted_domain[0] > intZero)
				{
					return true;
				}
				else
				{
					resultsCompo.push_back(tmvCompo);
					domains.push_back(contracted_domain);
					return true;
				}
			}

			t += step_exp_table[1].sup();

			if(bPrint)
			{
				printf("mode: %s,\t", modeNames[mode].c_str());
				printf("time = %f,\t", t);
				printf("step = %f,\t", step_exp_table[1].sup());
				printf("order = %d\n", order);
			}

			newStep = step_exp_table[1].sup() * LAMBDA_UP;
			if(newStep > step - THRESHOLD_HIGH)
			{
				newStep = 0;
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			return false;
		}
	}

	return true;
}

bool HybridSystem::reach_continuous_non_polynomial_taylor(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
		const double step, const double miniStep, const double time, const vector<int> & orders, const int globalMaxOrder,
		const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
		vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;
	Interval intZero;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	resultsCompo.clear();

	TaylorModelVec tmvTemp;
	initFp.composition_normal(tmvTemp, step_exp_table);
	resultsCompo.push_back(tmvTemp);
	domains.push_back(initFp.domain);

	Flowpipe newFlowpipe, currentFlowpipe = initFp;

	double newStep = 0;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_non_polynomial_taylor(newFlowpipe, strOdes[mode], precondition, step_exp_table, step_end_exp_table, newStep, miniStep, orders, globalMaxOrder, estimation, uncertainties[mode], uncertainty_centers[mode]);

		if(bvalid)
		{
			// over-approximate the flowpipe/invariant intersection
			TaylorModelVec tmvCompo;
			newFlowpipe.composition_normal(tmvCompo, step_exp_table);

			vector<Interval> contracted_domain = newFlowpipe.domain;
			vector<bool> local_boundary_intersected;
			int type = contract_interval_arithmetic(tmvCompo, contracted_domain, invariants[mode], local_boundary_intersected);

			if(type >= 0)
			{
				// collect the intersected invariant boundary
				if(invariant_boundary_intersected.size() != local_boundary_intersected.size())
				{
					invariant_boundary_intersected = local_boundary_intersected;
				}
				else
				{
					for(int i=0; i<local_boundary_intersected.size(); ++i)
					{
						if(local_boundary_intersected[i])
						{
							invariant_boundary_intersected[i] = true;
						}
					}
				}
			}

			switch(type)
			{
			case -1:	// invariant violated
				return true;
			case 0:		// domain is not contracted
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 1: 	// time interval is not contracted
				newFlowpipe.normalize();
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 2: 	// time interval is contracted
				if(contracted_domain[0] > intZero)
				{
					return true;
				}
				else
				{
					resultsCompo.push_back(tmvCompo);
					domains.push_back(contracted_domain);
					return true;
				}
			}

			t += step_exp_table[1].sup();

			if(bPrint)
			{
				printf("mode: %s,\t", modeNames[mode].c_str());
				printf("time = %f,\t", t);
				printf("step = %f,\t", step_exp_table[1].sup());
				printf("orders:\t");
				int num = orders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), orders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), orders[num]);
			}

			newStep = step_exp_table[1].sup() * LAMBDA_UP;
			if(newStep > step - THRESHOLD_HIGH)
			{
				newStep = 0;
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			return false;
		}
	}

	return true;
}


// adaptive orders and fixed step sizes
bool HybridSystem::reach_continuous_non_polynomial_taylor(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
		const double step, const double time, const int order, const int maxOrder, const int precondition, const vector<Interval> & estimation,
		const bool bPrint, const vector<string> & stateVarNames, vector<bool> & invariant_boundary_intersected,
		const vector<string> & modeNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;
	Interval intZero;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*maxOrder);

	resultsCompo.clear();

	TaylorModelVec tmvTemp;
	initFp.composition_normal(tmvTemp, step_exp_table);
	resultsCompo.push_back(tmvTemp);
	domains.push_back(initFp.domain);

	Flowpipe newFlowpipe, currentFlowpipe = initFp;

	int newOrder = order;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		bool bvalid = currentFlowpipe.advance_non_polynomial_taylor(newFlowpipe, strOdes[mode], precondition, step_exp_table, step_end_exp_table, newOrder, maxOrder, estimation, uncertainties[mode], uncertainty_centers[mode]);

		if(bvalid)
		{
			// over-approximate the flowpipe/invariant intersection
			TaylorModelVec tmvCompo;
			newFlowpipe.composition_normal(tmvCompo, step_exp_table);

			vector<Interval> contracted_domain = newFlowpipe.domain;
			vector<bool> local_boundary_intersected;
			int type = contract_interval_arithmetic(tmvCompo, contracted_domain, invariants[mode], local_boundary_intersected);

			if(type >= 0)
			{
				// collect the intersected invariant boundary
				if(invariant_boundary_intersected.size() != local_boundary_intersected.size())
				{
					invariant_boundary_intersected = local_boundary_intersected;
				}
				else
				{
					for(int i=0; i<local_boundary_intersected.size(); ++i)
					{
						if(local_boundary_intersected[i])
						{
							invariant_boundary_intersected[i] = true;
						}
					}
				}
			}

			switch(type)
			{
			case -1:	// invariant violated
				return true;
			case 0:		// domain is not contracted
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 1: 	// time interval is not contracted
				newFlowpipe.normalize();
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 2: 	// time interval is contracted
				if(contracted_domain[0] > intZero)
				{
					return true;
				}
				else
				{
					resultsCompo.push_back(tmvCompo);
					domains.push_back(contracted_domain);
					return true;
				}
			}

			t += step;

			if(bPrint)
			{
				printf("mode: %s,\t", modeNames[mode].c_str());
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("order = %d\n", newOrder);
			}

			if(newOrder > order)
			{
				--newOrder;
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			return false;
		}
	}

	return true;
}

bool HybridSystem::reach_continuous_non_polynomial_taylor(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
		const double step, const double time, const vector<int> & orders, const vector<int> & maxOrders, const int globalMaxOrder,
		const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
		vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames) const
{
	vector<Interval> step_exp_table, step_end_exp_table;
	Interval intZero;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	resultsCompo.clear();

	TaylorModelVec tmvTemp;
	initFp.composition_normal(tmvTemp, step_exp_table);
	resultsCompo.push_back(tmvTemp);
	domains.push_back(initFp.domain);

	Flowpipe newFlowpipe, currentFlowpipe = initFp;

	vector<int> newOrders = orders;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int localMaxOrder = newOrders[0];
		for(int i=1; i<newOrders.size(); ++i)
		{
			if(localMaxOrder < newOrders[i])
				localMaxOrder = newOrders[i];
		}

		bool bvalid = currentFlowpipe.advance_non_polynomial_taylor(newFlowpipe, strOdes[mode], precondition, step_exp_table, step_end_exp_table, newOrders, localMaxOrder, maxOrders, estimation, uncertainties[mode], uncertainty_centers[mode]);

		if(bvalid)
		{
			// over-approximate the flowpipe/invariant intersection
			TaylorModelVec tmvCompo;
			newFlowpipe.composition_normal(tmvCompo, step_exp_table);

			vector<Interval> contracted_domain = newFlowpipe.domain;
			vector<bool> local_boundary_intersected;
			int type = contract_interval_arithmetic(tmvCompo, contracted_domain, invariants[mode], local_boundary_intersected);

			if(type >= 0)
			{
				// collect the intersected invariant boundary
				if(invariant_boundary_intersected.size() != local_boundary_intersected.size())
				{
					invariant_boundary_intersected = local_boundary_intersected;
				}
				else
				{
					for(int i=0; i<local_boundary_intersected.size(); ++i)
					{
						if(local_boundary_intersected[i])
						{
							invariant_boundary_intersected[i] = true;
						}
					}
				}
			}

			switch(type)
			{
			case -1:	// invariant violated
				return true;
			case 0:		// domain is not contracted
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 1: 	// time interval is not contracted
				newFlowpipe.normalize();
				currentFlowpipe = newFlowpipe;
				resultsCompo.push_back(tmvCompo);
				domains.push_back(contracted_domain);
				break;
			case 2: 	// time interval is contracted
				if(contracted_domain[0] > intZero)
				{
					return true;
				}
				else
				{
					resultsCompo.push_back(tmvCompo);
					domains.push_back(contracted_domain);
					return true;
				}
			}

			t += step;

			if(bPrint)
			{
				printf("mode: %s,\t", modeNames[mode].c_str());
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("orders:\t");
				int num = orders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), newOrders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), newOrders[num]);
			}

			for(int i=0; i<newOrders.size(); ++i)
			{
				if(newOrders[i] > orders[i])
					--newOrders[i];
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			return false;
		}
	}

	return true;
}


// hybrid reachability

void HybridSystem::reach_hybrid(list<list<TaylorModelVec> > & resultsCompo, list<list<vector<Interval> > > & domains, list<int> & modeIDs, list<TreeNode *> & traceNodes,
		TreeNode * & traceTree, const vector<int> & integrationSchemes, const double step, const double miniStep,
		const double time, const int orderType, const vector<int> & orders, const vector<int> & maxOrders, const int globalMaxOrder,
		const bool bAdaptiveSteps, const bool bAdaptiveOrders, const int maxJmps, const int precondition, const vector<Interval> & estimation,
		const vector<vector<int> > & aggregType, const vector<vector<vector<RowVector> > > aggregationTemplate_candidates, const vector<RowVector> default_aggregation_template,
		const vector<vector<Matrix> > & weightTab, const vector<vector<vector<bool> > > & linear_auto, const vector<vector<vector<RowVector> > > & template_auto,
		const bool bPrint, const vector<string> & stateVarNames, const vector<string> & modeNames, const vector<string> & tmVarNames) const
{
	list<int> modeQueue;
	list<Flowpipe> flowpipeQueue;
	list<double> timePassedQueue;
	list<int> jumpsExecutedQueue;

	Interval intZero;
	int rangeDim = initialSet.tmv.tms.size();

	modeQueue.push_back(initialMode);
	flowpipeQueue.push_back(initialSet);
	timePassedQueue.push_back(0);
	jumpsExecutedQueue.push_back(0);

	// mode trace
	list<TreeNode *> nodeQueue;
	traceTree = new TreeNode(0, initialMode, intZero);
	nodeQueue.push_back(traceTree);

	for(; modeQueue.size() != 0;)
	{
		int initMode = modeQueue.front();
		Flowpipe initFp = flowpipeQueue.front();
		double timePassed = timePassedQueue.front();
		int jumpsExecuted = jumpsExecutedQueue.front();
		TreeNode *node = nodeQueue.front();

		modeQueue.pop_front();
		flowpipeQueue.pop_front();
		timePassedQueue.pop_front();
		jumpsExecutedQueue.pop_front();
		nodeQueue.pop_front();

		list<TaylorModelVec> mode_flowpipes;
		list<vector<Interval> > mode_domains;

		bool bvalid;
		vector<bool> invariant_boundary_intersected;

		switch(integrationSchemes[initMode])
		{
		case LOW_DEGREE:
		{
			switch(orderType)
			{
			case UNIFORM:
				if(bAdaptiveSteps)
				{
					bvalid = reach_continuous_low_degree(mode_flowpipes, mode_domains, initMode, initFp, step, miniStep, time-timePassed, orders[0], precondition, estimation, bPrint, stateVarNames, invariant_boundary_intersected, modeNames);
				}
				else if(bAdaptiveOrders)
				{
					bvalid = reach_continuous_low_degree(mode_flowpipes, mode_domains, initMode, initFp, step, time-timePassed, orders[0], maxOrders[0], precondition, estimation, bPrint, stateVarNames, invariant_boundary_intersected, modeNames);
				}
				else
				{
					bvalid = reach_continuous_low_degree(mode_flowpipes, mode_domains, initMode, initFp, step, time-timePassed, orders[0], precondition, estimation, bPrint, stateVarNames, invariant_boundary_intersected, modeNames);
				}
				break;
			case MULTI:
				if(bAdaptiveSteps)
				{
					bvalid = reach_continuous_low_degree(mode_flowpipes, mode_domains, initMode, initFp, step, miniStep, time-timePassed, orders, globalMaxOrder, precondition, estimation, bPrint, stateVarNames, invariant_boundary_intersected, modeNames);
				}
				else if(bAdaptiveOrders)
				{
					bvalid = reach_continuous_low_degree(mode_flowpipes, mode_domains, initMode, initFp, step, time-timePassed, orders, maxOrders, globalMaxOrder, precondition, estimation, bPrint, stateVarNames, invariant_boundary_intersected, modeNames);
				}
				else
				{
					bvalid = reach_continuous_low_degree(mode_flowpipes, mode_domains, initMode, initFp, step, time-timePassed, orders, globalMaxOrder, precondition, estimation, bPrint, stateVarNames, invariant_boundary_intersected, modeNames);
				}
				break;
			}
			break;
		}

		case HIGH_DEGREE:
		{
			switch(orderType)
			{
			case UNIFORM:
				if(bAdaptiveSteps)
				{
					bvalid = reach_continuous_high_degree(mode_flowpipes, mode_domains, initMode, initFp, step, miniStep, time-timePassed, orders[0], precondition, estimation, bPrint, stateVarNames, invariant_boundary_intersected, modeNames);
				}
				else if(bAdaptiveOrders)
				{
					bvalid = reach_continuous_high_degree(mode_flowpipes, mode_domains, initMode, initFp, step, time-timePassed, orders[0], maxOrders[0], precondition, estimation, bPrint, stateVarNames, invariant_boundary_intersected, modeNames);
				}
				else
				{
					bvalid = reach_continuous_high_degree(mode_flowpipes, mode_domains, initMode, initFp, step, time-timePassed, orders[0], precondition, estimation, bPrint, stateVarNames, invariant_boundary_intersected, modeNames);
				}
				break;
			case MULTI:
				if(bAdaptiveSteps)
				{
					bvalid = reach_continuous_high_degree(mode_flowpipes, mode_domains, initMode, initFp, step, miniStep, time-timePassed, orders, globalMaxOrder, precondition, estimation, bPrint, stateVarNames, invariant_boundary_intersected, modeNames);
				}
				else if(bAdaptiveOrders)
				{
					bvalid = reach_continuous_high_degree(mode_flowpipes, mode_domains, initMode, initFp, step, time-timePassed, orders, maxOrders, globalMaxOrder, precondition, estimation, bPrint, stateVarNames, invariant_boundary_intersected, modeNames);
				}
				else
				{
					bvalid = reach_continuous_high_degree(mode_flowpipes, mode_domains, initMode, initFp, step, time-timePassed, orders, globalMaxOrder, precondition, estimation, bPrint, stateVarNames, invariant_boundary_intersected, modeNames);
				}
				break;
			}
			break;
		}

		case NONPOLY_TAYLOR:
		{
			switch(orderType)
			{
			case UNIFORM:
				if(bAdaptiveSteps)
				{
					bvalid = reach_continuous_non_polynomial_taylor(mode_flowpipes, mode_domains, initMode, initFp, step, miniStep, time-timePassed, orders[0], precondition, estimation, bPrint, stateVarNames, invariant_boundary_intersected, modeNames);
				}
				else if(bAdaptiveOrders)
				{
					bvalid = reach_continuous_non_polynomial_taylor(mode_flowpipes, mode_domains, initMode, initFp, step, time-timePassed, orders[0], maxOrders[0], precondition, estimation, bPrint, stateVarNames, invariant_boundary_intersected, modeNames);
				}
				else
				{
					bvalid = reach_continuous_non_polynomial_taylor(mode_flowpipes, mode_domains, initMode, initFp, step, time-timePassed, orders[0], precondition, estimation, bPrint, stateVarNames, invariant_boundary_intersected, modeNames);
				}
				break;
			case MULTI:
				if(bAdaptiveSteps)
				{
					bvalid = reach_continuous_non_polynomial_taylor(mode_flowpipes, mode_domains, initMode, initFp, step, miniStep, time-timePassed, orders, globalMaxOrder, precondition, estimation, bPrint, stateVarNames, invariant_boundary_intersected, modeNames);
				}
				else if(bAdaptiveOrders)
				{
					bvalid = reach_continuous_non_polynomial_taylor(mode_flowpipes, mode_domains, initMode, initFp, step, time-timePassed, orders, maxOrders, globalMaxOrder, precondition, estimation, bPrint, stateVarNames, invariant_boundary_intersected, modeNames);
				}
				else
				{
					bvalid = reach_continuous_non_polynomial_taylor(mode_flowpipes, mode_domains, initMode, initFp, step, time-timePassed, orders, globalMaxOrder, precondition, estimation, bPrint, stateVarNames, invariant_boundary_intersected, modeNames);
				}
				break;
			}
			break;
		}
		}

		list<TaylorModelVec>::iterator tmvIter = mode_flowpipes.begin();
		list<vector<Interval> >::iterator doIter = mode_domains.begin();

		resultsCompo.push_back(mode_flowpipes);
		domains.push_back(mode_domains);
		modeIDs.push_back(initMode);
		traceNodes.push_back(node);

		if(!bvalid)
		{
			return;
		}

		if(jumpsExecuted == maxJmps)
		{
			if(bPrint)
			{
				printf("Maximum jump depth is reached.\n");
			}
			continue;
		}

		vector<vector<Interval> > step_exp_tables;
		vector<Interval> step_exp_table;

		for(; tmvIter!=mode_flowpipes.end(); ++tmvIter, ++doIter)
		{
			if(step_exp_table.size() == 0 || step_exp_table[1] != (*doIter)[0])
			{
				construct_step_exp_table(step_exp_table, (*doIter)[0], globalMaxOrder);
			}

			step_exp_tables.push_back(step_exp_table);
		}

		// over-approximate the intersection for each jump
		for(int i=0; i<transitions[initMode].size(); ++i)
		{
			vector<TaylorModelVec> intersection_flowpipes;
			vector<vector<Interval> > intersection_domains;

			list<TaylorModelVec>::iterator tmvIter = mode_flowpipes.begin();
			list<vector<Interval> >::iterator doIter = mode_domains.begin();
			double newTimePassed = 0;
			bool brecorded = false;

			if(bPrint)
			{
				printf("Dealing with the jump from %s to %s ...\n", modeNames[initMode].c_str(), modeNames[transitions[initMode][i].targetID].c_str());
			}

			// collect the intersected flowpipes

			vector<bool> guard_boundary_intersected;

			step_exp_table = step_exp_tables[0];

			Interval triggeredTime;

			for(int k=0; tmvIter!=mode_flowpipes.end(); ++tmvIter, ++doIter)
			{
				TaylorModelVec tmvIntersection = *tmvIter;
				vector<Interval> doIntersection = *doIter;

				if(step_exp_table[1] != (*doIter)[0])
				{
					step_exp_table = step_exp_tables[++k];
				}

				vector<bool> local_boundary_intersected;
				int type = contract_interval_arithmetic(tmvIntersection, doIntersection, transitions[initMode][i].guard, local_boundary_intersected);

				if(type >= 0 && aggregType[initMode][transitions[initMode][i].targetID] == PARA_AGGREG)
				{
					// collect the intersected guard boundary
					if(guard_boundary_intersected.size() != local_boundary_intersected.size())
					{
						guard_boundary_intersected = local_boundary_intersected;
					}
					else
					{
						for(int j=0; j<local_boundary_intersected.size(); ++j)
						{
							if(local_boundary_intersected[j])
							{
								guard_boundary_intersected[j] = true;
							}
						}
					}
				}

				if(type != -1)
				{
					if(!brecorded)
					{
						brecorded = true;

						triggeredTime = doIntersection[0];
						triggeredTime.add_assign(newTimePassed);
						newTimePassed += doIntersection[0].inf();
					}
					else
					{
						// compute the time interval when the jump is triggered
						triggeredTime.setSup(triggeredTime.sup() + doIntersection[0].sup());
					}

					intersection_flowpipes.push_back(tmvIntersection);
					intersection_domains.push_back(doIntersection);
				}
				else
				{
					if(!brecorded)
					{
						newTimePassed += doIntersection[0].sup();
					}
				}
			}

			vector<bool> boundary_intersected = invariant_boundary_intersected;
			for(int j=0; j<guard_boundary_intersected.size(); ++j)
			{
				boundary_intersected.push_back(guard_boundary_intersected[j]);
			}

			// aggregate the intersections
			if(intersection_flowpipes.size() > 0)
			{
				Flowpipe fpAggregation;
				TaylorModelVec tmvAggregation;
				vector<Interval> doAggregation;

				switch(aggregType[initMode][transitions[initMode][i].targetID])
				{
				case INTERVAL_AGGREG:
				{
					aggregate_flowpipes_by_interval(tmvAggregation, doAggregation, intersection_flowpipes, intersection_domains);

					break;
				}
				case PARA_AGGREG:
				{
					aggregate_flowpipes_by_Parallelotope(tmvAggregation, doAggregation, intersection_flowpipes, intersection_domains,
							invariants[initMode], transitions[initMode][i], boundary_intersected,
							aggregationTemplate_candidates[initMode][transitions[initMode][i].targetID], default_aggregation_template,
							weightTab, linear_auto, template_auto, globalMaxOrder, rangeDim);

					break;
				}
				}

				// contract the aggregation regarding to the guard and invariant
				vector<PolynomialConstraint> constraints = transitions[initMode][i].guard;
				for(int j=0; j<invariants[initMode].size(); ++j)
				{
					constraints.push_back(invariants[initMode][j]);
				}

				vector<bool> bVecDummy;
				contract_interval_arithmetic(tmvAggregation, doAggregation, constraints, bVecDummy);

				//reset map
				TaylorModelVec tmvImage;
				transitions[initMode][i].resetMap.reset(tmvImage, tmvAggregation, doAggregation, globalMaxOrder);

				int type = contract_interval_arithmetic(tmvImage, doAggregation, invariants[transitions[initMode][i].targetID], bVecDummy);

				if(type == -1)
				{
					if(bPrint)
					{
						printf("No intersection detected.\n");
					}

					continue;
				}
				else
				{
					int rangeDim = tmvImage.tms.size();
					Matrix coefficients(rangeDim, rangeDim+1);
					for(int i=0; i<rangeDim; ++i)
					{
						coefficients.set(1, i, i+1);
					}
					TaylorModelVec tmvTemp(coefficients);

					fpAggregation.tmv = tmvTemp;
					fpAggregation.tmvPre = tmvImage;
					fpAggregation.domain = doAggregation;
				}

				timePassed += newTimePassed;
				if(timePassed < time - THRESHOLD_HIGH)
				{
					modeQueue.push_back(transitions[initMode][i].targetID);
					flowpipeQueue.push_back(fpAggregation);
					timePassedQueue.push_back(timePassed);
					jumpsExecutedQueue.push_back(jumpsExecuted+1);

					TreeNode *child = new TreeNode(transitions[initMode][i].jumpID, transitions[initMode][i].targetID, triggeredTime);
					child->parent = node;
					node->children.push_back(child);

					nodeQueue.push_back(child);
				}
			}

			if(bPrint)
			{
				printf("Done.\n");
			}
		}

		step_exp_tables.clear();
	}
}























// class HybridReachability

HybridReachability::HybridReachability()
{
	traceTree = NULL;
	numOfJumps = 0;
}

HybridReachability::~HybridReachability()
{
	outputAxes.clear();
	orders.clear();
	maxOrders.clear();
	aggregationType.clear();
	default_aggregation_template.clear();
	aggregationTemplate_candidates.clear();
	aggregationTemplate_TaylorModel.clear();
	weightTab.clear();
	flowpipesCompo.clear();
	domains.clear();
	modeIDs.clear();
	stateVarTab.clear();
	stateVarNames.clear();
	tmVarTab.clear();
	tmVarNames.clear();
	modeTab.clear();
	modeNames.clear();
	unsafeSet.clear();
	traceNodes.clear();
	integrationSchemes.clear();

	delete traceTree;
}

void HybridReachability::dump(FILE *fp) const
{
	fprintf(fp,"state var ");
	for(int i=0; i<stateVarNames.size()-1; ++i)
	{
		fprintf(fp, "%s,", stateVarNames[i].c_str());
	}
	fprintf(fp, "%s\n\n", stateVarNames[stateVarNames.size()-1].c_str());

	// dump the constraints for the state space
	for(int i=0; i<system.invariants.size(); ++i)
	{
		fprintf(fp, "%s\n{\n", modeNames[i].c_str());

		for(int j=0; j<system.invariants[i].size(); ++j)
		{
			system.invariants[i][j].dump(fp, stateVarNames);
		}

		fprintf(fp, "}\n\n");
	}

	// dump the computation tree
	fprintf(fp, "computation paths\n{\n\n");

	string strEmpty;
	traceTree->dump(fp, strEmpty, modeNames);

	fprintf(fp, "}\n\n");

	switch(plotFormat)
	{
	case PLOT_GNUPLOT:
		switch(plotSetting)
		{
		case PLOT_INTERVAL:
			fprintf(fp, "gnuplot interval %s , %s\n\n", stateVarNames[outputAxes[0]].c_str(), stateVarNames[outputAxes[1]].c_str());
			break;
		case PLOT_OCTAGON:
			fprintf(fp, "gnuplot octagon %s , %s\n\n", stateVarNames[outputAxes[0]].c_str(), stateVarNames[outputAxes[1]].c_str());
			break;
		case PLOT_GRID:
			fprintf(fp, "gnuplot grid %d %s , %s\n\n", numSections, stateVarNames[outputAxes[0]].c_str(), stateVarNames[outputAxes[1]].c_str());
			break;
		}
		break;
	case PLOT_MATLAB:
		switch(plotSetting)
		{
		case PLOT_INTERVAL:
			fprintf(fp, "matlab interval %s , %s\n\n", stateVarNames[outputAxes[0]].c_str(), stateVarNames[outputAxes[1]].c_str());
			break;
		case PLOT_OCTAGON:
			fprintf(fp, "matlab octagon %s , %s\n\n", stateVarNames[outputAxes[0]].c_str(), stateVarNames[outputAxes[1]].c_str());
			break;
		case PLOT_GRID:
			fprintf(fp, "matlab grid %d %s , %s\n\n", numSections, stateVarNames[outputAxes[0]].c_str(), stateVarNames[outputAxes[1]].c_str());
			break;
		}
		break;
	}

	fprintf(fp, "output %s\n\n", outputFileName);

	if(bSafetyChecking)
	{
		// dump the unsafe set
		fprintf(fp, "unsafe set\n{\n");

		for(int i=0; i<bVecUnderCheck.size(); ++i)
		{
			if(bVecUnderCheck[i])
			{
				fprintf(fp, "%s\n{\n", modeNames[i].c_str());

				for(int j=0; j<unsafeSet[i].size(); ++j)
				{
					unsafeSet[i][j].dump(fp, stateVarNames);
				}

				fprintf(fp, "}\n\n");
			}
		}

		fprintf(fp, "}\n\n");
	}

	int rangeDim = system.hfOdes.size();

	list<list<TaylorModelVec> >::const_iterator fpIter = flowpipesCompo.begin();
	list<list<vector<Interval> > >::const_iterator fpdoIter = domains.begin();
	list<int>::const_iterator modeIter = modeIDs.begin();

	list<TaylorModelVec>::const_iterator tmvIter;
	list<vector<Interval> >::const_iterator doIter;

	fprintf(fp, "hybrid flowpipes\n{\n");

	for(; fpIter!=flowpipesCompo.end(); ++fpIter, ++fpdoIter, ++modeIter)
	{
		tmvIter = fpIter->begin();
		doIter = fpdoIter->begin();
		vector<string> newNames;

		if(doIter->size() != tmVarNames.size())
		{
			string tVar("local_t");
			newNames.push_back(tVar);

			// rename the variables
			for(int i=1; i<doIter->size(); ++i)
			{
				char name[NAME_SIZE];
				sprintf(name, "%s%d", local_var_name, i);
				string strName(name);
				newNames.push_back(strName);
			}
		}
		else
		{
			// do not change the name
			newNames = tmVarNames;
		}

		fprintf(fp, "%s\n{\n", modeNames[*modeIter].c_str());

		fprintf(fp, "tm var ");
		for(int i=1; i<newNames.size()-1; ++i)
		{
			fprintf(fp, "%s,", newNames[i].c_str());
		}
		fprintf(fp, "%s\n\n", newNames[newNames.size()-1].c_str());

		for(; tmvIter!=fpIter->end(); ++tmvIter, ++doIter)
		{
			fprintf(fp, "{\n");
			tmvIter->dump_interval(fp, stateVarNames, newNames);

			for(int i=0; i<doIter->size(); ++i)
			{
				fprintf(fp, "%s in ", newNames[i].c_str());
				(*doIter)[i].dump(fp);
				fprintf(fp, "\n");
			}

			fprintf(fp, "}\n\n");
		}

		fprintf(fp, "}\n\n");
	}

	fprintf(fp, "}\n");
}

void HybridReachability::run()
{
	// normalize the candidate vectors
	for(int i=0; i<aggregationTemplate_candidates.size(); ++i)
	{
		for(int j=0; j<aggregationTemplate_candidates[i].size(); ++j)
		{
			for(int k=0; k<aggregationTemplate_candidates[i][j].size(); ++k)
			{
				aggregationTemplate_candidates[i][j][k].normalize();
			}
		}
	}

	set_default_template();

	constructWeightTab();


/*
	// ==== test begin====
	// print out the jumps
	int num = 0;
	for(int i=0; i<system.transitions.size(); ++i)
	{
		for(int j=0; j<system.transitions[i].size(); ++j)
		{
			int start = system.transitions[i][j].startID;
			int end = system.transitions[i][j].targetID;

			printf("start:  %s,\tend:  %s\n", modeNames[start].c_str(), modeNames[end].c_str());

			switch(aggregationType[start][end])
			{
			case INTERVAL_AGGREG:
				printf("interval aggregation\n");
				break;
			case PARA_AGGREG:
				printf("paralleletope aggregation {\n");

				for(int k=0; k<aggregationTemplate_candidates[start][end].size(); ++k)
				{
					aggregationTemplate_candidates[start][end][k].dump(stdout);
				}

				printf("}\n");
				break;
			}

			++num;
		}
	}

	printf("total:  %d  jump(s)\n", num);

	exit(0);

	// ==== test end ====
*/

	compute_factorial_rec(globalMaxOrder+1);
	compute_power_4(globalMaxOrder+1);
	compute_double_factorial(2*globalMaxOrder);


	system.reach_hybrid(flowpipesCompo, domains, modeIDs, traceNodes, traceTree, integrationSchemes, step, miniStep, time, orderType, orders, maxOrders, globalMaxOrder, bAdaptiveSteps, bAdaptiveOrders,
			maxJumps, precondition, estimation, aggregationType, aggregationTemplate_candidates, default_aggregation_template, weightTab,
			linear_auto, template_auto, bPrint, stateVarNames, modeNames, tmVarNames);
}

void HybridReachability::plot_2D() const
{
	char filename[NAME_SIZE+10];

	switch(plotFormat)
	{
	case PLOT_GNUPLOT:
		sprintf(filename, "%s%s.plt", outputDir, outputFileName);
		break;
	case PLOT_MATLAB:
		sprintf(filename, "%s%s.m", outputDir, outputFileName);
		break;
	}

	FILE *fpPlotting = fopen(filename, "w");

	if(fpPlotting == NULL)
	{
		printf("Can not create the plotting file.\n");
		exit(1);
	}

	printf("Generating the plotting file...\n");
	switch(plotFormat)
	{
	case PLOT_GNUPLOT:
		plot_2D_GNUPLOT(fpPlotting);
		break;
	case PLOT_MATLAB:
		plot_2D_MATLAB(fpPlotting);
		break;
	}
	printf("Done.\n");

	fclose(fpPlotting);
}

void HybridReachability::plot_2D_GNUPLOT(FILE *fp) const
{
	switch(plotSetting)
	{
	case PLOT_INTERVAL:
		plot_2D_interval_GNUPLOT(fp);
		break;
	case PLOT_OCTAGON:
		plot_2D_octagon_GNUPLOT(fp);
		break;
	case PLOT_GRID:
		plot_2D_grid_GNUPLOT(fp);
		break;
	}
}

void HybridReachability::plot_2D_interval_GNUPLOT(FILE *fp) const
{
	fprintf(fp, "set terminal postscript\n");

	char filename[NAME_SIZE+10];
	sprintf(filename, "%s%s.eps", imageDir, outputFileName);
	fprintf(fp, "set output '%s'\n", filename);

	fprintf(fp, "set style line 1 linecolor rgb \"blue\"\n");
	fprintf(fp, "set autoscale\n");
	fprintf(fp, "unset label\n");
	fprintf(fp, "set xtic auto\n");
	fprintf(fp, "set ytic auto\n");
	fprintf(fp, "set xlabel \"%s\"\n", stateVarNames[outputAxes[0]].c_str());
	fprintf(fp, "set ylabel \"%s\"\n", stateVarNames[outputAxes[1]].c_str());
	fprintf(fp, "plot '-' notitle with lines ls 1\n");

	list<list<TaylorModelVec> >::const_iterator fpIter = flowpipesCompo.begin();
	list<list<vector<Interval> > >::const_iterator fpdoIter = domains.begin();
	list<int>::const_iterator modeIter = modeIDs.begin();

	list<TaylorModelVec>::const_iterator tmvIter;
	list<vector<Interval> >::const_iterator doIter;

	for(; fpIter!=flowpipesCompo.end(); ++fpIter, ++fpdoIter, ++modeIter)
	{
		tmvIter = fpIter->begin();
		doIter = fpdoIter->begin();

		for(; tmvIter != fpIter->end(); ++tmvIter, ++doIter)
		{
			vector<Interval> box;
			tmvIter->intEval(box, *doIter);

			// contract the interval according to the invariant
			vector<Interval> new_domain;
			vector<bool> bVecTemp;
			TaylorModelVec tmvInterval(box, new_domain);
			int type = contract_interval_arithmetic(tmvInterval, new_domain, system.invariants[*modeIter], bVecTemp);

			if(type < 0)
			{
				continue;
			}

			tmvInterval.intEval(box, new_domain);

			Interval X(box[outputAxes[0]]), Y(box[outputAxes[1]]);

			// output the vertices
			fprintf(fp, "%lf %lf\n", X.inf(), Y.inf());
			fprintf(fp, "%lf %lf\n", X.sup(), Y.inf());
			fprintf(fp, "%lf %lf\n", X.sup(), Y.sup());
			fprintf(fp, "%lf %lf\n", X.inf(), Y.sup());
			fprintf(fp, "%lf %lf\n", X.inf(), Y.inf());
			fprintf(fp, "\n\n");
		}
	}

	fprintf(fp, "e\n");
}

void HybridReachability::plot_2D_octagon_GNUPLOT(FILE *fp) const
{
	int x = outputAxes[0];
	int y = outputAxes[1];

	Interval intZero;

	int rangeDim = stateVarNames.size();
	Matrix output_poly_temp(8, rangeDim);

	output_poly_temp.set(1, 0, x);
	output_poly_temp.set(1, 1, y);
	output_poly_temp.set(-1, 2, x);
	output_poly_temp.set(-1, 3, y);
	output_poly_temp.set(1/sqrt(2), 4, x);
	output_poly_temp.set(1/sqrt(2), 4, y);
	output_poly_temp.set(1/sqrt(2), 5, x);
	output_poly_temp.set(-1/sqrt(2), 5, y);
	output_poly_temp.set(-1/sqrt(2), 6, x);
	output_poly_temp.set(1/sqrt(2), 6, y);
	output_poly_temp.set(-1/sqrt(2), 7, x);
	output_poly_temp.set(-1/sqrt(2), 7, y);

	// Construct the 2D template matrix.
	int rows = 8;
	int cols = rangeDim;

	Matrix sortedTemplate(rows, cols);
	RowVector rowVec(cols);
	list<RowVector> sortedRows;
	list<RowVector>::iterator iterp, iterq;

	output_poly_temp.row(rowVec, 0);
	sortedRows.push_back(rowVec);

	bool bInserted;

	// Sort the row vectors in the template by anti-clockwise order (only in the x-y space).
	for(int i=1; i<rows; ++i)
	{
		iterp = sortedRows.begin();
		iterq = iterp;
		++iterq;
		bInserted = false;

		for(; iterq != sortedRows.end();)
		{
			double tmp1 = output_poly_temp.get(i,x) * iterp->get(y) - output_poly_temp.get(i,y) * iterp->get(x);
			double tmp2 = output_poly_temp.get(i,x) * iterq->get(y) - output_poly_temp.get(i,y) * iterq->get(x);

			if(tmp1 < 0 && tmp2 > 0)
			{
				output_poly_temp.row(rowVec, i);
				sortedRows.insert(iterq, rowVec);
				bInserted = true;
				break;
			}
			else
			{
				++iterp;
				++iterq;
			}
		}

		if(!bInserted)
		{
			output_poly_temp.row(rowVec, i);
			sortedRows.push_back(rowVec);
		}
	}

	iterp = sortedRows.begin();
	for(int i=0; i<rows; ++i, ++iterp)
	{
		for(int j=0; j<cols; ++j)
		{
			sortedTemplate.set(iterp->get(j), i, j);
		}
	}

	ColVector b(rows);

	fprintf(fp, "set terminal postscript\n");

	char filename[NAME_SIZE+10];
	sprintf(filename, "%s%s.eps", imageDir, outputFileName);
	fprintf(fp, "set output '%s'\n", filename);

	fprintf(fp, "set style line 1 linecolor rgb \"blue\"\n");
	fprintf(fp, "set autoscale\n");
	fprintf(fp, "unset label\n");
	fprintf(fp, "set xtic auto\n");
	fprintf(fp, "set ytic auto\n");
	fprintf(fp, "set xlabel \"%s\"\n", stateVarNames[outputAxes[0]].c_str());
	fprintf(fp, "set ylabel \"%s\"\n", stateVarNames[outputAxes[1]].c_str());
	fprintf(fp, "plot '-' notitle with lines ls 1\n");

	vector<Interval> step_exp_table;
	Interval I(1);
	step_exp_table.push_back(I);

	// Compute the intersections of two facets.
	// The vertices are ordered clockwisely.

	gsl_matrix *C = gsl_matrix_alloc(2,2);
	gsl_vector *d = gsl_vector_alloc(2);
	gsl_vector *vertex = gsl_vector_alloc(2);

	list<list<TaylorModelVec> >::const_iterator fpIter = flowpipesCompo.begin();
	list<list<vector<Interval> > >::const_iterator fpdoIter = domains.begin();
	list<int>::const_iterator modeIter = modeIDs.begin();

	list<TaylorModelVec>::const_iterator tmvIter;
	list<vector<Interval> >::const_iterator doIter;

	for(; fpIter!=flowpipesCompo.end(); ++fpIter, ++fpdoIter, ++modeIter)
	{
		tmvIter = fpIter->begin();
		doIter = fpdoIter->begin();

		for(; tmvIter != fpIter->end(); ++tmvIter, ++doIter)
		{
			int rangeDim = tmvIter->tms.size();
			vector<Interval> box;
			tmvIter->intEval(box, *doIter);

			vector<Interval> new_domain;
			vector<bool> bVecTemp;
			TaylorModelVec tmvInterval(box, new_domain);
			int type = contract_interval_arithmetic(tmvInterval, new_domain, system.invariants[*modeIter], bVecTemp);

			if(type < 0)
			{
				continue;
			}

			tmvInterval.intEval(box, new_domain);

			// the box template
			b.set(box[x].sup(), 0);
			b.set(box[y].sup(), 2);
			b.set(-box[x].inf(), 4);
			b.set(-box[y].inf(), 6);

			// consider the other vectors
			Matrix other_vectors(rangeDim, rangeDim+1);
			for(int i=0; i<rangeDim; ++i)
			{
				if(i != x && i != y)
				{
					other_vectors.set(1, i, i+1);
				}
			}

			other_vectors.set(1/sqrt(2), x, x+1);
			other_vectors.set(1/sqrt(2), x, y+1);
			other_vectors.set(1/sqrt(2), y, x+1);
			other_vectors.set(-1/sqrt(2), y, y+1);

			TaylorModelVec tmv_other_vectors(other_vectors);

			for(int i=0; i<rangeDim; ++i)
			{

				RowVector rowVecTemp(rangeDim);

				for(int j=0; j<rangeDim; ++j)
				{
					rowVecTemp.set(other_vectors.get(i,j+1), j);
				}

				Interval intTemp = rho(*tmvIter, rowVecTemp, *doIter);
				new_domain[i+1].setSup(intTemp);

				rowVecTemp.neg_assign();

				intTemp = rho(*tmvIter, rowVecTemp, *doIter);
				intTemp.inv_assign();
				new_domain[i+1].setInf(intTemp);
			}

			type = contract_interval_arithmetic(tmv_other_vectors, new_domain, system.invariants[*modeIter], bVecTemp);

			if(type < 0)
			{
				continue;
			}


			RowVector template_vector_1(rangeDim);
			template_vector_1.set(1/sqrt(2), x);
			template_vector_1.set(1/sqrt(2), y);

			double sp = (rho(tmv_other_vectors, template_vector_1, new_domain)).sup();
			double sp2 = (1/sqrt(2))*b.get(0) + (1/sqrt(2))*b.get(2);
			b.set(sp, 1);


			RowVector template_vector_3(rangeDim);
			template_vector_3.set(-1/sqrt(2), x);
			template_vector_3.set(1/sqrt(2), y);

			sp = (rho(tmv_other_vectors, template_vector_3, new_domain)).sup();
			b.set(sp, 3);


			RowVector template_vector_5(rangeDim);
			template_vector_5.set(-1/sqrt(2), x);
			template_vector_5.set(-1/sqrt(2), y);

			sp = (rho(tmv_other_vectors, template_vector_5, new_domain)).sup();
			b.set(sp, 5);


			RowVector template_vector_7(rangeDim);
			template_vector_7.set(1/sqrt(2), x);
			template_vector_7.set(-1/sqrt(2), y);

			sp = (rho(tmv_other_vectors, template_vector_7, new_domain)).sup();
			b.set(sp, 7);

			Polyhedron polyTemplate(sortedTemplate, b);
			polyTemplate.tightenConstraints();

			double f1, f2;

			list<LinearConstraint>::iterator iterp, iterq;
			iterp = iterq = polyTemplate.constraints.begin();
			++iterq;

			for(; iterq != polyTemplate.constraints.end(); ++iterp, ++iterq)
			{
				gsl_matrix_set(C, 0, 0, iterp->A[x].midpoint());
				gsl_matrix_set(C, 0, 1, iterp->A[y].midpoint());
				gsl_matrix_set(C, 1, 0, iterq->A[x].midpoint());
				gsl_matrix_set(C, 1, 1, iterq->A[y].midpoint());

				gsl_vector_set(d, 0, iterp->B.midpoint());
				gsl_vector_set(d, 1, iterq->B.midpoint());

				gsl_linalg_HH_solve(C, d, vertex);

				double v1 = gsl_vector_get(vertex, 0);
				double v2 = gsl_vector_get(vertex, 1);

				if(iterp == polyTemplate.constraints.begin())
				{
					f1 = v1;
					f2 = v2;
				}

				fprintf(fp, "%lf %lf\n", v1, v2);
			}

			iterp = polyTemplate.constraints.begin();
			--iterq;

			gsl_matrix_set(C, 0, 0, iterp->A[x].midpoint());
			gsl_matrix_set(C, 0, 1, iterp->A[y].midpoint());
			gsl_matrix_set(C, 1, 0, iterq->A[x].midpoint());
			gsl_matrix_set(C, 1, 1, iterq->A[y].midpoint());

			gsl_vector_set(d, 0, iterp->B.midpoint());
			gsl_vector_set(d, 1, iterq->B.midpoint());

			gsl_linalg_HH_solve(C, d, vertex);

			double v1 = gsl_vector_get(vertex, 0);
			double v2 = gsl_vector_get(vertex, 1);

			fprintf(fp, "%lf %lf\n", v1, v2);

			fprintf(fp, "%lf %lf\n", f1, f2);
			fprintf(fp, "\n\n");
		}
	}

	fprintf(fp, "e\n");
}

void HybridReachability::plot_2D_grid_GNUPLOT(FILE *fp) const
{
	fprintf(fp, "set terminal postscript\n");

	char filename[NAME_SIZE+10];
	sprintf(filename, "%s%s.eps", imageDir, outputFileName);
	fprintf(fp, "set output '%s'\n", filename);

	fprintf(fp, "set style line 1 linecolor rgb \"blue\"\n");
	fprintf(fp, "set autoscale\n");
	fprintf(fp, "unset label\n");
	fprintf(fp, "set xtic auto\n");
	fprintf(fp, "set ytic auto\n");
	fprintf(fp, "set xlabel \"%s\"\n", stateVarNames[outputAxes[0]].c_str());
	fprintf(fp, "set ylabel \"%s\"\n", stateVarNames[outputAxes[1]].c_str());
	fprintf(fp, "plot '-' notitle with lines ls 1\n");

	list<list<TaylorModelVec> >::const_iterator fpIter = flowpipesCompo.begin();
	list<list<vector<Interval> > >::const_iterator fpdoIter = domains.begin();
	list<int>::const_iterator modeIter = modeIDs.begin();

	list<TaylorModelVec>::const_iterator tmvIter;
	list<vector<Interval> >::const_iterator doIter;

	vector<Interval> step_exp_table;
	Interval I(1);
	step_exp_table.push_back(I);

	for(; fpIter!=flowpipesCompo.end(); ++fpIter, ++fpdoIter, ++modeIter)
	{
		tmvIter = fpIter->begin();
		doIter = fpdoIter->begin();
		int domainDim = doIter->size();

		for(; tmvIter != fpIter->end(); ++tmvIter, ++doIter)
		{
			// decompose the domain
			list<vector<Interval> > grids;

			gridBox(grids, *doIter, numSections);

			// Transform the Taylor model into a Horner form
			vector<HornerForm> tmvHF;
			vector<Interval> remainders;
			int rangeDim = tmvIter->tms.size();

			for(int i=0; i<rangeDim; ++i)
			{
				HornerForm hfTemp;
				Interval intTemp;
				tmvIter->tms[i].toHornerForm(hfTemp, intTemp);
				tmvHF.push_back(hfTemp);
				remainders.push_back(intTemp);
			}

			// evaluate the images from all of the grids
			list<vector<Interval> >::const_iterator gIter = grids.begin();
			for(; gIter!=grids.end(); ++gIter)
			{
				vector<Interval> box;

				for(int i=0; i<rangeDim; ++i)
				{
					Interval intTemp;
					tmvHF[i].intEval(intTemp, *gIter);
					intTemp += remainders[i];
					box.push_back(intTemp);
				}

				// contract the interval according to the invariant
				vector<Interval> new_domain;
				vector<bool> bVecTemp;
				TaylorModelVec tmvInterval(box, new_domain);
				int type = contract_interval_arithmetic(tmvInterval, new_domain, system.invariants[*modeIter], bVecTemp);

				if(type < 0)
				{
					continue;
				}

				tmvInterval.intEval(box, new_domain);

				Interval X(box[outputAxes[0]]), Y(box[outputAxes[1]]);

				// output the vertices
				fprintf(fp, "%lf %lf\n", X.inf(), Y.inf());
				fprintf(fp, "%lf %lf\n", X.sup(), Y.inf());
				fprintf(fp, "%lf %lf\n", X.sup(), Y.sup());
				fprintf(fp, "%lf %lf\n", X.inf(), Y.sup());
				fprintf(fp, "%lf %lf\n", X.inf(), Y.inf());
				fprintf(fp, "\n\n");
			}
		}
	}

	fprintf(fp, "e\n");
}

void HybridReachability::plot_2D_MATLAB(FILE *fp) const
{
	switch(plotSetting)
	{
	case PLOT_INTERVAL:
		plot_2D_interval_MATLAB(fp);
		break;
	case PLOT_OCTAGON:
		plot_2D_octagon_MATLAB(fp);
		break;
	case PLOT_GRID:
		plot_2D_grid_MATLAB(fp);
		break;
	}
}

void HybridReachability::plot_2D_interval_MATLAB(FILE *fp) const
{
	list<list<TaylorModelVec> >::const_iterator fpIter = flowpipesCompo.begin();
	list<list<vector<Interval> > >::const_iterator fpdoIter = domains.begin();
	list<int>::const_iterator modeIter = modeIDs.begin();

	list<TaylorModelVec>::const_iterator tmvIter;
	list<vector<Interval> >::const_iterator doIter;

	for(; fpIter!=flowpipesCompo.end(); ++fpIter, ++fpdoIter, ++modeIter)
	{
		tmvIter = fpIter->begin();
		doIter = fpdoIter->begin();

		for(; tmvIter != fpIter->end(); ++tmvIter, ++doIter)
		{
			vector<Interval> box;
			tmvIter->intEval(box, *doIter);

			// contract the interval according to the invariant
			vector<Interval> new_domain;
			vector<bool> bVecTemp;
			TaylorModelVec tmvInterval(box, new_domain);
			int type = contract_interval_arithmetic(tmvInterval, new_domain, system.invariants[*modeIter], bVecTemp);

			if(type < 0)
			{
				continue;
			}

			tmvInterval.intEval(box, new_domain);

			Interval X(box[outputAxes[0]]), Y(box[outputAxes[1]]);

			// output the vertices
			fprintf(fp,"plot( [%lf,%lf,%lf,%lf,%lf] , [%lf,%lf,%lf,%lf,%lf] , 'b');\nhold on;\nclear;\n",
					X.inf(), X.sup(), X.sup(), X.inf(), X.inf(), Y.inf(), Y.inf(), Y.sup(), Y.sup(), Y.inf());
		}
	}
}

void HybridReachability::plot_2D_octagon_MATLAB(FILE *fp) const
{
	int x = outputAxes[0];
	int y = outputAxes[1];

	Interval intZero;

	int rangeDim = stateVarNames.size();
	Matrix output_poly_temp(8, rangeDim);

	output_poly_temp.set(1, 0, x);
	output_poly_temp.set(1, 1, y);
	output_poly_temp.set(-1, 2, x);
	output_poly_temp.set(-1, 3, y);
	output_poly_temp.set(1/sqrt(2), 4, x);
	output_poly_temp.set(1/sqrt(2), 4, y);
	output_poly_temp.set(1/sqrt(2), 5, x);
	output_poly_temp.set(-1/sqrt(2), 5, y);
	output_poly_temp.set(-1/sqrt(2), 6, x);
	output_poly_temp.set(1/sqrt(2), 6, y);
	output_poly_temp.set(-1/sqrt(2), 7, x);
	output_poly_temp.set(-1/sqrt(2), 7, y);

	// Construct the 2D template matrix.
	int rows = 8;
	int cols = rangeDim;

	Matrix sortedTemplate(rows, cols);
	RowVector rowVec(cols);
	list<RowVector> sortedRows;
	list<RowVector>::iterator iterp, iterq;

	output_poly_temp.row(rowVec, 0);
	sortedRows.push_back(rowVec);

	bool bInserted;

	// Sort the row vectors in the template by anti-clockwise order (only in the x-y space).
	for(int i=1; i<rows; ++i)
	{
		iterp = sortedRows.begin();
		iterq = iterp;
		++iterq;
		bInserted = false;

		for(; iterq != sortedRows.end();)
		{
			double tmp1 = output_poly_temp.get(i,x) * iterp->get(y) - output_poly_temp.get(i,y) * iterp->get(x);
			double tmp2 = output_poly_temp.get(i,x) * iterq->get(y) - output_poly_temp.get(i,y) * iterq->get(x);

			if(tmp1 < 0 && tmp2 > 0)
			{
				output_poly_temp.row(rowVec, i);
				sortedRows.insert(iterq, rowVec);
				bInserted = true;
				break;
			}
			else
			{
				++iterp;
				++iterq;
			}
		}

		if(!bInserted)
		{
			output_poly_temp.row(rowVec, i);
			sortedRows.push_back(rowVec);
		}
	}

	iterp = sortedRows.begin();
	for(int i=0; i<rows; ++i, ++iterp)
	{
		for(int j=0; j<cols; ++j)
		{
			sortedTemplate.set(iterp->get(j), i, j);
		}
	}

	ColVector b(rows);

	vector<Interval> step_exp_table;
	Interval I(1);
	step_exp_table.push_back(I);

	// Compute the intersections of two facets.
	// The vertices are ordered clockwisely.

	gsl_matrix *C = gsl_matrix_alloc(2,2);
	gsl_vector *d = gsl_vector_alloc(2);
	gsl_vector *vertex = gsl_vector_alloc(2);

	list<list<TaylorModelVec> >::const_iterator fpIter = flowpipesCompo.begin();
	list<list<vector<Interval> > >::const_iterator fpdoIter = domains.begin();
	list<int>::const_iterator modeIter = modeIDs.begin();

	list<TaylorModelVec>::const_iterator tmvIter;
	list<vector<Interval> >::const_iterator doIter;

	for(; fpIter!=flowpipesCompo.end(); ++fpIter, ++fpdoIter, ++modeIter)
	{
		tmvIter = fpIter->begin();
		doIter = fpdoIter->begin();

		for(; tmvIter != fpIter->end(); ++tmvIter, ++doIter)
		{
			int rangeDim = tmvIter->tms.size();
			vector<Interval> box;
			tmvIter->intEval(box, *doIter);

			vector<Interval> new_domain;
			vector<bool> bVecTemp;
			TaylorModelVec tmvInterval(box, new_domain);
			int type = contract_interval_arithmetic(tmvInterval, new_domain, system.invariants[*modeIter], bVecTemp);

			if(type < 0)
			{
				continue;
			}

			tmvInterval.intEval(box, new_domain);

			// the box template
			b.set(box[x].sup(), 0);
			b.set(box[y].sup(), 2);
			b.set(-box[x].inf(), 4);
			b.set(-box[y].inf(), 6);

			// consider the other vectors
			Matrix other_vectors(rangeDim, rangeDim+1);
			for(int i=0; i<rangeDim; ++i)
			{
				if(i != x && i != y)
				{
					other_vectors.set(1, i, i+1);
				}
			}

			other_vectors.set(1/sqrt(2), x, x+1);
			other_vectors.set(1/sqrt(2), x, y+1);
			other_vectors.set(1/sqrt(2), y, x+1);
			other_vectors.set(-1/sqrt(2), y, y+1);

			TaylorModelVec tmv_other_vectors(other_vectors);

			for(int i=0; i<rangeDim; ++i)
			{

				RowVector rowVecTemp(rangeDim);

				for(int j=0; j<rangeDim; ++j)
				{
					rowVecTemp.set(other_vectors.get(i,j+1), j);
				}

				Interval intTemp = rho(*tmvIter, rowVecTemp, *doIter);
				new_domain[i+1].setSup(intTemp);

				rowVecTemp.neg_assign();

				intTemp = rho(*tmvIter, rowVecTemp, *doIter);
				intTemp.inv_assign();
				new_domain[i+1].setInf(intTemp);
			}

			type = contract_interval_arithmetic(tmv_other_vectors, new_domain, system.invariants[*modeIter], bVecTemp);

			if(type < 0)
			{
				continue;
			}


			RowVector template_vector_1(rangeDim);
			template_vector_1.set(1/sqrt(2), x);
			template_vector_1.set(1/sqrt(2), y);

			double sp = (rho(tmv_other_vectors, template_vector_1, new_domain)).sup();
			double sp2 = (1/sqrt(2))*b.get(0) + (1/sqrt(2))*b.get(2);
			b.set(sp, 1);


			RowVector template_vector_3(rangeDim);
			template_vector_3.set(-1/sqrt(2), x);
			template_vector_3.set(1/sqrt(2), y);

			sp = (rho(tmv_other_vectors, template_vector_3, new_domain)).sup();
			b.set(sp, 3);


			RowVector template_vector_5(rangeDim);
			template_vector_5.set(-1/sqrt(2), x);
			template_vector_5.set(-1/sqrt(2), y);

			sp = (rho(tmv_other_vectors, template_vector_5, new_domain)).sup();
			b.set(sp, 5);


			RowVector template_vector_7(rangeDim);
			template_vector_7.set(1/sqrt(2), x);
			template_vector_7.set(-1/sqrt(2), y);

			sp = (rho(tmv_other_vectors, template_vector_7, new_domain)).sup();
			b.set(sp, 7);

			Polyhedron polyTemplate(sortedTemplate, b);
			polyTemplate.tightenConstraints();

			double f1, f2;

			list<LinearConstraint>::iterator iterp, iterq;
			iterp = iterq = polyTemplate.constraints.begin();
			++iterq;

			vector<double> vertices_x, vertices_y;

			for(; iterq != polyTemplate.constraints.end(); ++iterp, ++iterq)
			{
				gsl_matrix_set(C, 0, 0, iterp->A[x].midpoint());
				gsl_matrix_set(C, 0, 1, iterp->A[y].midpoint());
				gsl_matrix_set(C, 1, 0, iterq->A[x].midpoint());
				gsl_matrix_set(C, 1, 1, iterq->A[y].midpoint());

				gsl_vector_set(d, 0, iterp->B.midpoint());
				gsl_vector_set(d, 1, iterq->B.midpoint());

				gsl_linalg_HH_solve(C, d, vertex);

				double v1 = gsl_vector_get(vertex, 0);
				double v2 = gsl_vector_get(vertex, 1);

				if(iterp == polyTemplate.constraints.begin())
				{
					f1 = v1;
					f2 = v2;
				}

				vertices_x.push_back(v1);
				vertices_y.push_back(v2);
//				fprintf(fp, "%lf %lf\n", v1, v2);
			}

			iterp = polyTemplate.constraints.begin();
			--iterq;

			gsl_matrix_set(C, 0, 0, iterp->A[x].midpoint());
			gsl_matrix_set(C, 0, 1, iterp->A[y].midpoint());
			gsl_matrix_set(C, 1, 0, iterq->A[x].midpoint());
			gsl_matrix_set(C, 1, 1, iterq->A[y].midpoint());

			gsl_vector_set(d, 0, iterp->B.midpoint());
			gsl_vector_set(d, 1, iterq->B.midpoint());

			gsl_linalg_HH_solve(C, d, vertex);

			double v1 = gsl_vector_get(vertex, 0);
			double v2 = gsl_vector_get(vertex, 1);

			vertices_x.push_back(v1);
			vertices_y.push_back(v2);
			vertices_x.push_back(f1);
			vertices_y.push_back(f2);

			fprintf(fp, "plot( ");

			fprintf(fp, "[ ");
			for(int i=0; i<vertices_x.size()-1; ++i)
			{
				fprintf(fp, "%lf , ", vertices_x[i]);
			}
			fprintf(fp, "%lf ] , ", vertices_x.back());

			fprintf(fp, "[ ");
			for(int i=0; i<vertices_y.size()-1; ++i)
			{
				fprintf(fp, "%lf , ", vertices_y[i]);
			}
			fprintf(fp, "%lf ] , ", vertices_y.back());

			fprintf(fp, "'b');\nhold on;\nclear;\n");

//			fprintf(fp, "%lf %lf\n", v1, v2);
//			fprintf(fp, "%lf %lf\n", f1, f2);
//			fprintf(fp, "\n\n");
		}
	}
}

void HybridReachability::plot_2D_grid_MATLAB(FILE *fp) const
{
	list<list<TaylorModelVec> >::const_iterator fpIter = flowpipesCompo.begin();
	list<list<vector<Interval> > >::const_iterator fpdoIter = domains.begin();
	list<int>::const_iterator modeIter = modeIDs.begin();

	list<TaylorModelVec>::const_iterator tmvIter;
	list<vector<Interval> >::const_iterator doIter;

	vector<Interval> step_exp_table;
	Interval I(1);
	step_exp_table.push_back(I);

	for(; fpIter!=flowpipesCompo.end(); ++fpIter, ++fpdoIter, ++modeIter)
	{
		tmvIter = fpIter->begin();
		doIter = fpdoIter->begin();
		int domainDim = doIter->size();

		for(; tmvIter != fpIter->end(); ++tmvIter, ++doIter)
		{
			// decompose the domain
			list<vector<Interval> > grids;

			gridBox(grids, *doIter, numSections);

			// Transform the Taylor model into a Horner form
			vector<HornerForm> tmvHF;
			vector<Interval> remainders;
			int rangeDim = tmvIter->tms.size();

			for(int i=0; i<rangeDim; ++i)
			{
				HornerForm hfTemp;
				Interval intTemp;
				tmvIter->tms[i].toHornerForm(hfTemp, intTemp);
				tmvHF.push_back(hfTemp);
				remainders.push_back(intTemp);
			}

			// evaluate the images from all of the grids
			list<vector<Interval> >::const_iterator gIter = grids.begin();
			for(; gIter!=grids.end(); ++gIter)
			{
				vector<Interval> box;

				for(int i=0; i<rangeDim; ++i)
				{
					Interval intTemp;
					tmvHF[i].intEval(intTemp, *gIter);
					intTemp += remainders[i];
					box.push_back(intTemp);
				}

				// contract the interval according to the invariant
				vector<Interval> new_domain;
				vector<bool> bVecTemp;
				TaylorModelVec tmvInterval(box, new_domain);
				int type = contract_interval_arithmetic(tmvInterval, new_domain, system.invariants[*modeIter], bVecTemp);

				if(type < 0)
				{
					continue;
				}

				tmvInterval.intEval(box, new_domain);

				Interval X(box[outputAxes[0]]), Y(box[outputAxes[1]]);

				// output the vertices
				fprintf(fp,"plot( [%lf,%lf,%lf,%lf,%lf] , [%lf,%lf,%lf,%lf,%lf] , 'b');\nhold on;\nclear;\n",
						X.inf(), X.sup(), X.sup(), X.inf(), X.inf(), Y.inf(), Y.inf(), Y.sup(), Y.sup(), Y.inf());
			}
		}
	}
}

bool HybridReachability::declareStateVar(const string & vName)
{
	map<string,int>::const_iterator iter;

	if((iter = stateVarTab.find(vName)) == stateVarTab.end())
	{
		stateVarTab[vName] = stateVarNames.size();
		stateVarNames.push_back(vName);
		return true;
	}
	else
	{
		return false;
	}
}

int HybridReachability::getIDForStateVar(const string & vName) const
{
	map<string,int>::const_iterator iter;
	if((iter = stateVarTab.find(vName)) == stateVarTab.end())
	{
		return -1;
	}

	return iter->second;
}

bool HybridReachability::getStateVarName(string & vName, int id) const
{
	if(id>=0 && id<stateVarNames.size())
	{
		vName = stateVarNames[id];
		return true;
	}
	else
	{
		return false;
	}
}

bool HybridReachability::declareTMVar(const string & vName)
{
	map<string,int>::const_iterator iter;

	if((iter = tmVarTab.find(vName)) == tmVarTab.end())
	{
		tmVarTab[vName] = tmVarNames.size();
		tmVarNames.push_back(vName);
		return true;
	}
	else
	{
		return false;
	}
}

int HybridReachability::getIDForTMVar(const string & vName) const
{
	map<string,int>::const_iterator iter;
	if((iter = tmVarTab.find(vName)) == tmVarTab.end())
	{
		return -1;
	}

	return iter -> second;
}

bool HybridReachability::getTMVarName(string & vName, const int id) const
{
	if(id>=0 && id<tmVarNames.begin()->size())
	{
		vName = (*tmVarNames.begin())[id];
		return true;
	}
	else
	{
		return false;
	}
}

bool HybridReachability::declareMode(const string & mName, const TaylorModelVec & ode, const vector<Interval> & uncertainties, const vector<PolynomialConstraint> & inv, const int integrationScheme)
{
	map<string,int>::const_iterator iter;
	Interval intZero;

	int rangeDim = uncertainties.size();
	TaylorModelVec tmvEmpty;
	vector<HornerForm> hfsEmpty;
	vector<Interval> uncertainty_centers_empty;

	if((iter = modeTab.find(mName)) == modeTab.end())
	{
		int modeID = modeNames.size();
		modeTab[mName] = modeID;
		modeNames.push_back(mName);

		system.modes.push_back(modeID);

		if(modeID > system.odes.size())
		{
			for(int i=system.odes.size(); i<modeID; ++i)
			{
				system.odes.push_back(tmvEmpty);
				system.hfOdes.push_back(hfsEmpty);
				system.uncertainty_centers.push_back(uncertainty_centers_empty);
			}
		}

		TaylorModelVec tmvTemp = ode;

		vector<HornerForm> hfOde;
		vector<Interval> uncertainties_centered = uncertainties;

		for(int i=0; i<uncertainties_centered.size(); ++i)
		{
			Interval M;
			uncertainties_centered[i].remove_midpoint(M);

			if(!M.subseteq(intZero))
			{
				TaylorModel tmTemp(M, rangeDim+1);
				tmvTemp.tms[i].add_assign(tmTemp);
			}

			HornerForm hf;
			tmvTemp.tms[i].expansion.toHornerForm(hf);
			hfOde.push_back(hf);
		}

		system.odes.push_back(tmvTemp);
		system.hfOdes.push_back(hfOde);
		system.invariants.push_back(inv);
		system.uncertainties.push_back(uncertainties_centered);

		integrationSchemes.push_back(integrationScheme);

		return true;
	}
	else
	{
		return false;
	}
}

bool HybridReachability::declareMode(const string & mName, const vector<string> & strOde, const vector<Interval> & uncertainties, const vector<PolynomialConstraint> & inv, const int integrationScheme)
{
	map<string,int>::const_iterator iter;
	Interval intZero;

	vector<string> odeEmpty;

	if((iter = modeTab.find(mName)) == modeTab.end())
	{
		int modeID = modeNames.size();
		modeTab[mName] = modeID;
		modeNames.push_back(mName);

		system.modes.push_back(modeID);

		if(modeID > system.strOdes.size())
		{
			for(int i=system.strOdes.size(); i<modeID; ++i)
			{
				system.strOdes.push_back(odeEmpty);
			}
		}

		system.strOdes.push_back(strOde);
		system.invariants.push_back(inv);

		vector<Interval> uncertainties_centered = uncertainties;
		vector<Interval> uncertainty_centers;

		for(int i=0; i<uncertainties_centered.size(); ++i)
		{
			Interval M;
			uncertainties_centered[i].remove_midpoint(M);
			uncertainty_centers.push_back(M);
		}

		system.uncertainties.push_back(uncertainties_centered);
		system.uncertainty_centers.push_back(uncertainty_centers);

		integrationSchemes.push_back(integrationScheme);

		return true;
	}
	else
	{
		return false;
	}
}

int HybridReachability::getIDForMode(const string & mName) const
{
	map<string,int>::const_iterator iter;
	if((iter = modeTab.find (mName)) == modeTab.end())
	{
		return -1;
	}

	return iter -> second;
}

bool HybridReachability::getModeName(string & mName, int id) const
{
	if(id>=0 && id<modeNames.size())
	{
		mName = modeNames[id];
		return true;
	}
	else
	{
		return false;
	}
}

void HybridReachability::declareTrans(const int start, const int end, const vector<PolynomialConstraint> & guard, const ResetMap & reset, const int aggregType, const vector<vector<double> > & candidates)
{
	DiscTrans transition(++numOfJumps, start, end, guard, reset);

	system.transitions[start].push_back(transition);
	aggregationType[start][end] = aggregType;

	switch(aggregType)
	{
	case INTERVAL_AGGREG:
		// we do nothing
		break;
	case PARA_AGGREG:
	{
		vector<RowVector> rowVecs;
		for(int i=0; i<candidates.size(); ++i)
		{
			RowVector rowVec(candidates[i].size());

			for(int j=0; j<candidates[i].size(); ++j)
			{
				rowVec.set(candidates[i][j], j);
			}

			rowVec.normalize();
			rowVecs.push_back(rowVec);
		}

		aggregationTemplate_candidates[start][end] = rowVecs;
		break;
	}
	}
}

void HybridReachability::declareTrans()
{
	vector<DiscTrans> transVec;

	for(int i=0; i<system.modes.size(); ++i)
	{
		system.transitions.push_back(transVec);
	}

	// template types
	vector<int> iVec;

	for(int i=0; i<system.modes.size(); ++i)
	{
		iVec.push_back(-1);
	}

	for(int i=0; i<system.modes.size(); ++i)
	{
		aggregationType.push_back(iVec);
	}

	// parallelotopic templates
	vector<RowVector> paraTemplate;
	vector<vector<RowVector> > paraTemplateVec;

	for(int i=0; i<system.modes.size(); ++i)
	{
		paraTemplateVec.push_back(paraTemplate);
	}

	for(int i=0; i<system.modes.size(); ++i)
	{
		aggregationTemplate_candidates.push_back(paraTemplateVec);
	}

	// Taylor model templates
	TaylorModelVec tmvTemp;
	vector<TaylorModelVec> tmvVec;

	for(int i=0; i<system.modes.size(); ++i)
	{
		tmvVec.push_back(tmvTemp);
	}

	for(int i=0; i<system.modes.size(); ++i)
	{
		aggregationTemplate_TaylorModel.push_back(tmvVec);
	}
}

void HybridReachability::initialConfig(const int modeID, const Flowpipe & initialSet)
{
	system.initialMode = modeID;
	system.initialSet = initialSet;
}

void HybridReachability::set_default_template()
{
	int rangeDim = system.initialSet.tmv.tms.size();

	default_aggregation_template.clear();

	for(int i=0; i<rangeDim; ++i)
	{
		RowVector rowVec(rangeDim);
		rowVec.set(1, i);
		default_aggregation_template.push_back(rowVec);
	}

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=i+1; j<rangeDim; ++j)
		{
			RowVector rowVec(rangeDim);
			rowVec.set(1/sqrt(2), i);
			rowVec.set(1/sqrt(2), j);
			default_aggregation_template.push_back(rowVec);
		}
	}

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=i+1; j<rangeDim; ++j)
		{
			RowVector rowVec(rangeDim);
			rowVec.set(1/sqrt(2), i);
			rowVec.set(-1/sqrt(2), j);
			default_aggregation_template.push_back(rowVec);
		}
	}
}

void HybridReachability::constructWeightTab()
{
	int num_modes = system.modes.size();
	int rangeDim = default_aggregation_template[0].size();

	// construct the template information table
	vector<RowVector> rowVecs_empty;
	vector<bool> bVec_empty;
	vector<vector<bool> > bVecVec_temp;
	vector<vector<RowVector> > rowVecVecs_temp;

	for(int i=0; i<num_modes; ++i)
	{
		bVecVec_temp.push_back(bVec_empty);
		rowVecVecs_temp.push_back(rowVecs_empty);
	}

	for(int i=0; i<num_modes; ++i)
	{
		linear_auto.push_back(bVecVec_temp);
		template_auto.push_back(rowVecVecs_temp);
	}

	// construct the default weight table
	vector<Matrix> matVec;
	Matrix matDummy(1);
	for(int i=0; i<num_modes; ++i)
	{
		matVec.push_back(matDummy);
	}

	for(int i=0; i<num_modes; ++i)
	{
		weightTab.push_back(matVec);
	}

	// construct the weighted table for default template
	int num_default = default_aggregation_template.size();
	Matrix defaultWeightTab(num_default, num_default);
	for(int i=0; i<num_default; ++i)
	{
		for(int j=i+1; j<num_default; ++j)
		{
			// we assume that all the candidates and default template vectors are normalized
			double weight = 1 - fabs(default_aggregation_template[i].innerProd(default_aggregation_template[j]));
			defaultWeightTab.set(weight, i, j);
			defaultWeightTab.set(weight, j, i);
		}
	}

	// construct the weighted table for every jump
	for(int i=0; i<num_modes; ++i)
	{
		int num_invariants = system.invariants[i].size();

		// extract the normal vectors from the linear invaiant constraints
		vector<bool> linear_inv;
		vector<RowVector> template_inv;

		for(int j=0; j<num_invariants; ++j)
		{
			if(system.invariants[i][j].p.degree() == 1)
			{
				RowVector rowVec(rangeDim);
				system.invariants[i][j].p.constraintCoefficients(rowVec);
				rowVec.normalize();
				template_inv.push_back(rowVec);
				linear_inv.push_back(true);
			}
			else
			{
				RowVector rowVec(rangeDim);
				template_inv.push_back(rowVec);
				linear_inv.push_back(false);
			}
		}

		Matrix tabInv(1, 1);

		if(num_invariants > 0)
		{
			// table for invariant constraints
			Matrix tabInv_update(num_invariants, num_invariants);
			tabInv = tabInv_update;

			for(int j=0; j<num_invariants; ++j)
			{
				for(int k=j+1; k<num_invariants; ++k)
				{
					if(linear_inv[j] && linear_inv[k])
					{
						double weight = 1 - fabs(template_inv[j].innerProd(template_inv[k]));
						tabInv.set(weight, j, k);
						tabInv.set(weight, k, j);
					}
					else
					{
						tabInv.set(-1, j, k);
						tabInv.set(-1, k, j);
					}
				}
			}

	//		tabInv.output(stdout);
		}

		for(int j=0; j<system.transitions[i].size(); ++j)
		{
			int targetID = system.transitions[i][j].targetID;

			if(aggregationType[i][targetID] != PARA_AGGREG)
			{
				continue;
			}

			int num_candidates = aggregationTemplate_candidates[i][targetID].size();
			int auto_start = num_candidates;
			int guard_start = auto_start + num_invariants;
			int default_start = guard_start + system.transitions[i][j].guard.size();
			int weightMatSize = default_start + num_default;

//			printf("auto_start: %d,\tguard_start: %d,\tflow_start: %d,\tdefault_start: %d,\tdefault size: %d,\tmatrix size: %d\n", auto_start, guard_start, flow_start, default_start, num_default, weightMatSize);

			linear_auto[i][targetID] = linear_inv;
			template_auto[i][targetID] = template_inv;

			for(int k=0; k<system.transitions[i][j].guard.size(); ++k)
			{
				if(system.transitions[i][j].guard[k].p.degree() == 1)
				{
					RowVector rowVec(rangeDim);
					system.transitions[i][j].guard[k].p.constraintCoefficients(rowVec);
					rowVec.normalize();
					template_auto[i][targetID].push_back(rowVec);
					linear_auto[i][targetID].push_back(true);
				}
				else
				{
					RowVector rowVec(rangeDim);
					template_auto[i][targetID].push_back(rowVec);
					linear_auto[i][targetID].push_back(false);
				}
			}

			Matrix weightMat(weightMatSize, weightMatSize);
			// construct the whole weighted table
			for(int k=0; k<weightMatSize; ++k)
			{
				for(int m=k+1; m<weightMatSize; ++m)
				{
					if(k < auto_start)
					{
						// the first vector is a candidate

						if(m < auto_start)
						{
							// the second vector is a candidate
							double weight = 1 - fabs(aggregationTemplate_candidates[i][targetID][k].innerProd(aggregationTemplate_candidates[i][targetID][m]));
							weightMat.set(weight, k, m);
							weightMat.set(weight, m, k);
						}
						else if(m >= auto_start && m < default_start)
						{
							// the second vector is from the guard
							int posm = m - auto_start;
							if(linear_auto[i][targetID][posm])
							{
								double weight = 1 - fabs(aggregationTemplate_candidates[i][targetID][k].innerProd(template_auto[i][targetID][posm]));
								weightMat.set(weight, k, m);
								weightMat.set(weight, m, k);
							}
							else
							{
								weightMat.set(-1, k, m);
								weightMat.set(-1, m, k);
							}
						}
						else
						{
							// the second vector is from the default template
							double weight = 1 - fabs(aggregationTemplate_candidates[i][targetID][k].innerProd(default_aggregation_template[m - default_start]));
							weightMat.set(weight, k, m);
							weightMat.set(weight, m, k);
						}
					}
					else if(k >= auto_start && k < guard_start)
					{
						// the first vector is from the invariant
						int posk = k - auto_start;

						if(m >= auto_start && m < guard_start)
						{
							// the second vector is from the invariant
							int posm = m - auto_start;
							weightMat.set(tabInv.get(posk, posm), k, m);
							weightMat.set(tabInv.get(posm, posk), m, k);
						}
						else if(m >= guard_start && m < default_start)
						{
							// the second vector is from the guard
							int posm = m - auto_start;
							if(linear_auto[i][targetID][posk] && linear_auto[i][targetID][posm])
							{
								double weight = 1 - fabs(template_auto[i][targetID][posk].innerProd(template_auto[i][targetID][posm]));
								weightMat.set(weight, k, m);
								weightMat.set(weight, m, k);
							}
							else
							{
								weightMat.set(-1, k, m);
								weightMat.set(-1, m, k);
							}
						}
						else
						{
							// the second vector is from the default template
							if(linear_auto[i][targetID][posk])
							{
								double weight = 1 - fabs(template_auto[i][targetID][posk].innerProd(default_aggregation_template[m - default_start]));
								weightMat.set(weight, k, m);
								weightMat.set(weight, m, k);
							}
							else
							{
								weightMat.set(-1, k, m);
								weightMat.set(-1, m, k);
							}
						}
					}
					else if(k >= guard_start && k < default_start)
					{
						// the first vector is from the guard
						int posk = k - auto_start;

						if(m >= guard_start && m < default_start)
						{
							// the second vector is from the guard
							int posm = m - auto_start;
							if(linear_auto[i][targetID][posk] && linear_auto[i][targetID][posm])
							{
								double weight = 1 - fabs(template_auto[i][targetID][posk].innerProd(template_auto[i][targetID][posm]));
								weightMat.set(weight, k, m);
								weightMat.set(weight, m, k);
							}
							else
							{
								weightMat.set(-1, k, m);
								weightMat.set(-1, m, k);
							}
						}
						else
						{
							// the second vector is from the default template
							if(linear_auto[i][targetID][posk])
							{
								double weight = 1 - fabs(template_auto[i][targetID][posk].innerProd(default_aggregation_template[m - default_start]));
								weightMat.set(weight, k, m);
								weightMat.set(weight, m, k);
							}
							else
							{
								weightMat.set(-1, k, m);
								weightMat.set(-1, m, k);
							}
						}
					}
					else
					{
						// the first vector is from the default template
						int posk = k - default_start;
						int posm = m - default_start;

						double weight = 1 - fabs(default_aggregation_template[posk].innerProd(default_aggregation_template[posm]));
						weightMat.set(weight, k, m);
						weightMat.set(weight, m, k);
					}
				}
			}

			weightTab[i][targetID] = weightMat;

//			weightMat.output(stdout); exit(0);// test work
		}
	}
}

int HybridReachability::safetyChecking()
{
	int rangeDim = flowpipesCompo.front().front().tms.size();
	Interval intZero;

	bool bDumpCounterexamples = true;

	int mkres = mkdir(counterexampleDir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for counterexamples.\n");
		bDumpCounterexamples = false;
	}

	char filename_counterexamples[NAME_SIZE+10];
	FILE *fpDumpCounterexamples;

	if(bDumpCounterexamples)
	{
		sprintf(filename_counterexamples, "%s%s%s", counterexampleDir, outputFileName, str_counterexample_dumping_name_suffix);
		fpDumpCounterexamples = fopen(filename_counterexamples, "w");
	}

	// we first check the intersection of the linear invariant constraints and linear unsafe constraints
	for(int i=0; i<unsafeSet.size(); ++i)
	{
		list<LinearConstraint> lcs;

		for(int j=0; j<system.invariants[i].size(); ++j)
		{
			if(system.invariants[i][j].p.degree() == 1)
			{
				vector<Interval> intRowVec;

				for(int k=0; k<rangeDim; ++k)
				{
					intRowVec.push_back(intZero);
				}

				system.invariants[i][j].p.constraintCoefficients(intRowVec);
				LinearConstraint lc(intRowVec, system.invariants[i][j].B);
				lcs.push_back(lc);
			}
		}

		for(int j=0; j<unsafeSet[i].size(); ++j)
		{
			if(unsafeSet[i][j].p.degree() == 1)
			{
				vector<Interval> intRowVec;

				for(int k=0; k<rangeDim; ++k)
				{
					intRowVec.push_back(intZero);
				}

				unsafeSet[i][j].p.constraintCoefficients(intRowVec);
				LinearConstraint lc(intRowVec, unsafeSet[i][j].B);
				lcs.push_back(lc);
			}
		}

		Polyhedron P(lcs);

		if(bVecUnderCheck[i] && P.empty())
		{
			bVecUnderCheck[i] = false;
		}
	}

	// the main procedure of safety checking
	list<list<TaylorModelVec> >::const_iterator fpIter = flowpipesCompo.begin();
	list<list<vector<Interval> > >::const_iterator fpdoIter = domains.begin();
	list<int>::const_iterator modeIter = modeIDs.begin();
	list<TreeNode *>::const_iterator nodeIter = traceNodes.begin();

	vector<Interval> step_exp_table;

	list<TaylorModelVec>::const_iterator tmvIter;
	list<vector<Interval> >::const_iterator doIter;

	int maxOrder = 0;
	int result = SAFE;

	for(; fpIter!=flowpipesCompo.end(); ++fpIter, ++fpdoIter, ++modeIter, ++nodeIter)
	{
		if(!bVecUnderCheck[*modeIter])
		{
			continue;
		}

		tmvIter = fpIter->begin();
		doIter = fpdoIter->begin();

		if(unsafeSet[*modeIter].size() == 0)
		{
			return UNSAFE;
		}

		int domainDim = doIter->size();

		list<TaylorModelVec> flowpipe_counterexamples;
		list<vector<Interval> > counterexample_domains;
		list<Interval> localTimes;
		Interval localTime;

		for(; tmvIter!=fpIter->end(); ++tmvIter, ++doIter)
		{

			int tmp = maxOrder;
			for(int i=0; i<tmvIter->tms.size(); ++i)
			{
				int order = tmvIter->tms[i].expansion.degree();
				if(maxOrder < order)
				{
					maxOrder = order;
				}
			}

			if(step_exp_table.size() == 0 || step_exp_table[1] != (*doIter)[0] || maxOrder > tmp)
			{
				construct_step_exp_table(step_exp_table, (*doIter)[0], 2*maxOrder);
			}

			bool bsafe = false;

			vector<Interval> tmvPolyRange;
			tmvIter->polyRangeNormal(tmvPolyRange, step_exp_table);

			for(int i=0; i<unsafeSet[*modeIter].size(); ++i)
			{
				TaylorModel tmTemp;

				// interval evaluation on the constraint
				unsafeSet[*modeIter][i].hf.insert_normal(tmTemp, *tmvIter, tmvPolyRange, step_exp_table, domainDim);

				Interval intTemp;
				tmTemp.intEvalNormal(intTemp, step_exp_table);

				if(intTemp > unsafeSet[*modeIter][i].B)
				{
					// no intersection with the unsafe set
					bsafe = true;
					break;
				}
				else
				{
					continue;
				}
			}

			if(!bsafe)
			{
				// collect the skeptical counterexamples

				if(bDumpCounterexamples)
				{
					flowpipe_counterexamples.push_back(*tmvIter);
					counterexample_domains.push_back(*doIter);
					localTimes.push_back(localTime);
				}

				result = UNKNOWN;
			}

			localTime += (*doIter)[0];
		}

		if(bDumpCounterexamples && flowpipe_counterexamples.size() > 0)
		{
			// dump the skeptical counterexamples
			dump_potential_counterexample(fpDumpCounterexamples, flowpipe_counterexamples, counterexample_domains, *nodeIter, localTimes);
		}
	}

	if(bDumpCounterexamples)
	{
		fclose(fpDumpCounterexamples);
	}

	return result;
}

unsigned long HybridReachability::numOfFlowpipes() const
{
	list<list<TaylorModelVec> >::const_iterator fpIter = flowpipesCompo.begin();
	unsigned long sum = 0;

	for(; fpIter!=flowpipesCompo.end(); ++fpIter)
	{
		sum += (unsigned long)(fpIter->size());
	}

	return (sum-1);
}

void HybridReachability::dump_potential_counterexample(FILE *fp, const list<TaylorModelVec> & flowpipes, const list<vector<Interval> > & domains, TreeNode * const node, const list<Interval> & localTimes) const
{
	// dump the flowpipes
	fprintf(fp, "%s\n{\n", modeNames[node->modeID].c_str());

	list<TaylorModelVec>::const_iterator fpIter = flowpipes.begin();
	list<vector<Interval> >::const_iterator doIter = domains.begin();
	list<Interval>::const_iterator timeIter = localTimes.begin();

	for(; fpIter!=flowpipes.end(); ++fpIter, ++doIter, ++timeIter)
	{
		fprintf(fp, "starting time %lf\n{\n", timeIter->sup());

		fpIter->dump_interval(fp, stateVarNames, tmVarNames);

		for(int i=0; i<doIter->size(); ++i)
		{
			fprintf(fp, "%s in ", tmVarNames[i].c_str());
			(*doIter)[i].dump(fp);
			fprintf(fp, "\n");
		}

		fprintf(fp, "}\n\n");
	}

	fprintf(fp, "computation path\n{\n\n");

	TreeNode *iterator = node;

	string strTrace;
	char buffer[NAME_SIZE];

	for(; iterator->parent != NULL; iterator = iterator->parent)
	{
		string strLocalTime;
		iterator->localTime.toString(strLocalTime);

		sprintf(buffer, "( %d ", iterator->jumpID);
		string strJumpID(buffer);

		strTrace = ' ' + strJumpID + ',' + ' ' + strLocalTime + ')' + ' ' + '-' + '>' + ' ' + modeNames[iterator->modeID] + strTrace;
	}

	strTrace = modeNames[iterator->modeID] + strTrace;

	fprintf(fp, "%s;\n\n}\n\n}\n\n", strTrace.c_str());
}

































// class FactorTab

FactorTab::FactorTab()
{
	index = 0;
}

FactorTab::FactorTab(const int index_input, const Interval & factor_input, const Interval & intercept_input)
{
	index = index_input;
	factor = factor_input;
	intercept = intercept_input;
}

FactorTab::~FactorTab()
{
}

bool compareFactor(const FactorTab & a, const FactorTab & b)
{
	if(a.factor > b.factor)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool compareIntercept(const FactorTab & a, const FactorTab & b)
{
	if(a.intercept < b.intercept)
	{
		return true;
	}
	else
	{
		return false;
	}
}























void generateNodeSeq(list<TreeNode *> & result, TreeNode *root)
{
	list<TreeNode *> queue;
	queue.push_back(root);

	result.clear();

	for(; queue.size() > 0;)
	{
		TreeNode *current = queue.front();
		queue.pop_front();

		result.push_back(current);

		if(current->children.size() == 0)
		{
			continue;
		}

		list<TreeNode *>::iterator iter = current->children.begin();
		for(; iter != current->children.end(); ++iter)
		{
			queue.push_back(*iter);
		}
	}
}


































void aggregate_flowpipes_by_interval(TaylorModelVec & tmvAggregation, vector<Interval> & doAggregation, const vector<TaylorModelVec> & flowpipes, const vector<vector<Interval> > & domains)
{
	vector<Interval> intVecAggregation;

	flowpipes[0].intEval(intVecAggregation, domains[0]);

	Interval I1, I2;

	for(int k=0; k<flowpipes.size(); ++k)
	{
		vector<Interval> intVecTemp;
		flowpipes[k].intEval(intVecTemp, domains[k]);

		for(int i=0; i<intVecAggregation.size(); ++i)
		{

			intVecAggregation[i].inf(I1);
			intVecTemp[i].inf(I2);

			if(I2 <= I1)
			{
				intVecAggregation[i].setInf(I2);
			}

			intVecAggregation[i].sup(I1);
			intVecTemp[i].sup(I2);

			if(I2 >= I1)
			{
				intVecAggregation[i].setSup(I2);
			}
		}
	}

	TaylorModelVec tmvTemp(intVecAggregation, doAggregation);

	tmvAggregation = tmvTemp;
}

bool vector_selection(FactorTab & lst_selected, list<FactorTab> & lst_unselected, Matrix & matTemplate, const vector<RowVector> & rowVecs, int & rank)
{
	lst_unselected.sort(compareFactor);

	list<FactorTab> candidates;
	bool bvalid = false;


/*
	// =============== test begin ==================
	printf("Candidates:\n");
	list<FactorTab>::iterator testIter = lst_unselected.begin();
	for(; testIter!=lst_unselected.end(); ++testIter)
	{
		printf("vector: ");
		rowVecs[testIter->index].dump(stdout);
		printf("\tintercept: %lf, factor: %lf\n\n", testIter->intercept.midpoint(), testIter->factor.midpoint());
	}
	// =============== test end ==================
*/


	Interval intZero;

	for(; !bvalid && lst_unselected.size()!=0;)
	{
		list<FactorTab>::iterator facIter = lst_unselected.begin();

		Interval factor = facIter->factor;
		candidates.push_back(*facIter);
		facIter = lst_unselected.erase(facIter);

		for(; facIter!=lst_unselected.end(); )
		{
			if(facIter->factor < intZero)
			{
				facIter = lst_unselected.erase(facIter);
			}
			else if(facIter->factor.within(factor, THRESHOLD_LOW))
			{
				candidates.push_back(*facIter);
				facIter = lst_unselected.erase(facIter);
			}
			else
			{
				break;
			}
		}

		for(; !bvalid && candidates.size()!=0; )
		{
			facIter = candidates.begin();
			Interval min_intercept = facIter->intercept;
			list<FactorTab>::iterator iter_selected = facIter;

			++facIter;

			for(; facIter!=candidates.end(); ++facIter)
			{
				if(min_intercept > facIter->intercept)
				{
					min_intercept = facIter->intercept;
					iter_selected = facIter;
				}
			}

			lst_selected = *iter_selected;
			candidates.erase(iter_selected);

			bvalid = check_validity(matTemplate, rowVecs[lst_selected.index], rank);
		}
	}

	if(bvalid)
	{
		++rank;

		// insert the unselected elements back
		list<FactorTab>::iterator facIter = candidates.begin();

		for(; facIter!=candidates.end(); ++facIter)
		{
			lst_unselected.push_back(*facIter);
		}


/*
		// =============== test begin ==================
		printf("selected: ");
		rowVecs[lst_selected.index].dump(stdout);
		printf("\tintercept: %lf, factor: %lf\n\n\n", lst_selected.intercept.midpoint(), lst_selected.factor.midpoint());
		// =============== test end ==================
*/




		return true;	// one element is selected
	}
	else
	{
		return false;	// nothing in the unselected list is selectable
	}
}

bool check_validity(Matrix & matTemplate, const RowVector & rowVec, const int rank)
{
	int num = rowVec.size();
	for(int i=0; i<num; ++i)
	{
		matTemplate.set(rowVec.get(i), rank, i);
	}

	int r = matTemplate.rank();

	if(r == rank+1)
	{
		return true;
	}
	else
	{
		return false;
	}
}

void aggregate_flowpipes_by_Parallelotope(TaylorModelVec & tmvAggregation, vector<Interval> & doAggregation, const vector<TaylorModelVec> & flowpipes,
		const vector<vector<Interval> > & domains, const vector<PolynomialConstraint> & invariant, const DiscTrans & jump, vector<bool> & boundary_intersected,
		const vector<RowVector> & template_candidates, const vector<RowVector> & template_default, const vector<vector<Matrix> > & weightTab,
		const vector<vector<vector<bool> > > & linear_auto, const vector<vector<vector<RowVector> > > & template_auto, const int globalMaxOrder, const int rangeDim)
{
	int num_invariants = invariant.size();
	int num_candidates = template_candidates.size();

	int auto_start = num_candidates;
	int guard_start = auto_start + num_invariants;
	int default_start = guard_start + jump.guard.size();

	ColVector parallelotope_b(2*rangeDim);

	int startID = jump.startID;
	int targetID = jump.targetID;
	Matrix weightMat = weightTab[startID][targetID];

	vector<Interval> rhoPos;
	vector<Interval> rhoNeg;

	Interval intZero, intInvalid(INVALID), intOne(1), intUnit(-1,1);

	list<FactorTab> lst_unselected;
	list<FactorTab> lst_selected;

	Matrix paraTemplate(rangeDim, rangeDim);
	int num_selected = 0;

	vector<Interval> step_exp_table;
	Interval intStep;

	vector<vector<Interval> > step_exp_tables;

	for(int i=0; i<flowpipes.size(); ++i)
	{
		if(step_exp_table.size() == 0 || intStep != domains[i][0])
		{
			construct_step_exp_table(step_exp_table, domains[i][0], 2*globalMaxOrder);
			intStep = domains[i][0];
		}

		step_exp_tables.push_back(step_exp_table);
	}

	// 1: we first consider the user specified template vectors

	if(template_candidates.size() > 0)
	{
		for(int i=0; i<template_candidates.size(); ++i)
		{
			FactorTab facT(i, intZero, intInvalid);
			lst_unselected.push_back(facT);
			rhoPos.push_back(intInvalid);
			rhoNeg.push_back(intInvalid);
		}

		// 1.1: compute the intercepts

		step_exp_table = step_exp_tables[0];

		for(int i=0, k=0; i<flowpipes.size(); ++i)
		{
			if(step_exp_table[1] != domains[i][0])
			{
				step_exp_table = step_exp_tables[++k];
			}

			for(int j=0; j<template_candidates.size(); ++j)
			{
				Interval tmp1 = rhoNormal(flowpipes[i], template_candidates[j], step_exp_table);

				RowVector rowVec = template_candidates[j];
				rowVec.neg_assign();

				Interval tmp2 = rhoNormal(flowpipes[i], rowVec, step_exp_table);

				if(tmp1 >= rhoPos[j])	// tmp1.up > rhoPos[j].up
				{
					rhoPos[j] = tmp1;
				}

				if(tmp2 >= rhoNeg[j])	// tmp2.up > rhoNeg[j].up
				{
					rhoNeg[j] = tmp2;
				}
			}
		}

		// 1.2: select the vector with the smallest intercept
		list<FactorTab>::iterator facIter = lst_unselected.begin();
		list<FactorTab>::iterator min_facIter = facIter;

		facIter->intercept = rhoPos[0] + rhoNeg[0];
		++facIter;

		for(int i=1; facIter!=lst_unselected.end(); ++facIter, ++i)
		{
			facIter->intercept = rhoPos[i] + rhoNeg[i];

			if(facIter->intercept < min_facIter->intercept)
			{
				min_facIter = facIter;
			}
		}

		int index_selected = min_facIter->index;


/*
		// =============== test begin ==================
		printf("Candidates:\n");
		list<FactorTab>::iterator testIter = lst_unselected.begin();
		for(; testIter!=lst_unselected.end(); ++testIter)
		{
			printf("vector: ");
			template_candidates[testIter->index].dump(stdout);
			printf("\tintercept: %lf, factor: %lf\n\n", testIter->intercept, testIter->factor);
		}

		printf("selected: ");
		template_candidates[index_selected].dump(stdout);
		printf("\tintercept: %lf, factor: %lf\n\n\n", min_facIter->intercept, min_facIter->factor);
		// =============== test end ==================
*/



		lst_selected.push_back(*min_facIter);

		parallelotope_b.set(rhoPos[index_selected].sup(), num_selected);
		parallelotope_b.set(rhoNeg[index_selected].sup(), num_selected+rangeDim);

		lst_unselected.erase(min_facIter);

		for(int i=0; i<rangeDim; ++i)
		{
			paraTemplate.set(template_candidates[index_selected].get(i), num_selected, i);
		}

		++num_selected;

		if(num_selected < rangeDim)
		{
			// 1.3: update the factor table according to the selected vectors

			for(facIter=lst_unselected.begin(); facIter!=lst_unselected.end(); ++facIter)
			{
				facIter->factor.set( weightMat.get(index_selected, facIter->index) );
			}

			// 1.4: consider the remaining vectors

			for(; lst_unselected.size() != 0 && num_selected < rangeDim; )
			{
				FactorTab vector_selected;
				bool bselected = vector_selection(vector_selected, lst_unselected, paraTemplate, template_candidates, num_selected);

				// update the factors
				if(bselected)
				{
					parallelotope_b.set(rhoPos[vector_selected.index].sup(), num_selected-1);
					parallelotope_b.set(rhoNeg[vector_selected.index].sup(), num_selected-1+rangeDim);

					lst_selected.push_back(vector_selected);

					facIter = lst_unselected.begin();
					for(; facIter!=lst_unselected.end(); ++facIter)
					{
						facIter->factor.mul_assign( weightMat.get(facIter->index, vector_selected.index) );
					}
				}
				else
				{
					break;
				}
			}
		}
	}

	// 2: we consider the auto selected template vectors

	if(num_selected < rangeDim)
	{
/*
		// evaluate the template vector from the flowpipes
		int mid = flowpipes.size() / 2;
		vector<Interval> domainCenterPoint;
		for(int i=0; i<domains[mid].size(); ++i)
		{
			Interval I(domains[mid][i].midpoint());
			domainCenterPoint.push_back(I);
		}

		vector<Interval> flowCenterPoint;
		flowpipes[mid].intEval(flowCenterPoint, domainCenterPoint);

		Interval intZero;
		flowCenterPoint.insert(flowCenterPoint.begin(), intZero);

		bool bvalid_template_flow = false;
		Interval intSparsity(-SPARSITY, SPARSITY);

		RowVector template_flow(rangeDim);
		for(int i=0; i<odeHF.size(); ++i)
		{
			Interval I;
			odeHF[i].intEval(I, flowCenterPoint);
			template_flow.set(I.midpoint(), i);

			if(!bvalid_template_flow && !I.subseteq(intSparsity))
			{
				bvalid_template_flow = true;
			}
		}

		if(bvalid_template_flow)
		{
			template_flow.normalize();
		}

		// update the weight matrix

		if(bvalid_template_flow)
		{
			for(int i=0; i<weightMat.cols(); ++i)
			{
				if(i < auto_start)
				{
					// candidate vector
					double weight = 1 - fabs(template_candidates[i].innerProd(template_flow));
					weightMat.set(weight, i, flow_start);
					weightMat.set(weight, flow_start, i);
				}
				else if(i >= auto_start && i < flow_start)
				{
					// auto selected vector
					int posi = i - auto_start;
					if(linear_auto[startID][targetID][posi])
					{
						double weight = 1 - fabs(template_auto[startID][targetID][posi].innerProd(template_flow));
						weightMat.set(weight, i, flow_start);
						weightMat.set(weight, flow_start, i);
					}
					else
					{
						weightMat.set(-1, i, flow_start);
						weightMat.set(-1, flow_start, i);
					}
				}
				else if(i >= flow_start && i < default_start)
				{
					// flow vector
					double weight = 1 - fabs(template_flow.innerProd(template_flow));
					weightMat.set(weight, i, flow_start);
					weightMat.set(weight, flow_start, i);
				}
				else
				{
					// default vector
					double weight = 1 - fabs(template_default[i - default_start].innerProd(template_flow));
					weightMat.set(weight, i, flow_start);
					weightMat.set(weight, flow_start, i);
				}
			}
		}

		local_template_auto.push_back(template_flow);
*/


		rhoPos.clear();
		rhoNeg.clear();
		lst_unselected.clear();

		if(template_auto[startID][targetID].size() > 0)
		{
			for(int i=0; i<template_auto[startID][targetID].size(); ++i)
			{
				if(linear_auto[startID][targetID][i] && boundary_intersected[i])
				{
					FactorTab facT(i, intOne, intInvalid);
					lst_unselected.push_back(facT);
				}

				rhoPos.push_back(intInvalid);
				rhoNeg.push_back(intInvalid);
			}

			// 2.1: compute the support functions

			step_exp_table = step_exp_tables[0];

			for(int i=0, k=0; i<flowpipes.size(); ++i)
			{
				if(step_exp_table[1] != domains[i][0])
				{
					step_exp_table = step_exp_tables[++k];
				}

				list<FactorTab>::iterator vectorIter = lst_unselected.begin();
				for(; vectorIter!=lst_unselected.end(); ++vectorIter)
				{
					Interval tmp1 = rhoNormal(flowpipes[i], template_auto[startID][targetID][vectorIter->index], step_exp_table);

					RowVector rowVec = template_auto[startID][targetID][vectorIter->index];
					rowVec.neg_assign();

					Interval tmp2 = rhoNormal(flowpipes[i], rowVec, step_exp_table);

					if(tmp1 >= rhoPos[vectorIter->index])
					{
						rhoPos[vectorIter->index] = tmp1;
					}

					if(tmp2 >= rhoNeg[vectorIter->index])
					{
						rhoNeg[vectorIter->index] = tmp2;
					}
				}
			}

			// 2.2: compute the intercepts
			list<FactorTab>::iterator facIter = lst_unselected.begin();

			for(; facIter!=lst_unselected.end(); ++facIter)
			{
				facIter->intercept = rhoPos[facIter->index] + rhoNeg[facIter->index];
			}

			// 2.3: update the factor table
			list<FactorTab>::iterator selectedIter = lst_selected.begin();
			for(; selectedIter!=lst_selected.end(); ++selectedIter)
			{
				for(facIter=lst_unselected.begin(); facIter!=lst_unselected.end(); ++facIter)
				{
					facIter->factor.mul_assign( weightMat.get(selectedIter->index, facIter->index+auto_start) );
				}
			}

			// 2.4: consider the remaining vectors

			for(; lst_unselected.size() != 0 && num_selected < rangeDim; )
			{
				FactorTab vector_selected;
				bool bselected = vector_selection(vector_selected, lst_unselected, paraTemplate, template_auto[startID][targetID], num_selected);

				// update the factors
				if(bselected)
				{
					parallelotope_b.set(rhoPos[vector_selected.index].sup(), num_selected-1);
					parallelotope_b.set(rhoNeg[vector_selected.index].sup(), num_selected-1+rangeDim);

					vector_selected.index += auto_start;
					lst_selected.push_back(vector_selected);

					facIter = lst_unselected.begin();
					for(; facIter!=lst_unselected.end(); ++facIter)
					{
						facIter->factor.mul_assign( weightMat.get(facIter->index+auto_start, vector_selected.index) );
					}
				}
				else
				{
					break;
				}
			}
		}
	}

	// 3: if the vectors are not enough, we find the remaining ones in the default template

	if(num_selected < rangeDim)
	{
		rhoPos.clear();
		rhoNeg.clear();
		lst_unselected.clear();

		if(template_default.size() > 0)
		{
			for(int i=0; i<template_default.size(); ++i)
			{
				FactorTab facT(i, intOne, intInvalid);
				lst_unselected.push_back(facT);
				rhoPos.push_back(intInvalid);
				rhoNeg.push_back(intInvalid);
			}

			// 3.1: compute the support functions

			step_exp_table = step_exp_tables[0];

			for(int i=0, k=0; i<flowpipes.size(); ++i)
			{
				if(step_exp_table[1] != domains[i][0])
				{
					step_exp_table = step_exp_tables[++k];
				}

				for(int j=0; j<template_default.size(); ++j)
				{
					Interval tmp1 = rhoNormal(flowpipes[i], template_default[j], step_exp_table);

					RowVector rowVec = template_default[j];
					rowVec.neg_assign();

					Interval tmp2 = rhoNormal(flowpipes[i], rowVec, step_exp_table);

					if(tmp1 >= rhoPos[j])
					{
						rhoPos[j] = tmp1;
					}

					if(tmp2 >= rhoNeg[j])
					{
						rhoNeg[j] = tmp2;
					}
				}
			}

			// 3.2: compute the intercepts
			list<FactorTab>::iterator facIter = lst_unselected.begin();

			for(int i=0; facIter!=lst_unselected.end(); ++facIter, ++i)
			{
				facIter->intercept = rhoPos[i] + rhoNeg[i];
			}

			// 3.3: update the factor table
			list<FactorTab>::iterator selectedIter = lst_selected.begin();
			for(; selectedIter!=lst_selected.end(); ++selectedIter)
			{
				for(facIter=lst_unselected.begin(); facIter!=lst_unselected.end(); ++facIter)
				{
					facIter->factor.mul_assign( weightMat.get(selectedIter->index, facIter->index+default_start) );
				}
			}

			// 3.4: consider the remaining vectors
			for(; lst_unselected.size() != 0 && num_selected < rangeDim; )
			{
				FactorTab vector_selected;
				bool bselected = vector_selection(vector_selected, lst_unselected, paraTemplate, template_default, num_selected);

				// update the factors
				if(bselected)
				{
					parallelotope_b.set(rhoPos[vector_selected.index].sup(), num_selected-1);
					parallelotope_b.set(rhoNeg[vector_selected.index].sup(), num_selected-1+rangeDim);

					vector_selected.index += default_start;
					lst_selected.push_back(vector_selected);

					facIter = lst_unselected.begin();
					for(; facIter!=lst_unselected.end(); ++facIter)
					{
						facIter->factor.mul_assign( weightMat.get(facIter->index+default_start, vector_selected.index) );
					}
				}
				else
				{
					break;
				}
			}
		}
	}

	// 4: we use the template parallelotope to over-approximate the flowpipe union

	Parallelotope paraAggregation(paraTemplate, parallelotope_b);

	// converse the parallelotope to a Taylor model and construct a normalized domain for it
	paraAggregation.toTaylorModel(tmvAggregation);

	doAggregation.clear();

	doAggregation.push_back(intZero);

	for(int i=0; i<rangeDim; ++i)
	{
		doAggregation.push_back(intUnit);
	}
}






