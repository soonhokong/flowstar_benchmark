/*---
  Flow*: A Taylor Model Based Flowpipe analyzer.
  Authors: Xin Chen, Erika Abraham and Sriram Sankaranarayanan.
  Email: Xin Chen <xin.chen@cs.rwth-aachen.de> if you have questions or comments.
  
  The code is released as is under the GNU General Public License (GPL). Please consult the file LICENSE.txt for
  further information.
---*/

#ifndef MODELPARSER_H_
#define MODELPARSER_H_

#include "Hybrid.h"

extern int lineNum;

extern mpfr_prec_t intervalNumPrecision;

extern ContinuousReachability continuousProblem;
extern HybridReachability hybridProblem;

extern ParseSetting parseSetting;
extern ParseResult parseResult;

extern int yyparse();

extern vector<Interval> gUncertainties;

void parseError(const char *str, int lnum);

#endif /* MODELPARSER_H_ */
