%{
/*---
  Flow*: A Taylor Model Based Flowpipe analyzer.
  Authors: Xin Chen, Erika Abraham and Sriram Sankaranarayanan.
  Email: Xin Chen <xin.chen@cs.rwth-aachen.de> if you have questions or comments.
  
  The code is released as is under the GNU General Public License (GPL). Please consult the file LICENSE.txt for
  further information.
---*/

#include "modelParser.h"
#include "modelParser.tab.h"
%}

delim [ \t\r]
line [\n]
whitespace {delim}+
digit [0-9]
letter [a-zA-Z_]
ident {letter}({letter}|{digit}|".")*
number ("-"?)(({digit}+)|({digit}*"."{digit}*)|({digit}+"e"(({digit}+)|("-"{digit}+)))|({digit}*"."{digit}*"e"(({digit}+)|("-"{digit}+)))|({digit}*"."{digit}*"e"("+"({digit}+)|("-"{digit}+))))


%%
"\'" { return '\'';}
"+" {return '+';}
"*" {return '*';}
"-" {return '-';}
"=" {return EQ;}
">=" {return GEQ;}
"<=" {return LEQ;}
"," {return ',';}
";" {return ';';}
":" {return ':';}
"(" {return '(';}
")" {return ')';}
"{" {return '{';}
"}" {return '}';}
"[" {return '[';}
"]" {return ']';}
":=" {return ASSIGN;}
"^" {return '^';}
"/" {return '/';}
">" {return '>';}
"mode" {return MODE;}
"init" {return INIT;}
"in" {return BELONGSTO;}
"poly ode 1" {return POLYODE1;}
"poly ode 2" {return POLYODE2;}
"visualize" {return VISUALIZE;}
"end" {return END;}
"setting" {return SETTING;}
"continuous reachability" {return CONTINUOUS;}
"hybrid reachability" {return HYBRID;}
"fixed steps" {return FIXEDST;}
"adaptive steps" {return ADAPTIVEST;}
"fixed orders" {return FIXEDORD;}
"adaptive orders" {return ADAPTIVEORD;}
"min" {return MIN;}
"max" {return MAX;}
"remainder estimation" {return REMEST;}
"interval" {return INTERVAL;}
"octagon" {return OCTAGON;}
"grid" {return GRID;}
"QR precondition" {return QRPRECOND;}
"identity precondition" {return IDPRECOND;}
"time" {return TIME;}
"modes" {return MODES;}
"jumps" {return JUMPS;}
"inv" {return INV;}
"guard" {return GUARD;}
"reset" {return RESET;}
"start" {return START;}
"max jumps" {return MAXJMPS;}
"output" {return OUTPUT;}
"print on" {return PRINTON;}
"print off" {return PRINTOFF;}
"unsafe set" {return UNSAFESET;}
"state var" {return STATEVAR;}
"tm var" {return TMVAR;}
"taylor model" {return TM;}
"parallelotope aggregation" {return PARAAGGREG;}
"interval aggregation" {return INTAGGREG;}
"taylor model aggregation" {return TMAGGREG;}
"continuous flowpipes" {return CONTINUOUSFLOW;}
"hybrid flowpipes" {return HYBRIDFLOW;}
"taylor picard" {return TAYLOR_PICARD;}
"taylor remainder" {return TAYLOR_REMAINDER;}
"taylor polynomial" {return TAYLOR_POLYNOMIAL;}
"exp" {return EXP;}
"sin" {return SIN;}
"cos" {return COS;}
"log" {return LOG;}
"sqrt" {return SQRT;}
"nonpoly ode" {return NPODE_TAYLOR;}
"cutoff" {return CUTOFF;}
"precision" {return PRECISION;}
"gnuplot" {return GNUPLOT;}
"matlab" {return MATLAB;}
"computation paths" {return COMPUTATIONPATHS;}


{number} { yylval.dblVal = atof( (char *)yytext ); return NUM; }

{ident}	{ yylval.identifier = new string(yytext); return IDENT; }

{whitespace} 

{line} { lineNum++; }

"#" {	/* Comment line */
	int c;
	c = yyinput();
	while(c!='\n' && c!=EOF)
	{
		c = yyinput();
	}
	
	if(c == '\n')
	{
		++lineNum;
	}
}

.  { printf("Unknown Character in line %d : %s -- Ignored\n", lineNum, yytext); }

%%


int yywrap()
{
	return 1;
}

void parseODE()
{
	int length = parseSetting.strODE.size();
	char *strODE = new char[length+5];

	sprintf(strODE, "%s", parseSetting.strODE.c_str());

	YY_BUFFER_STATE current = YY_CURRENT_BUFFER;

	YY_BUFFER_STATE yybs = yy_scan_string(strODE);

	yy_switch_to_buffer(yybs);
    yyparse();
    yy_delete_buffer(yybs);

	yy_switch_to_buffer(current);

	delete [] strODE;
}