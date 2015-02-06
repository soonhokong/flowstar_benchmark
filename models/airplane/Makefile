FLOWSTAR_BIN:=../../flowstar-1.2.0/flowstar
GNUPLOT_BIN:=gnuplot
EPSPDF_BIN=epstopdf

all: airplane_single_unsat airplane_single_sat

airplane_single_unsat:
	$(FLOWSTAR_BIN) < $@.model
	if [ -e outputs/$@.plt ] ; then $(GNUPLOT_BIN) < outputs/$@.plt; fi
	if [ -e images/$@.eps ] ; then $(EPSPDF_BIN) images/$@.eps; fi

airplane_single_sat:
	$(FLOWSTAR_BIN) < $@.model
	if [ -e outputs/$@.plt ] ; then $(GNUPLOT_BIN) < outputs/$@.plt; fi
	if [ -e images/$@.eps ] ; then $(EPSPDF_BIN) images/$@.eps; fi

airplane_single_nl_sat:
	$(FLOWSTAR_BIN) < $@.model
	if [ -e outputs/$@.plt ] ; then $(GNUPLOT_BIN) < outputs/$@.plt; fi
	if [ -e images/$@.eps ] ; then $(EPSPDF_BIN) images/$@.eps; fi

nonlinear_test:
	$(FLOWSTAR_BIN) < $@.model
	if [ -e outputs/$@.plt ] ; then $(GNUPLOT_BIN) < outputs/$@.plt; fi
	if [ -e images/$@.eps ] ; then $(EPSPDF_BIN) images/$@.eps; fi
