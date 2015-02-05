FLOWSTAR_BIN:=../../flowstar-1.2.0/flowstar 
GNUPLOT_BIN:=gnuplot
EPSPDF_BIN=epstopdf

all: thermostat_double_unsat thermostat_double_sat

thermostat_double_unsat:
	$(FLOWSTAR_BIN) < $@.model
	if [ -e outputs/$@.plt ] ; then $(GNUPLOT_BIN) < outputs/$@.plt; fi
	if [ -e images/$@.eps ] ; then $(EPSPDF_BIN) images/$@.eps; fi

thermostat_double_sat:
	$(FLOWSTAR_BIN) < $@.model
	if [ -e outputs/$@.plt ] ; then $(GNUPLOT_BIN) < outputs/$@.plt; fi
	if [ -e images/$@.eps ] ; then $(EPSPDF_BIN) images/$@.eps; fi

