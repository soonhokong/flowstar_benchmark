FLOWSTAR_BIN:=../../flowstar-1.2.0/flowstar
GNUPLOT_BIN:=gnuplot
EPSPDF_BIN=epstopdf

all: water_double_unsat water_double_sat

water_double_unsat:
	$(FLOWSTAR_BIN) < $@.model
	if [ -e outputs/$@.plt ] ; then $(GNUPLOT_BIN) < outputs/$@.plt; fi
	if [ -e images/$@.eps ] ; then $(EPSPDF_BIN) images/$@.eps; fi

water_double_sat:
	$(FLOWSTAR_BIN) < $@.model
	if [ -e outputs/$@.plt ] ; then $(GNUPLOT_BIN) < outputs/$@.plt; fi
	if [ -e images/$@.eps ] ; then $(EPSPDF_BIN) images/$@.eps; fi
