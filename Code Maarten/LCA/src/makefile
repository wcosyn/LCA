HDRS = newcoef.h\
	   recmosh.h\
	   pair.h\
	   nucleus.h\
	   nucleuspp.h\
	   nucleusnn.h\
	   nucleusnp.h\
	   shell.h\
	   wavefunctionp.h\
	   operator_virtual.h\
	   density_rel.h\
	   rms_ob.h\
	   norm_tb.h\
	   nucleusall.h\
	   norm_ob.h\
	   density_ob3.h\
	   operator_virtual_ob.h\
	   density_ob_integrand3.h\
	   density_ob_integrand_cf.h\
	   correlation_functions.h\
	   triplet.h\
	   threebodycoef.h\
	   density_rel_integrand2.h\
	   speedy.h\
	   kinenergy_ob.h\
	   paircoef.h\
	   tripletcoef.h\
	   threej.h\
	   density_td.h\
	   wsnucleus.h\
	   wsexpansion.h\
	   wsnucleusnp.h\
	   wsnucleuspp.h\
	   wspair.h\
	   wswf.h\

	   
	   
OBJS =     newcoef.o\
	   recmosh.o\
	   pair.o\
	   nucleus.o\
	   nucleuspp.o\
	   nucleusnn.o\
	   nucleusnp.o\
	   shell.o\
	   wavefunctionp.o\
	   operator_virtual.o\
	   density_rel.o\
	   rms_ob.o\
	   norm_tb.o\
	   nucleusall.o\
	   norm_ob.o\
	   density_ob3.o\
	   operator_virtual_ob.o\
	   density_ob_integrand3.o\
	   density_ob_integrand_cf.o\
	   correlation_functions.o\
	   triplet.o\
	   threebodycoef.o\
	   density_rel_integrand2.o\
	   speedy.o\
	   kinenergy_ob.o\
	   paircoef.o\
	   tripletcoef.o\
	   threej.o\
	   density_td.o\
	   wsnucleus.o\
	   wsexpansion.o\
	   wsnucleusnp.o\
	   wsnucleuspp.o\
	   wspair.o\
	   wswf.o\


#MAINWF = maintest.o\

MAINWS = mainwswf.o\

MAINDENS = maindens.o\

MAINOP = mainop.o\

RFLAGS=-Wl,-rpath,'$$ORIGIN/lib'

LFLAGS=`gsl-config --libs` -lpcuba -lhcuba -L../lib -fopenmp

CFLAGS= `gsl-config --cflags` -I../include -fopenmp

OFLAGS= -std=c++11 -O3
	   
EFILE =	../ws\
	../op


all: ws op
clean:; @rm -f *.o
wf: ../wf
ws: ../ws
dens: ../dens
op: ../op

#../wf : $(HDRS) $(OBJS) $(MAINWF)
#	@ echo "Linking wf ..."
#	g++ -Wall ${OFLAGS}  -o ../wf $(OBJS) $(MAINWF) ${LFLAGS} ${RFLAGS}

../ws : $(HDRS) $(OBJS) $(MAINWS)
	@ echo "Linking ws ..."
	g++ -Wall ${OFLAGS} -o ../ws $(OBJS) $(MAINWS) ${LFLAGS} ${RFLAGS}

../dens : $(HDRS) $(OBJS) $(MAINDENS)
	@ echo "Linking dens ..."
	g++ -Wall ${OFLAGS} -o ../dens $(OBJS) $(MAINDENS) ${LFLAGS} ${RFLAGS}

../op 	: $(HDRS) $(OBJS) $(MAINOP)
	@ echo "Linking dens ..."
	g++ -Wall ${OFLAGS} -o ../op $(OBJS) $(MAINOP) ${LFLAGS} ${RFLAGS}
	
$(OBJS) : $(HDRS)

#$(MAINWF) : $(HDRS)

$(MAINWS) : $(HDRS)

$(MAINDENS) : $(HDRS)

%.o	: %.cpp
	@ echo "Compiling..."
	g++ -Wall ${OFLAGS} -c $*.cpp  ${CFLAGS}
	
# $(EFILE) : $(OBJS) $(HDRS)
# 	@ echo "Linking $(EFILE)..."
# 	g++ -Wall -o $(EFILE) $(OBJS)  -lgsl -lm 
# 
# $(OBJS) : $(HDRS)
# NOTES
# Implicit rule for .o files: 
# $(CXX) $(CPPFLAGS) $(CXXFLAGS) -c
#
# One file can be the target of several rules.
# All the prerequisites mentioned in all the rules are merged into one list of prerequisites for the target.
# If the target is older than any prerequisite from any rule, the recipe is executed. 
# There can only be one recipe to be executed for a file.
# If more than one rule gives a recipe for the same file, make uses the last one given and prints an error message.
