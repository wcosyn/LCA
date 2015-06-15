# A Makefile with auto dep generation
# Note that the build will fail if header
# files are renamed. This is because the
# *.d dependency files will be out of date
# and point to non-existing header files
# We can solve this automatically but to
# keep the code nice and simple just do
# $> make clean (or rm *.d)
# $> make
# If header files have been renamed
#
# @author Camille Colle

CC=g++
WFLAGS=-Wall -Wextra #-Wshadow #-Weffc++ # You can expect to see a buttload of warnings when enabling all these flags
CFLAGS= -O3 -std=c++11 -fPIC $(WFLAGS) `gsl-config --cflags` -fopenmp
# -fext-numeric-literals is there to allow conversion between C and C++ complex numbers
#
# -Wl means pass arguments to runtime linker. -Wl is followed by comma seperated list with arguments to be passed to runtime linker
#  This is needed because the -L, -l flags only link @compile time but not runtime! If the libraries are not in the standard library paths
#  or LD_LIBRARY_PATH is not set we will get in trouble because the dynamic linker will (ld) will not be able to find our libs @runtime. 
#  To avoid using LD_LIBRARY_PATH we use rpath...
LIBS=-L../lib -Wl,-rpath=../lib `gsl-config --libs` 
INC_DIR=../include
INC=-I$(INC_DIR)

TARGETS=mainwswf mainop
SRC_DIR=../src
# object files for building, no main()'s here !!!
SRC:=$(wildcard $(SRC_DIR)/*.cpp)
OBJS:=$(notdir $(filter-out $(TARGETS:=.o), $(SRC:.cpp=.o))) # all cpp files in current directory
DEPDIR:=dep
DEPS=$(patsubst %.o, $(DEPDIR)/%.d, $(OBJS)) # make dep files


.DEFAULT: all
all: $(TARGETS) $(SOLIBS)


-include $(DEPS)


# if object_deps gets called from rule with target $@ will contain the target's name
define object_deps
	@tput setf 3; tput bold
	@printf "\n\nobject_deps called with arguments $@\n\n";
	@cat $(patsubs %,%.cpp,$@) | grep -E "\<.*\>"
	@tput sgr0;
endef

$(SOLIBS): %.so: %.o
	@tput setaf 4; tput bold
	@printf "Building shared object $@ : $^ \n"
	@tput sgr0
	$(CC) -shared $< $(CFLAGS) $(INC) $(LIBS) -o $@
	@printf "\n"

# dependency of each target is all the object files and it's own %.o file.
# Technically too much dependencies... (targets will generally not depend on all $(OBJS) files)
# But it's better to have too much dependencies then too few.
$(TARGETS): %: $(SRC_DIR)/%.o $(OBJS)
	@printf "wildcard test: $(SRC:.cpp=.o)\n"
	@printf "\n\nOBJS VARIABLE IS NOW:\n$(OBJS)\n\n"
	@printf "\npattern:\n$(TARGETS:=.o)\n\n"
	@printf "filter test:\n$(filter-out $(TARGETS).o, $(OBJS))\n\n"
	@tput setaf 1; tput bold
	@printf "Linking executable $@  :  $^ \n"
	@tput sgr0
	$(CC) $^ $(CFLAGS) $(INC) $(LIBS) -o $@
	@printf "\n"
	@mkdir -p ../bin
	mv $@ ../bin/


%.o : %.cpp
	@echo DEPS INCLUDED ARE $(DEPS)
	@mkdir -p $(DEPDIR)
	@tput setf 6
	@printf "Building dependencies for $<\n"
	@tput sgr0
	$(CC) -MM $< $(CFLAGS) $(INC) > $(patsubst %.o,$(DEPDIR)/%.d, $(notdir $@))
	@tput setf 6
	@printf "Building  $@ :  $^ \n"
	@tput sgr0
	$(CC) -c $< $(CFLAGS) $(INC) -o $@

.PHONY: clean

clean:
	rm *.o $(DEPDIR)/*.d 
