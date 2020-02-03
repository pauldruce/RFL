################################################################################
#
# Makefile to compile and link C++ programs
#
# Version valid for Linux machines
#
# "make" compiles and links the specified main programs and modules
# using the specified libraries (if any), and produces the executables
# 
# "make clean" removes all files created by "make"
#
# Dependencies on included files are automatically taken care of
#
################################################################################

all: rmxeq mkdep mkxeq
.PHONY: all


# main programs and required modules 

MAIN = duav main

SOURCE = action clifford derivative hamiltonian metropolis HMC_core HMC misc MMC_core MMC leapfrog

# search path for modules

VPATH = ../../source


# additional include directories

INCPATH = ../../include /home/pmxmd10/gsl/include


# additional libraries to be included 
 
LIBS = gsl openblas armadillo

LIBPATH = /home/pmxmd10/gsl/lib


# scheduling and optimization options (such as -DSSE -DSSE2 -DP4)
 
CFLAGS = -std=c++11 -pedantic -fstrict-aliasing -Wall -Wno-long-long -Werror 

############################## do not change ###################################

SHELL=/bin/bash

CC=g++

PGMS= $(MAIN) $(SOURCE)

INCDIRS = $(addprefix -I,$(INCPATH))

OBJECTS = $(addsuffix .o,$(SOURCE))

LDFLAGS = $(addprefix -L,$(LIBPATH)) $(addprefix -l,$(LIBS))

-include $(addsuffix .d,$(PGMS))


# rule to make dependencies

$(addsuffix .d,$(PGMS)): %.d: %.cpp Makefile
	@ $(CC) -MM $(INCDIRS) $< -o $@


# rule to compile source programs

$(addsuffix .o,$(PGMS)): %.o: %.cpp Makefile
	$(CC) $< -c $(CFLAGS) $(INCDIRS) -o $@


# rule to link object files

$(MAIN): %: %.o $(OBJECTS) Makefile
	$(CC) $< $(OBJECTS) $(CFLAGS) $(LDFLAGS) -o $@


# produce executables

mkxeq: $(MAIN)


# remove old executables and old error log file

rmxeq:
	@ -rm -f $(MAIN); \
        echo "delete old executables"		


# make dependencies

mkdep:  $(addsuffix .d,$(PGMS))
	@ echo "generate tables of dependencies"


# clean directory 

clean:
	@ -rm -rf *.d *.o .tmp $(MAIN)
.PHONY: clean

################################################################################
