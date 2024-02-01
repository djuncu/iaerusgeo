# Version production?
LIVR = 0

# 32 ou 64 bits?
PROC = $(shell uname -i)
ifeq ($(PROC), x86_64)
	ENV = 64
else
        ENV = 32
endif

# RedHat Release?
REDHAT = $(shell /usr/bin/lsb_release -sr)
RH_VER = $(word 1, $(subst ., ,$(REDHAT)))

# Options d'optimisation
ifneq ($(LIVR),0)
        COPTLVL=-O2
        FOPTLVL=-O2
else
        COPTLVL=-g
        FOPTLVL=-ggdb -fbounds-check
endif

S1 = src
S = C_F90
NR = $(S)/nr

FC = gfortran
CC = gcc
F2PY = f2py
#3
PYTHON_SO_SUFFIX := $(shell python3-config --extension-suffix)

SRC =	$(NR)/nrtype.f90		\
	$(NR)/nrutil.f90		\
	$(NR)/nr.f90			\
	$(NR)/ludcmp.f90		\
	$(NR)/lubksb.f90		\
	$(S)/py_ifc.f90			\
	$(S)/algoconf.f90		\
	$(S)/solzenith.f90		\
	$(S)/julday.f90			\
	$(S)/brdfmodels.f90		\
	$(S)/invers.f90			\
	$(S)/LM.f90			\
	$(S)/flotsam_interface.f90	\
	$(S)/calculate_ref_tol.f90	\
	$(S)/start_mod.f90		\
	$(S)/aerosol_parameters.f90	\
	$(S)/start_exe.f90	        \
	$(S)/smoothing.f90


TMPOBJ = $(SRC:.c=.o)
objects = $(TMPOBJ:.f90=.o)

# Chemins des binaires, includes et bibliotheques installees
MAIN_ENV_BASE = /usr
#MAIN_ENV_BASE = /home/six
#MAIN_ENV_DIR = $(MAIN_ENV_BASE)/ops/env$(ENV)_rhel$(RH_VER)
MAIN_ENV_DIR = $(MAIN_ENV_BASE)

MOD = AerusL3
LIB = lib$(MOD).so

CFLAGS = $(COPTLVL) -std=c99 -fPIC -m$(ENV)
FFLAGS = $(FOPTLVL) -funroll-loops -I$(MAIN_ENV_DIR)/include -I$(MAIN_ENV_DIR)/lib -ffree-line-length-none -fPIC -m$(ENV) -fcheck=all -fmax-array-constructor=137000 -fallow-invalid-boz -fallow-argument-mismatch
# double check if -L$(S1) is necessary
F2PYFLAGS = -L$(MAIN_ENV_DIR)/lib -L. -l$(MOD) -L$(S1)

FLOTSAM_DIR = /cnrm/vegeo/juncud/software/flotsam-compiled
FLOTSAM_LIB_FLAGS = -I$(FLOTSAM_DIR)/include -L$(FLOTSAM_DIR)/lib -Wl,-rpath,$(FLOTSAM_DIR)/lib -lflotsam
F_FLOTSAM_FLAGS = -O2 -funroll-loops -fstrength-reduce -frerun-loop-opt -fcheck=all -shared -o

# the compilation is certainly ﬁckle....
# first, NPY_DISTUTILS_APPEND_FLAGS=1 was required, now it breaks the compilation. 
# (actually since numpy 1.18 "1" should be default anyway)
# Not sure why, but have to keep an eye on it.
#unexport LD_LIBRARY_FLAGS

# double check if export statement is necessary
export LDFLAGS="-Wl,-rpath,$(S1) -Wl,-rpath,."
#export -e NPY_DISTUTILS_APPEND_FLAGS=1
#export FC=gfortran # << probably unnecessary

$(MOD): $(LIB)
	# NOTE juncud:24Nov2023: unexpectedly this starts to crash without any obvious changes to the last time I tried compiling.
	#   the issue seems to be with calling f2py with the -c flag.
	#   might have something to do with this warning? "[...]python3.11/site-packages/numpy/f2py/f2py2e.py:686: VisibleDeprecationWarning:
	#   distutils has been deprecated since NumPy 1.26.Use the Meson backend instead, or generate wrappers without -c and use a custom build script"

	#$(F2PY) -c -m $@ --fcompiler=gnu95 $(S)/py_ifc.f90 $(S)/start_mod.f90 $(S)/start_exe_flotsam.f90 $(S)/smoothing.f90 --f90flags="$(FFLAGS) $(FLOTSAM_LIB_FLAGS)" $(F2PYFLAGS) $(F2PY_FLOTSAM_FLAGS)
	
	# this needs to be commented if using the makef2py script....
	#f2py -c -m $@ $(S)/py_ifc.f90 $(S)/start_mod.f90 $(S)/start_exe_flotsam.f90 $(S)/start_exe.f90 $(S)/smoothing.f90 --f90flags="$(FFLAGS)" $(F2PYFLAGS)
	#f2py -m $@ $(S)/py_ifc.f90 $(S)/start_mod.f90 $(S)/start_exe_flotsam.f90 $(S)/start_exe.f90 $(S)/smoothing.f90 #$(F2PYFLAGS)

$(LIB):
	$(FC) -cpp -fPIC -shared -o $@ $(SRC) $(FFLAGS) $(FLOTSAM_LIB_FLAGS)

all   :
	make clear
	make

clean:
	rm -f $(LIB) *.mod *~ */*~ */*/*~

clear:
	make clean
	rm -f *.a *.pyc */*.pyc */*/*.pyc *.so $(MOD)$(PYTHON_SO_SUFFIX)
