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

#S = src
S = C_F90
NR = $(S)/nr

FC = gfortran
CC = gcc

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
	$(S)/calculate_ref_tol.f90	\
	$(S)/start_mod.f90		\
	$(S)/aerosol_parameters.f90	\
	$(S)/start_exe.f90		\
	$(S)/smoothing.f90
TMPOBJ = $(SRC:.c=.o)
objects = $(TMPOBJ:.f90=.o)

# Chemins des binaires, includes et bibliotheques installees
MAIN_ENV_BASE = /usr
#MAIN_ENV_BASE = /home/six
#MAIN_ENV_DIR = $(MAIN_ENV_BASE)/ops/env$(ENV)_rhel$(RH_VER)
MAIN_ENV_DIR = $(MAIN_ENV_BASE)

MOD = AerusL3
LIB = lib$(MOD).a

CFLAGS = $(COPTLVL) -std=c99 -fPIC -m$(ENV)
FFLAGS = $(FOPTLVL) -funroll-loops -I$(MAIN_ENV_DIR)/include -I$(MAIN_ENV_DIR)/lib -ffree-line-length-none -fPIC -m$(ENV) -fcheck=all -fmax-array-constructor=137000 -fallow-invalid-boz -fallow-argument-mismatch
F2PYFLAGS = -L$(MAIN_ENV_DIR)/lib -L. -l$(MOD)

$(MOD): $(LIB)
	f2py -c -m $@ --fcompiler=gnu95 --f90flags="$(FFLAGS)" $(S)/py_ifc.f90 $(S)/start_mod.f90 $(S)/start_exe.f90 $(S)/smoothing.f90 $(F2PYFLAGS)

$(LIB): $(objects)
	ar -cr $@ $(objects) *.mod

%.o : %.f90
	$(FC) -c $< -o $@ $(FFLAGS)

all   :
	make clear
	make

clean:  
	rm -f $(LIB) *.mod $(objects) *~ */*~ */*/*~

clear:
	make clean
	rm -f *.a *.pyc */*.pyc */*/*.pyc $(MOD).so
