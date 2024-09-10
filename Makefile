
FBJS_calc_carbchem_user_simple = \
                    ../ctoaster.carrotcake/src/common/gem_cmn.f90 \
                    ../ctoaster.carrotcake/src/common/gem_util.f90 \
                    ../ctoaster.carrotcake/src/common/gem_carbchem.f90 \
                    calc_carbchem_user_simple.f90

#### fortran compiler command name
FC = gfortran

#### COMPILER FLAGS
FLAGS = -fdefault-real-8 -O2 -O3 -funroll-loops -msse
#FLAGS += -g -ffpe-trap=zero,overflow,invalid -O0 -Wall -fbounds-check

calc_carbchem_user: 
	$(FC) $(FLAGS) $(FBJS_calc_carbchem_user_simple) -o testcarrot.exe
clean: 
	rm -f *.o *.mod *.exe
