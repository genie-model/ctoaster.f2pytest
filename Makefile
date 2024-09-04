
FBJS_calc_carbchem = \
                    ../../cgenie.muffin/genie-main/src/fortran/cmngem/gem_cmn.f90 \
                    ../../cgenie.muffin/genie-main/src/fortran/cmngem/gem_util.f90 \
                    ../../cgenie.muffin/genie-main/src/fortran/cmngem/gem_carbchem.f90 \
                    calc_carbchem.f90

#### fortran compiler command name
FC = gfortran

#### COMPILER FLAGS
FLAGS = -fdefault-real-8 -O2 -O3 -funroll-loops -msse
#FLAGS += -g -ffpe-trap=zero,overflow,invalid -O0 -Wall -fbounds-check

calc_carbchem: 
	$(FC) $(FLAGS) $(FBJS_calc_carbchem) -o calc_carbchem.exe
clean: 
	rm -f *.o *.mod *.exe
