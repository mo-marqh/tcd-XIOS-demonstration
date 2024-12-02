# Make file for the resample demonstartion XIOS programme
# Targets provided our detailed below...
#
# all: (default) Build the context programme
# clean: Delete all final products and working files
# run: run the programme
#
# Environment Variables expected by this MakeFile:
#
# FC: mpif90
# FCFLAGS: -g & include files for netcdf & xios
# LD_FLAGS: for xios, netcdf, netcdff, stfc++
# LD_LIBRARY_PATH: for netCDF & XIOS libs
# XIOS_BINDIR: The directory for XIOS binary files

.PHONY: all, clean, run

OBJ = custom_type_mod.o context_def_mod.o

all: context_def_test

custom_type_mod.o: custom_type_mod.F90
	$(FC) $(FCFLAGS) -c $<

context_def_mod.o: context_def_mod.F90 custom_type_mod.o
	$(FC) $(FCFLAGS) -c $<

context_def_test.o: context_def_test.F90 custom_type_mod.o context_def_mod.o
	$(FC) $(FCFLAGS) -c $<

# fortran linking
context_def_test: context_def_test.o $(OBJ)
	$(FC) -o context_def_test.exe context_def_test.o $(OBJ) $(LDFLAGS) \
		&& ln -fs $(XIOS_BINDIR)/xios_server.exe .

run:
	mpiexec -n 1 ./context_def_test.exe : -n 1 ./xios_server.exe

clean:
	rm -f *.exe *.o *.mod *.MOD *.out *.err *.nc
