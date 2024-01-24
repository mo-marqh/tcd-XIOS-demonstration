FCFLAGS = -g -ffree-line-length-none
NETCDF_LIB_DIR64 = $(NETCDF_DIR)/lib64

FC = mpif90       # compiler driver for MPI programs


# compiler flag, includes
FCFLAGS += -I$(XIOS_INCDIR) $(NETCDF_INCDIR)

# loader flags
LDFLAGS = \
	-L$(XIOS_LIBDIR) \
	$(NETCDF_LIBDIR) \
	-lxios \
	-lnetcdf \
	-lnetcdff \
	-lstdc++ 

all: nc resample

# fortran compilation
%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

# fortran linking
resample: resample.o
	$(FC) -o resample.exe resample.o $(LDFLAGS) \
		&& ln -fs $(XIOS_BINDIR)/xios_server.exe .

nc:
	ncgen -k nc4 -o axis_input.nc axis_input.cdl

run:
	mpiexec -n 1 ./resample.exe : -n 1 ./xios_server.exe

# cleanup
clean:
	rm -f *.exe *.o *.mod *.MOD *.out *.err *.nc
