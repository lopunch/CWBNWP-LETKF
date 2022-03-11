SHELL = 	/bin/sh

CMD     = cwb_letkf.exe
CPP     = cpp -C -P -traditional -Wno-invalid-pp-token -ffreestanding

#FC      = mpifrtpx -X08 -Nquickdbg 
#FC      = mpifrtpx -X08 -O2 -Kautoobjstack 
 FC      = mpifrtpx -X08 -Kfast
CPPFLAGS = -DREAL64 #-DNC4
FCFLAGS  =
 LDFLAGS = -SSL2 -I/package/fx1000/netcdf-4.7.4/include \
                 -L/package/fx1000/netcdf-4.7.4/lib -lnetcdff -lnetcdf \
                 -L/package/fx1000/hdf5-1.10.7/lib -lhdf5 -lhdf5_hl -lhdf5 -lz

#-----------------------------------------------------------------------------

RM = rm -f
AR = ar -r

####MODULES = config.mod \
####		  param.mod  \
####		  grid.mod \
####		  netcdf_io.mod \
####		  gts_omboma.mod \
####		  simulated_radar.mod \
####		  kdtree.mod \
####		  eigen.mod \
####		  mpi_util.mod


LIBOBJ  = module_config.o \
          module_param.o  \
		  module_grid.o \
          module_netcdf_io.o \
		  module_gts_omboma.o \
		  module_radar.o \
		  module_kdtree2.o \
		  module_eigen.o \
		  module_mpi_util.o \
		  module_letkf_core.o \
		  module_localization.o \
		  module_projection.o

.SUFFIXES: .f90 .o .mod

all: $(CMD)

#-----------------------------------------------------------------------------

#$(CMD): $(LIBOBJ)
#	$(FC) -o $@ $(LIBOBJ) $(LDFLAGS) 
$(CMD): libMOD.a cwb_letkf.o
	$(RM)    $@
	$(FC) -o $@ cwb_letkf.o libMOD.a $(LDFLAGS)


libMOD.a: $(LIBOBJ)
	@echo " "
	$(RM) $@
	$(AR) $@ $(LIBOBJ)

#-----------------------------------------------------------------------------
.f90.o:
	@echo " "
	$(RM) $@
	$(CPP) $(CPPFLAGS) $*.f90 > $*.F90
	$(FC) -c $(FCFLAGS) $*.F90 $(LDFLAGS)

#-----------------------------------------------------------------------------
clean:
	$(RM) *.o *.mod *.a *.F90

#------dependency
#cwb_letkf.o : module_mpi_util.o module_param.o module_config.o module_grid.o module_eigen.o module_localization.o module_gts_omboma.o module_radar.o module_letkf_core.o
cwb_letkf.o : libMOD.a
module_grid.o : module_param.o module_config.o module_netcdf_io.o module_mpi_util.o
module_gts_omboma.o : module_param.o module_config.o module_mpi_util.o module_projection.o
module_radar.o : module_param.o module_config.o module_mpi_util.o module_projection.o
module_letkf_core.o : module_config.o module_localization.o module_eigen.o module_grid.o module_gts_omboma.o module_radar.o module_mpi_util.o module_projection.o
module_mpi_util.o : module_config.o module_param.o
module_localization.o : module_kdtree2.o module_param.o module_config.o module_gts_omboma.o module_radar.o
module_projection.o : module_config.o module_param.o
