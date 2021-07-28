SHELL = 	/bin/sh

CMD     = cwb_letkf.exe

#FC      = mpifrtpx100  -Nquickdbg#-Haefosux 
#FC      = mpifrtpx100 -O2 -Kautoobjstack 
 FC      = mpifrtpx100 -Kfast,autoobjstack,temparraystack
FCFLAGS =
LDFLAGS = -SSL2 -L/package/fx100/netcdf-4.1.3/lib -lnetcdff -lnetcdf -lm -I/package/fx100/netcdf-4.1.3/include \
		        -L/package/fx100/hdf5-1.8.9/lib -lhdf5_hl -lhdf5 -lz
#LDFLAGS = -SSL2 -L/users/xa24/PkgFX100/netcdf-4.6.1/lib -lnetcdff -lnetcdf -I/users/xa24/PkgFX100/netcdf-4.6.1/include \
#		        -L/users/xa24/PkgFX100/hdf5-1.8.20/lib -lhdf5_hl -lhdf5 -lm -lz

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
		  module_kdtree.o \
		  module_eigen.o \
		  module_mpi_util.o \
		  module_letkf_core.o \
		  module_localization.o

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
	$(FC) -c $(FCFLAGS) $*.f90 $(LDFLAGS)

#-----------------------------------------------------------------------------
clean:
	$(RM) *.o *.mod *.a

#------dependency
cwb_letkf.o : libMOD.a
module_grid.o : module_param.o module_config.o module_netcdf_io.o module_mpi_util.o
module_gts_omboma.o : module_param.o module_config.o module_mpi_util.o
module_radar.o : module_param.o module_config.o module_mpi_util.o
module_letkf_core.o : module_config.o module_localization.o module_eigen.o module_grid.o module_gts_omboma.o module_radar.o module_mpi_util.o
module_mpi_util.o : module_config.o module_param.o
module_localization.o : module_kdtree.o module_param.o module_config.o module_gts_omboma.o module_radar.o
