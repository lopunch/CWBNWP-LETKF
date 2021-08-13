SHELL = 	/bin/sh

CMD     = cwb_letkf.exe
CPP     = cpp -C -P -traditional -Wno-invalid-pp-token -ffreestanding

#FC      = mpifrtpx100  -Nquickdbg#-Haefosux 
#FC      = mpifrtpx100 -O2 -Kautoobjstack 
 FC      = mpifrtpx100 -Kfast,autoobjstack,temparraystack
CPPFLAGS = #-DREAL64 #-DNC4
FCFLAGS  =
#LDFLAGS  = -SSL2 -L/package/fx100/netcdf-4.1.3/lib -lnetcdff -lnetcdf -lm -I/package/fx100/netcdf-4.1.3/include \
		         -L/package/fx100/hdf5-1.8.9/lib -lhdf5_hl -lhdf5 -lz
 LDFLAGS = -SSL2 -I/data2/datusers/xa09/fx100/netcdf-4.6.1/include \
                 -L/data2/datusers/xa09/fx100/netcdf-4.6.1/lib -lnetcdff -lnetcdf \
                 -L/data2/datusers/xa09/fx100/hdf5-1.8.22_static/lib -lhdf5 -lhdf5_hl -lhdf5 \
                 -L/users/xa09/fx100/zlib-1.2.8 -lz

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
		  module_localization.o \

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
module_gts_omboma.o : module_param.o module_config.o module_mpi_util.o
module_radar.o : module_param.o module_config.o module_mpi_util.o
module_letkf_core.o : module_config.o module_localization.o module_eigen.o module_grid.o module_gts_omboma.o module_radar.o module_mpi_util.o
module_mpi_util.o : module_config.o module_param.o
module_localization.o : module_kdtree.o module_param.o module_config.o module_gts_omboma.o module_radar.o
