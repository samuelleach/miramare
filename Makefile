# Commander Makefile for use with EITHER the planck
# module build system OR a set of stand-alone
# config files.  Default make target prints help.
#
# This makefile requires the use of GNU make or
# compatible tools that allow the substitution
# operator ":=".
#

TOPDIR := $(shell pwd)
export TOPDIR

ifdef TARGET 
	ifndef PREFIX
		PREFIX := ./$(TARGET)
	endif
	INSTALL := $(PREFIX)/miramare
	include $(TOPDIR)/config/config.$(TARGET)
else
	ifdef MIRAMARE
		include $(TOPDIR)/config/config.$(MIRAMARE)
		ifndef INSTALL
			INSTALL := $(TOPDIR)/install_$(MIRAMARE)
		endif
	else
		$(error MIRAMARE undefined)UNDEFINED
	endif
	ifndef MAKE
		export MAKE := make
	endif
	ifndef AR
		export AR := ar
	endif
	ifndef ARFLAGS
		export ARFLAGS := crv
	endif
	ifndef RANLIB
		export RANLIB := ranlib
	endif
	ifndef CTAR
		export CTAR := tar czvf
	endif
	ifndef F90
		export F90 := f90
	endif
	ifndef MPF90
		export MPF90 := mpif90
	endif
	ifndef F90FLAGS
		export F90FLAGS := -g -O2
	endif
	ifndef MPF77
		export MPF77 := mpif77
	endif
	ifndef FFLAGS
		export FFLAGS := -g -O2
	endif
	ifndef MPCC
		export MPCC := cc
	endif
	ifndef CFLAGS
		export CFLAGS := -g -O2
	endif
	ifndef LDFLAGS
		export LDFLAGS := -lm
	endif
	ifndef FORTRAN_UPPER
		export FORTRAN_UPPER := 0
	endif
	ifndef CFITSIO_LINK
		export CFITSIO_LINK := -L/usr/local/lib -lcfitsio
	endif
	ifndef LAPACK_LINK
		export LAPACK_LINK := -L/usr/local/lib -llapack -lblas
	endif
endif

export F90COMP	:= $(F90FLAGS) -I$(TOPDIR)/src/lenspix -I$(TOPDIR)/src/miramare $(LAPACK_INCLUDE) $(CFITSIO_INCLUDE) $(HEALPIX_INCLUDE)
export FCOMP	:= $(FFLAGS)   -I$(TOPDIR)/src/lenspix -I$(TOPDIR)/src/miramare $(LAPACK_INCLUDE) $(CFITSIO_INCLUDE) $(HEALPIX_INCLUDE)
export CCOMP	:= $(CFLAGS)   -I$(TOPDIR)/src/lenspix -I$(TOPDIR)/src/miramare $(LAPACK_INCLUDE) $(CFITSIO_INCLUDE) $(HEALPIX_INCLUDE)
export LINK	:= -L$(TOPDIR)/src/lenspix -lmiramare_lenspix $(CFITSIO_LINK) $(LAPACK_LINK) $(HEALPIX_LINK) $(LDFLAGS)

all:	sqmatrix\
	libmiramare_lenspix\
	newmodule\
	smoothing_nmaps smoothing_and_degrade \
	make_window add_white_noise ud_grade_mask combine_mask\
	lincom_map jackknife qucorr template_fit \
	cl2pf mpiBatch ud_grade_hitmap\
	maps_to_iqumap degrade_map conversion_factor\
	hitmap_to_regionfile miramare_estamp wmap_to_madam_invn fgres\
	degrade_and_smoothing_covmat map_to_dermap


help :
	@echo ' '
	@echo '  This Makefile is used to build Miramare in a way that is'
	@echo '  customized to your system.  You must export the MIRAMARE'
	@echo '  environment variable and set it to the name of your platform.'
	@echo '  Then you must create a config file named config/config.<$MIRAMARE>.'
	@echo '  I suggest copying the example config file and modifying the'
	@echo '  parameters for your platform.'
	@echo ' '
	@echo '  The following make targets are supported:'
	@echo ' '
	@echo '    make         : build everything'
	@echo '    make help    : print this help screen'
	@echo '    make install : install everything'
	@echo '    make clean   : remove build files'
	@echo '    make dist    : construct a date-stamped source tarball'
	@echo ' '

sqmatrix:
	@cd src/$@; $(MAKE)

libmiramare_lenspix:
	@cd src/lenspix; $(MAKE)

libmiramare: libmiramare_lenspix
	@cd src/miramare; $(MAKE)

miramare_estspec: libmiramare_lenspix
	@cd src/miramare; $(MAKE) $@

miramare_estamp: libmiramare_lenspix
	@cd src/miramare; $(MAKE)

miramare_getdist: libmiramare_lenspix
	@cd src/miramare; $(MAKE)

smoothing:
	@cd src/$@; $(MAKE)

newmodule: libmiramare_lenspix
	@cd src/$@; $(MAKE)

smoothing_nmaps: make_window libmiramare_lenspix
	@cd src/$@; $(MAKE)

add_white_noise: libmiramare_lenspix
	@cd src/$@; $(MAKE)

make_window: libmiramare_lenspix
	@cd src/$@; $(MAKE)

ud_grade_mask: libmiramare_lenspix
	@cd src/$@; $(MAKE)

ud_grade_hitmap: libmiramare_lenspix
	@cd src/$@; $(MAKE)

combine_mask: libmiramare_lenspix
	@cd src/$@; $(MAKE)

lincom_map: libmiramare_lenspix
	@cd src/$@; $(MAKE)

qucorr: libmiramare_lenspix
	@cd src/$@; $(MAKE)

jackknife: libmiramare_lenspix
	@cd src/$@; $(MAKE)

template_fit: libmiramare_lenspix
	@cd src/$@; $(MAKE)

maps_to_iqumap: libmiramare_lenspix
	@cd src/$@; $(MAKE)

degrade_map: libmiramare_lenspix
	@cd src/$@; $(MAKE)

conversion_factor: libmiramare_lenspix libmiramare
	@cd src/$@; $(MAKE)

hitmap_to_regionfile: libmiramare_lenspix
	@cd src/$@; $(MAKE)

wmap_to_madam_invn: libmiramare_lenspix
	@cd src/$@; $(MAKE)

smoothing_and_degrade: libmiramare_lenspix
	@cd src/$@; $(MAKE)

degrade_and_smoothing_covmat: libmiramare_lenspix
	@cd src/$@; $(MAKE)

map_to_dermap: libmiramare_lenspix
	@cd src/$@; $(MAKE)

fgres: libmiramare_lenspix
	@cd src/$@; $(MAKE)

cl2pf:
	@cd src/$@; $(MAKE)

mpiBatch:
	@cd src/$@; $(MAKE)

clean:
	@cd src/newmodule; $(MAKE) $@
	@cd src/sqmatrix; $(MAKE) $@
	@cd src/lenspix; $(MAKE) $@
	@cd src/miramare; $(MAKE) $@
	@cd src/smoothing; $(MAKE) $@
	@cd src/smoothing_nmaps; $(MAKE) $@
	@cd src/smoothing_and_degrade; $(MAKE) $@
	@cd src/make_window; $(MAKE) $@
	@cd src/add_white_noise; $(MAKE) $@
	@cd src/ud_grade_mask; $(MAKE) $@
	@cd src/ud_grade_hitmap; $(MAKE) $@
	@cd src/lincom_map; $(MAKE) $@
	@cd src/combine_mask; $(MAKE) $@
	@cd src/maps_to_iqumap; $(MAKE) $@
	@cd src/degrade_map; $(MAKE) $@
	@cd src/conversion_factor; $(MAKE) $@
	@cd src/hitmap_to_regionfile; $(MAKE) $@
	@cd src/wmap_to_madam_invn; $(MAKE) $@
	@cd src/fgres; $(MAKE) $@
	@cd src/degrade_and_smoothing_covmat; $(MAKE) $@
