MOBJS	:=	settings.o paramdef.o propose.o InstrumentData.o\
		CMBData.o miramare_utils.o foregroundscalings.o runparameters.o \
		foreground_simulations.o noise_simulations.o \
		likelihood_CMBFOREGROUND.o calclike.o conjgrad_wrapper.o \
		EstCovmat.o MCMC.o

all:  libmiramare.a miramare_estspec miramare_getdist miramare_estamp

miramare_estspec: libmiramare.a miramare_estspec.o
	$(MPF90) -o $@ miramare_estspec.o -L. -lmiramare $(LINK)
	cp $@ ../../bin/

miramare_estamp: libmiramare.a miramare_estamp.o
	$(MPF90) -o $@ miramare_estamp.o -L. -lmiramare $(LINK)
	cp $@ ../../bin/

miramare_getdist: libmiramare.a miramare_getdist.o
	$(MPF90) -o $@ miramare_getdist.o -L. -lmiramare $(LINK)
	cp $@ ../../bin/ 

libmiramare.a: $(MOBJS)
	$(AR) $(ARFLAGS) libmiramare.a $(MOBJS)
	$(RANLIB) libmiramare.a

%.o: %.F90
	$(MPF90) $(F90COMP) -c $<

%.o: %.f90
	$(MPF90) $(F90COMP) -c $<

%.o: %.f
	$(MPF77) $(FCOMP) -c $<

%.o: %.c
	$(MPCC) $(CCOMP) -c $< 

clean:
	@rm -f *.o *.mod *.MOD *.a *~ miramare_estspec miramare_estamp miramare_getdist