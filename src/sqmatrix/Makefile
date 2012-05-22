all: sqmatrix

OBJS	:= matrix.o fits_routines.o

sqmatrix : $(OBJS) sqmatrix.o
	$(F90) -o $@ $@.o $(OBJS) $(CFITSIO_LINK)
	cp $@ ../../bin/.

%.o : %.F90
	$(F90) $(F90COMP) -c $<

%.o : %.f90
	$(F90) $(F90COMP) -c $<

%.o : %.f
	$(MPF77) $(FCOMP) -c $<

%.o : %.c
	$(MPCC) $(CCOMP) -c $< 

clean :
	@rm -f *.o *.mod *.MOD *.a *~ sqmatrix
