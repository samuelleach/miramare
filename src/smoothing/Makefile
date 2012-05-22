all: smoothing

OBJS	:= 

smoothing : $(OBJS) smoothing.o
	$(F90) $(F90COMP) -o $@ smoothing.o $(LINK)
	cp $@ ../../bin/.

%.o : %.F90
	$(F90) $(F90COMP) -c $<

%.o : %.f90
	$(F90) $(F90COMP) -c $<

%.o : %.f
	$(F77) $(FCOMP) -c $<

%.o : %.c
	$(CC) $(CCOMP) -c $< 

clean :
	@rm -f *.o *.mod *.MOD *.a *~ smoothing
