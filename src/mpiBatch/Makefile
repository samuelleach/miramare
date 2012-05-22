MOBJS	:=

MODULE	:= mpiBatch
all: $(MODULE)

$(MODULE) : $(MODULE).o
	$(MPCC) -o $(MODULE) $(MODULE).o
	cp $@ ../../bin/.

%.o : %.F90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f
	$(F77) $(FCOMP) -c $<

#%.o : %.c
#	$(MPCC) -static_mpi -mt_mpi -c $< 
%.o : %.c
	$(MPCC) -c $< 

clean :
	@rm -f *.o *.mod *.MOD *.a *~ $(MODULE)
