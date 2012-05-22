MOBJS	:=

MODULE	:= cl2pf
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

%.o : %.c
	$(MPCC) $(CCOMP) -c $< 

clean :
	@rm -f *.o *.mod *.MOD *.a *~ $(MODULE)
