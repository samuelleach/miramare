MOBJS	:=

MODULE	:= jackknife
all: $(MODULE)

$(MODULE) : $(MODULE).o
	$(MPF90) -o $(MODULE) $(MODULE).o $(LINK)
	cp $@ ../../bin/.

%.o : %.F90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f
	$(F77) $(FCOMP) -c $<

%.o : %.c
	$(CC) $(CCOMP) -c $< 

clean :
	@rm -f *.o *.mod *.MOD *.a *~ $(MODULE)
