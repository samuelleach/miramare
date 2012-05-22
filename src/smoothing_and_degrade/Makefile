MOBJS	:= beamtools.o mask_maker.o

MODULE	:= smoothing_and_degrade
all: $(MODULE)

$(MODULE) : $(MOBJS) $(MODULE).o
	$(MPF90) -I../make_window -o $(MODULE) $(MOBJS) $(MODULE).o $(LINK) 
	cp $@ ../../bin/.

%.o : %.F90
	$(MPF90) -I../make_window $(F90COMP) -c $<

%.o : %.f90
	$(MPF90) -I../make_window $(F90COMP) -c $<

%.o : %.f
	$(MPF77) $(FCOMP) -c $<

%.o : %.c
	$(MPCC) $(CCOMP) -c $< 

clean :
	@rm -f *.o *.mod *.MOD *.a *~ $(MODULE)
