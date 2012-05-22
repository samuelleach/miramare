MOBJS	:= beamtools.o

all: make_window

make_window : $(MOBJS) make_window.o
	$(MPF90) -o $@ $(MOBJS) make_window.o $(LINK)
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
	@rm -f *.o *.mod *.MOD *.a *~ make_window
