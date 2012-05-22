OBJS	:=	toms760.o spin_alm_tools.o utils.o HealpixObj.o inifile.o \
		Matrix_utils.o HealpixVis.o mpi_stop.o HealpixTools.o


libmiramare_lenspix.a : $(OBJS)
	$(AR) $(ARFLAGS) libmiramare_lenspix.a $(OBJS)
	$(RANLIB) libmiramare_lenspix.a

%.o : %.F90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f
	$(MPF77) $(FCOMP) -c $<

%.o : %.c
	$(MPCC) $(CCOMP) -c $< 

clean :
	@rm -f *.o *.mod *.MOD *.a *~

