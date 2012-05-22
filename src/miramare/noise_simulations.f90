module noise_simulations

  use runparameters
  use HealpixObj
  use cmbdata

contains

  subroutine get_noise_simulation(WeightSet,MapSet,nfreq,runpar,writetocache)
    integer, intent(in) :: nfreq
    Type(HealpixMap) :: WeightSet(nfreq)
    Type(HealpixMap) :: MapSet(nfreq)
    type(RunParams) :: runpar
    logical, optional, intent(in) :: writetocache
    
    integer ff,ss,pp
    real dummy   
    Type(HealpixMap) :: noisemap
    character(LEN=512) filename
    character(LEN=120) paramstr
    logical do_sim
    
    call InitRandom(runpar%noise_rand_seed)
    dummy = gaussian1()

    call HealpixMap_Nullify(noisemap)
    call HealpixMap_Init(noisemap, MapSet(1)%npix, pol = runpar%want_pol,&
         nested=orderingisnest(runpar%ordering_in))
    do ff = 1, nfreq
       do_sim = .true.
       write (paramstr,*) ff
       filename = trim(runpar%cachedir)//'/'//trim(runpar%rootname)//&
            '_noise_'//trim(adjustl(paramstr))//'.fits'

       !------------------------------------------------------
       !Read in cached noise simulation or perform simulations
       !------------------------------------------------------
       if(present(writetocache)) then
          if(FileExists(filename)) then 
             if (runpar%feedback > 0) print *, 'Reading ',trim(filename)
             call HealpixMap_read(noisemap,trim(filename))            
             do_sim = .false.
          endif
       endif

       if(do_sim) then
          !Need to insert new simulations depending on whether weight matrix has 6 elements or not.
          if (runpar%Feedback>0) print *,'Simulating the noise for channel ',ff
          do ss = 1, runpar%nstokes
             do pp = 0, MapSet(1)%npix-1 
                noisemap%tqu(pp,runpar%stokes_list(ss)) = & 
                     Gaussian1()/sqrt(WeightSet(ff)%TQU(pp,runpar%stokes_list(ss)))
             end do
          end do
       endif
       
       !-------------------------------------------
       !Add (optionally rescaled) noise to the maps
       !-------------------------------------------
       do ss = 1, runpar%nstokes
          do pp = 0, Mapset(ff)%npix -1
             MapSet(ff)%TQU(pp,runpar%stokes_list(ss)) = MapSet(ff)%TQU(pp,runpar%stokes_list(ss)) + &
                  noisemap%tqu(pp,runpar%stokes_list(ss))*runpar%sim_noise_factor(ff)
          enddo
       end do

       !----------------------------------
       !Write noise maps to the disk cache
       !----------------------------------
       if(do_sim) then
          if (present(writetocache)) then
             if(writetocache) then
                if (runpar%feedback > 0) print *, 'Writing ',trim(filename)
                call HealpixMap_write(noisemap,'!'//trim(filename))            
             end if
          end if
       endif
    end do

    call HealpixMap_Free(noisemap)


 end subroutine get_noise_simulation
  

end module noise_simulations
