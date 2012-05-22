program conversion_factor

      use IniFile
      use mpi_stop
      use HealpixObj !Need to understand how to do without this.
      use foregroundscalings


      character(LEN=Ini_max_string_len) InputFile
      logical bad

      !Parameters
      character(LEN=100) :: from_unit
      character(LEN=100) :: to_unit
      real*8 freq
      logical verbose

      
      real*8 pref_val_in,pref_val_out
      real*8 factor

      character(len=*), parameter :: CODE = "CONVERSION_FACTOR"
      
      !----------------------
      !Read in parameter file
      !----------------------
      InputFile = GetParam(1)
      call Ini_Open(InputFile, 1, bad, .false.)
      if (bad) then 
         print *,CODE,': Unit conversion module.'
         print *,CODE,': Inspired by conversion_factor.pro in the PSM.'
         print *,'---------------------'
         print *,'Parameter file format'
         print *,'---------------------'
         print *,'verbose = T or F (default = F).'
         print *,'from_unit = K_CMB or K_RJ.'
         print *,'to_unit = K_CMB or K_RJ.'
         print *,'freq = Frequency in GHz..'
         print *,'---------------------'
         call DoStop(trim(concat(CODE,': Error opening parameter file')))
      endif
      Ini_fail_on_not_found = .false.
      if (InputFile == '') call DoStop('No parameter input file')      

      verbose = Ini_Read_logical('verbose',.false.)
      from_unit = Ini_Read_string('from_unit')
      to_unit = Ini_Read_string('to_unit')
      freq = Ini_Read_real('freq',23.)

      if(verbose) print *,CODE//': Running'
      if(verbose) print *,CODE//': Converting from'//trim(from_unit)//&
           ' to '//trim(to_unit)

      pref_val_in=1d0
      pref_val_out=1d0
      

      !QUICK AND DIRTY VERSION FOR NOW

!      select case(from_unit)
!         case('K_CMB'):
!         case('K_RJ'):            
!      end select
!
!      select case(to_unit)
!         case('K_CMB'):
!         case('K_RJ'):            
!      end select
 
      if ((trim(from_unit) .eq. 'K_CMB') .and. (trim(to_unit) .eq. 'K_RJ')) factor = S_CMB(freq)
      if ((trim(from_unit) .eq. 'K_RJ') .and. (trim(to_unit) .eq. 'K_CMB')) factor = 1./S_CMB(freq)
      if (trim(from_unit) .eq. trim(to_unit)) factor = 1.
      
      write (*,1) factor
1     format(1F10.5)
      

 
    end program conversion_factor
