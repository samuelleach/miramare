program make_window
! AUTHOR: S. Leach- based on write_cvbl.pro by J.P.Leahy

      use IniFile
      use mpi_stop
      use HealpixObj
      use HealpixVis
      use beamtools
      use alm_tools

      character(LEN=Ini_max_string_len) InputFile
      logical bad

      integer ll
      integer lmax_beam,lthreshold,lthreshold_varbeam,lthreshold2

      !Parameters
      logical verbose
      real(dp) input_fwhm_arcmin
      real(dp) target_fwhm_arcmin
      integer nside
      character(LEN=256) :: input_beam_file
      character(LEN=256) :: signal_beam_file  
      character(LEN=256) :: variance_beam_file

      !Data
      real(dp), pointer     :: input_beam(:)
      real(dp), allocatable :: gaussian_beam(:,:)
      real(dp), pointer     :: variance_beam(:)
      real(dp), pointer     :: signal_beam(:)
      

      character(len=*), parameter :: CODE = "MAKE_WINDOW"

      !----------------------
      !Read in parameter file
      !----------------------
      InputFile = GetParam(1)
      call Ini_Open(InputFile, 1, bad, .false.)
      if (bad) then
         print *,CODE,': Make smoothing window functions for signal, variance and hitmaps.'
         print *,'---------------------'
         print *,'Parameter file format'
         print *,'---------------------'
         print *,'Use either:'
         print *,'input_fwhm_arcmin = Input beam FWHM in arcmin.'
         print *,'or:'
         print *,'input_beam_file = Filename of two column beam B_l.'
         print *,'then:'
         print *,'verbose = T or F (default = F).'
         print *,'target_fwhm_arcmin = Target FWHM in arcmin for smoothing.'
         print *,'nside = nside value of variance and hit maps to be smoothed.'
         print *,'signal_beam_file = Name of output B_l for smoothing signal maps.'
         print *,'variance_beam_file = Name of output B_l for smoothing variance maps at nside.'
         print *,'---------------------'
         call DoStop('Error opening parameter file')
      endif
      Ini_fail_on_not_found = .false.
      if (InputFile == '') call DoStop('No parameter input file')      
      
      verbose            = ini_read_logical('verbose',.false.)
      input_fwhm_arcmin  = Ini_Read_Real('input_fwhm_arcmin',-1.)
      target_fwhm_arcmin = Ini_Read_Real('target_fwhm_arcmin',60.)
      nside              = Ini_Read_Int('nside',-1)
      input_beam_file    = Ini_Read_String('input_beam_file')
      signal_beam_file   = Ini_Read_String('signal_beam_file')
      variance_beam_file = Ini_Read_String('variance_beam_file')


      if(nside .eq. -1) then
        print *,trim(CODE),': ERROR: Must choose a value for nside (128, 256, 512 etc).'
        stop
      end if

      select case (input_fwhm_arcmin .eq. -1)
      case(.true.)
         if(verbose) print *,trim(CODE),': Reading ',trim(input_beam_file)
         call ReadBeam(input_beam,input_beam_file,lmax_beam)
      case (.false.)
         if(input_beam_file .ne. '') then
            print *,trim(CODE),& 
                 ': ERROR: Choose one of input_beam_file or input_fwhm_arcmin in parameter file.'
            stop
         end if
!         lmax_beam = 3*2048+1
         lmax_beam = 7000
         allocate(gaussian_beam(0:lmax_beam,3))
         allocate(input_beam(0:lmax_beam))
         call gaussbeam(input_fwhm_arcmin,lmax_beam,gaussian_beam)
         input_beam(:) = gaussian_beam(:,1)
         deallocate(gaussian_beam)
      end select
       
      allocate(gaussian_beam(0:lmax_beam,3))
      allocate(signal_beam(0:lmax_beam))

      !----------------------------------------------------------------
      !Make beam for signal smoothing (ratio of output and input beams)
      !----------------------------------------------------------------
      if(verbose) print *,trim(CODE),': Making signal beam'
      call gaussbeam(target_fwhm_arcmin,lmax_beam,gaussian_beam)
      do ll = 0, lmax_beam 
         signal_beam(ll) = input_beam(0)*gaussian_beam(ll,1)/input_beam(ll)
      end do
      call CleanBeam(signal_beam,lmax_beam,1d-6,lthreshold)

      !------------------
      ! Write signal beam
      !------------------
      if(signal_beam_file .ne. '') then
         if(verbose) print *,trim(CODE),': Signal beam lmax = '//trim(inttostr(lthreshold))
         if(verbose) print *,trim(CODE),': Writing signal beam: '//trim(signal_beam_file)
         call WriteBeam(signal_beam,signal_beam_file,lthreshold)
      end if


      !--------------------------------------
      ! Make beam for smoothing variance maps
      !--------------------------------------
      lthreshold_varbeam = 3*lthreshold
      if(lthreshold_varbeam .gt. 3*nside) then
        lthreshold_varbeam = 3*nside
      endif
      if(verbose) print *,trim(CODE),': Setting lmax of variance beam to '//trim(inttostr(lthreshold_varbeam))

      if(verbose) print *,trim(CODE),': Making variance beam'
      allocate(variance_beam(0:lthreshold_varbeam))
      call variancebeam(signal_beam,lmax_beam,nside,variance_beam,lthreshold_varbeam)
           
      !--------------------
      ! Write variance beam
      !--------------------t
      call CleanBeam(variance_beam,lthreshold_varbeam,4.d-4,lthreshold2)
      if(variance_beam_file .ne. '') then
         if(verbose) print *,trim(CODE),': Variance beam lmax = '//&
              trim(inttostr(lthreshold2))
         if(verbose) print *, trim(CODE),': Writing variance beam: ',&
              trim(variance_beam_file)
         call WriteBeam(variance_beam,trim(variance_beam_file),lthreshold2)
      end if

      !----------------------------
      ! Write inverse variance beam
      !----------------------------
      !Make beam for smoothing inverse variance maps eg hitmaps.
!      inversevariance_beam(:) = variance_beam(:)/variance_beam(0)**2
!      call CleanBeam(inversevariance_beam,lmax_beam,4.d-4,lthreshold)
!      if(inversevariance_beam_file .ne. '') then
!         if(verbose) print *, trim(CODE),': Writing inverse variance beam: ',&
!              trim(inversevariance_beam_file)
!         if(verbose) print *,trim(CODE),': lmax = '//&
!              trim(inttostr(lthreshold))
!         !Allow to add comment
!         call WriteBeam(inversevariance_beam,inversevariance_beam_file,lthreshold)
!      end if

!      deallocate(signal_beam)
!      deallocate(gaussian_beam)
!      deallocate(variance_beam)!,inversevariance_beam)

    end program Make_window
    
 
