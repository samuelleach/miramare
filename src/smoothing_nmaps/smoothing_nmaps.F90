program smoothing_nmaps

      use IniFile
      use mpi_stop
      use HealpixObj
      use HealpixVis
      use pix_tools
      use beamtools
      use alm_tools
      
      character(LEN=Ini_max_string_len) InputFile
      logical bad,want_pol

      type(HealpixMap)   :: map
      Type(HealpixAlm)   :: alm
      type(HealpixMap)   :: map_tmp
      Type(HealpixInfo)  :: HH
      type(HealpixMap)   :: map_out

      real undef      
      real(dp), pointer     :: beam(:)
      real(dp), allocatable :: gaussian_beam(:,:)

      !Parameters
      logical verbose
      character(LEN=256) :: infile
      character(LEN=256) :: outfile
      character(LEN=256) :: beam_file
      logical discard_offdiagonal
      real(dp) fwhm_arcmin
      real(dp) rescale_factor
      integer nlmax

      character(len=*), parameter :: CODE = "SMOOTHING_NMAPS"
      integer mm

      
#ifdef MPIPIX
      integer mpi_division_method,i     
      call mpi_init(i)
      mpi_division_method = 3
#endif

      
      undef= -1.63750e30
      !----------------------
      !Read in parameter file
      !----------------------
      InputFile = GetParam(1)
      call Ini_Open(InputFile, 1, bad, .false.)
      if (bad) then
         print *,CODE,': Scalar smoothing of Healpix maps (eg 6 column maps).'
         print *,'---------------------'
         print *,'Parameter file format'
         print *,'---------------------'
         print *,'verbose = T or F (default = F).'
         print *,'infile  = Healpix map filename.'
         print *,'outfile = Output Healpix map filename.'
         print *,'nlmax   = Maximum l for smoothing.'
         print *,'discard_offdiagonal   = T or F.'
         print *,'rescale_factor        = Multiplicative factor (default = 1.) '
         print *,'Then use either:'
         print *,'fwhm_arcmin = Beam FWHM in arcmin.'
         print *,'Or:'
         print *,'beam_file = Name of B_l file.'
         print *,'---------------------'
         call DoStop(trim(concat(CODE,': Error opening parameter file')))
      endif
      Ini_fail_on_not_found = .false.
      if (InputFile == '') call DoStop('No parameter input file')      

      verbose         = Ini_Read_logical('verbose',.false.)
      discard_offdiagonal = Ini_Read_logical('discard_offdiagonal',.false.)
      infile          = Ini_Read_String('infile')
      outfile         = Ini_Read_String('outfile')
      nlmax           = Ini_Read_int('nlmax',2500)
      fwhm_arcmin     = Ini_Read_Real('fwhm_arcmin',-1.)
      beam_file       = Ini_Read_String('beam_file')
      rescale_factor  = Ini_Read_Real('rescale_factor',1.)


      if(verbose) print *,trim(CODE),': running'

      if (beam_file .ne. '') then
         call ReadBeam(beam, beam_file, nlmax)
         if (fwhm_arcmin .ne. -1.) then
            print *,trim(CODE),': WARNING: Choose fwhm_arcmin or beam_file parameter.'
         endif
      endif

      if (fwhm_arcmin .ne. -1.) then
         allocate(gaussian_beam(0:nlmax,3))
         call gaussbeam(fwhm_arcmin,nlmax,gaussian_beam)
         beam(:)=gaussian_beam(:,1)
         deallocate(gaussian_beam)
      endif


      if(verbose) print *,trim(CODE),': reading ',trim(infile)
      call HealpixMap_read(map,trim(infile))
      where( abs(map%TQU(:,:)/undef -1.) .lt. 1.e-5 )  map%TQU(:,:) = 0.

      if (map%nmaps .gt. 1) then
         want_pol = .true.
      else
         want_pol = .false.
      endif

      if(verbose) print *,trim(CODE),': Initialise healpix'
      call HealpixInit(HH,map%nside, 2*nlmax,want_pol, w8dir='')

      call HealpixMap_init(map_tmp,map%npix,nmaps=1)
      if (HH%MpiID .eq. 0) then
    
      do mm=1,map%nmaps 
         if(verbose) print *,trim(CODE),': map2alm, map ',mm
         map_tmp%tqu(:,1) = map%tqu(:,mm)
         
         call HealpixMap2alm(HH,map_tmp,alm,nlmax)

         if(verbose) print *,trim(CODE),': smoothing alm ',mm
         call HealpixAlm_smooth_beam(alm,beam)

         if(verbose) print *,trim(CODE),': alm2map ',mm
         call HealpixAlm2Map(HH,alm, map_tmp, map%npix)
         if(verbose) print *,trim(CODE),': rescaling map by factor = ',rescale_factor
         map%tqu(:,mm) = map_tmp%tqu(:,1) *rescale_factor
      end do

      if(discard_offdiagonal .and. map%nmaps .gt. 1) then
         if(verbose) print *,trim(CODE),': Discarding off diagonal terms.'
         call HealpixMap_init(map_out,map%npix,nmaps=3)
         map_out%tqu(:,1) = map%tqu(:,1)
         map_out%tqu(:,2) = map%tqu(:,4)
         map_out%tqu(:,3) = map%tqu(:,6)
      else
         map_out = map
      endif


      if(verbose) print *,trim(CODE),': writing ',trim(outfile)
      call HealpixMap_write(map_out,trim(outfile),units='')
      endif
    
#ifdef MPIPIX
      call HealpixFree(HH)
      call mpi_finalize(i)
#endif
   

    end program smoothing_nmaps
