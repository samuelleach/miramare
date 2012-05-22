program wmap_to_madam_invn

      use IniFile
      use mpi_stop
      use HealpixObj
      use HealpixVis
      use pix_tools
      use MatrixUtils
      
      character(LEN=Ini_max_string_len) InputFile
      logical bad

      type(HealpixMap)   :: hitmap
      type(HealpixMap)   :: map

      real undef      

      !Parameters
      logical verbose
      character(LEN=100) :: infile
      character(LEN=100) :: outfile
      real(dp) wmap_sigma0_temp
      real(dp) wmap_sigma0_pol
      character(LEN=4) :: ordering_out

      character(len=*), parameter :: CODE = "WMAP_TO_MADAM_INVN"
      
      
      undef= -1.63750e30
      !----------------------
      !Read in parameter file
      !----------------------
      InputFile = GetParam(1)
      call Ini_Open(InputFile, 1, bad, .false.)
      if (bad) then
         print *,CODE,': Converts WMAP noise covariance information into 6 column Madam format.'
         print *,'---------------------'
         print *,'Parameter file format'
         print *,'---------------------'
         print *,'verbose = T or F (default = F).'
         print *,'infile  = WMAP map filename.'
         print *,'outfile = Output MADAM style inverse noise covariance filename.'
         print *,'wmap_sigma0_temp    = For rescaling WMAP hit counts to TEMP^-2 units.'
         print *,'wmap_sigma0_pol  = For rescaling WMAP hit counts to TEMP^-2 units (Stokes Q and U).'
         print *,'ordering_out     = NEST or RING (default NEST)'
         print *,'---------------------'
         call DoStop(trim(concat(CODE,': Error opening parameter file')))
      endif
      Ini_fail_on_not_found = .false.
      if (InputFile == '') call DoStop('No parameter input file')      

      verbose = Ini_Read_logical('verbose',.false.)
      infile = Ini_Read_String('infile')
      outfile = Ini_Read_String('outfile')
      ordering_out = Ini_Read_String('ordering_out')
      wmap_sigma0_temp = Ini_Read_real('wmap_sigma0_temp',1.)
      wmap_sigma0_pol = Ini_Read_real('wmap_sigma0_pol',1.)

      if(ordering_out .eq. '') ordering_out ='NEST'

      if(verbose) print *,trim(CODE),': running'

      if(verbose) print *,trim(CODE),': reading hitcounts from ',trim(infile)
      call HealpixMap_read(hitmap,trim(infile),extno=1)
      
      if (hitmap%nmaps .ne. 4) then
         print *,map%nmaps
         call DoStop(trim(concat(CODE,': Expecting 4 columns in extension 1.')))            
         stop
      end if

      call HealpixMap_init(map,hitmap%npix,nmaps=6)

      if(verbose) print *,trim(CODE),': Converting WMAP hitcounts to Madam format '
      map%tqu(:,1)= hitmap%tqu(:,1)/wmap_sigma0_temp**2
      map%tqu(:,2)= 0.
      map%tqu(:,3)= 0.
      map%tqu(:,4)= hitmap%tqu(:,2)/wmap_sigma0_pol**2
      map%tqu(:,5)= hitmap%tqu(:,3)/wmap_sigma0_pol**2
      map%tqu(:,6)= hitmap%tqu(:,4)/wmap_sigma0_pol**2


      map%ordering = ord_nest  !Hard wired for now. Seems problems reading in the ordering.
      if(verbose) print *,trim(CODE),': Reordering to ',&
           trim(ordering_out)!,' from ',map%ordering
      select case(trim(ordering_out))
         case('RING')
            call HealpixMap_ForceRing(map)
         case('NEST')
            call HealpixMap_ForceNest(map)
         case DEFAULT
            print *,trim(CODE),': ordering out must be NEST or RING'
       end select

      
      if(verbose) print *,trim(CODE),': writing ',trim(outfile)
      call HealpixMap_write(map,trim(outfile),units='')

      

    end program wmap_to_madam_invn
