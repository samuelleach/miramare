program ud_grade_hitmap

      use IniFile
      use mpi_stop
      use HealpixObj
      use HealpixVis
      use pix_tools
      
      character(LEN=Ini_max_string_len) InputFile
      logical bad

      type(HealpixMap)   :: hitmap
      type(HealpixMap)   :: hitmap_out

      real undef      

      !Parameters
      logical verbose
      character(LEN=256) :: infile
      character(LEN=256) :: outfile
      integer nside

      character(len=*), parameter :: CODE = "UD_GRADE_HITMAP"
      
      
      undef= -1.63750e30
      !----------------------
      !Read in parameter file
      !----------------------
      InputFile = GetParam(1)
      call Ini_Open(InputFile, 1, bad, .false.)
      if (bad) then
         print *,CODE,': Up or degrading of healpix hitmaps or N^-1 maps via summation.'
         print *,'---------------------'
         print *,'Parameter file format'
         print *,'---------------------'
         print *,'verbose = T or F (default = F).'
         print *,'infile = Healpix hitmap filename.'
         print *,'outfile = Output Healpix hitmap filename.'
         print *,'nside = Output Healpix hitmap nside value (default = 512)'
         print *,'---------------------'
         call DoStop(trim(concat(CODE,': Error opening parameter file')))
      endif
      Ini_fail_on_not_found = .false.
      if (InputFile == '') call DoStop('No parameter input file')      

      verbose = Ini_Read_logical('verbose',.false.)
      infile = Ini_Read_String('infile')
      outfile = Ini_Read_String('outfile')
      nside = Ini_Read_Int('nside',512)

      if(verbose) print *,trim(CODE),': running'


      if(verbose) print *,trim(CODE),': reading ',trim(infile)
      call HealpixMap_read(hitmap,trim(infile))
      where( abs(hitmap%TQU(:,:)/undef -1.) .lt. 1.e-5 )  hitmap%TQU(:,:) = 0.

      !Ud grade hitmap
      if(verbose) print *,trim(CODE),': accumulating hitmap from nside ',&
           hitmap%nside,' to ',nside
      call HealpixMap_udgrade(hitmap, hitmap_out, nside)
      
      ! Multiply by "accumulation factor"
      hitmap_out%TQU(:,:) = hitmap_out%TQU(:,:)*(float(hitmap%nside)/float(nside))**2

      if(verbose) print *,trim(CODE),': writing ',trim(outfile)
      call HealpixMap_write(hitmap_out,trim(outfile),units='')

      

    end program ud_grade_hitmap
