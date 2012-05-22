program maps_to_iqumap

      use IniFile
      use mpi_stop
      use HealpixObj
      use HealpixVis
      
      character(LEN=Ini_max_string_len) InputFile
      logical bad

      type(HealpixMap)   :: map
      type(HealpixMap)   :: map_out

      real undef      
      integer mm
      
      !Parameters
      integer nmaps
      logical verbose
      character(LEN=256) :: outfile
      character(256), allocatable :: infile(:)
      
      character(len=*), parameter :: CODE = "MAPS_TO_IQUMAP"

      undef= -1.6375e30
      !----------------------
      !Read in parameter file
      !----------------------
      InputFile = GetParam(1)
      call Ini_Open(InputFile, 1, bad, .false.)
      if (bad) then
         print *,CODE,': Take I,Q, and U (and Nobs map) and convert to IQU(Nobs) map.'
         print *,'---------------------'
         print *,'Parameter file format'
         print *,'---------------------'
         print *,'verbose = T or F (default = F).'
         print *,'nmaps = number of maps to assemble (default = 3).'
         print *,'outfile = output Healpix mask filename.'
         print *,'infile(1) = Healpix I map filename.'
         print *,'infile(2) = Healpix Q map filename.'
         print *,'infile(3) = Healpix U map filename.'
         print *,'infile(4) = Healpix Nobs map filename.'
         print *,'outfile = Healpix output map filename.'
         print *,'---------------------'
         call DoStop(trim(concat(CODE,': Error opening parameter file')))
      endif
      Ini_fail_on_not_found = .false.
      if (InputFile == '') call DoStop('No parameter input file')      

      verbose = Ini_Read_logical('verbose',.false.)
      nmaps= Ini_Read_Int('nmaps',3)
      outfile = Ini_Read_String('outfile')
      allocate(infile(nmaps))
      do mm=1,nmaps
        infile(mm) = Ini_Read_String_array('infile',mm)
      end do

      if(verbose) print *,CODE,': running'


      !Read in first map
      if(verbose) print *,CODE,': reading ',trim(infile(1))      
      call HealpixMap_read(map,trim(infile(1)))      

      call HealpixMap_init(map_out,nside=map%nside,nmaps=nmaps)
      map_out%tqu(:,1)=map%tqu(:,1)
      
      !Read in other maps
      do mm=2,nmaps
        if(verbose) print *,CODE,': reading ',trim(infile(mm))      
        call HealpixMap_read(map,trim(infile(mm)))      
        map_out%tqu(:,mm)= map%tqu(:,1)
      enddo

      if(verbose) print *,CODE,': writing ',trim(outfile)      
      call HealpixMap_write(map_out,trim(outfile),units='')
      

    end program maps_to_iqumap
