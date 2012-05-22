program combine_mask

      use IniFile
      use mpi_stop
      use HealpixObj
      use HealpixVis
      
      character(LEN=Ini_max_string_len) InputFile
      logical bad

      type(HealpixMap)   :: mask
      type(HealpixMap)   :: map_temp

      real undef!,threshold      
      integer mm,pp
      
      !Parameters
      integer nmaps
      logical verbose
      character(LEN=100) :: outfile
      character(256), allocatable :: infile(:)
      
      character(len=*), parameter :: CODE = "COMBINE_MASK"

      undef= -1.6375e30
      !----------------------
      !Read in parameter file
      !----------------------
      InputFile = GetParam(1)
      call Ini_Open(InputFile, 1, bad, .false.)
      if (bad) then
         print *,CODE,': Combination of masks (or hitmaps) to make union mask.'
         print *,'---------------------'
         print *,'Parameter file format'
         print *,'---------------------'
         print *,'verbose = T or F (default = F).'
         print *,'nmaps = number of maps to combine.'
         print *,'infile(n) = nth Healpix map filename.'
         print *,'outfile = output Healpix mask filename.'
!         print *,'threshold = mask cutoff - eg for converting hitmaps to masks (default = 1).'
         print *,'----------------------'
         call DoStop(trim(concat(CODE,': Error opening parameter file')))
      endif
      Ini_fail_on_not_found = .false.
      if (InputFile == '') call DoStop('No parameter input file')      

      verbose = Ini_Read_logical('verbose',.false.)
      nmaps   = Ini_Read_Int('nmaps',1)
      outfile = Ini_Read_String('outfile')
      allocate(infile(nmaps))
      do mm=1,nmaps
        infile(mm) = Ini_Read_String_array('infile',mm)
      end do
!      threshold = Ini_Read_real('threshold',undef)
!      if (threshold .eq. undef) threshold = 0.

      if(verbose) print *,CODE,': running'

      !-----------------
      !Read in first map
      !-----------------
      if(verbose) print *,CODE,': reading ',trim(infile(1))      
      call HealpixMap_read(map_temp,trim(infile(1)))   
      
      call Healpixmap_init(mask,nside=map_temp%nside,nmaps=1,pol=.false.)     
      mask%ordering = map_temp%ordering

      where( map_temp%TQU(:,1) .eq. undef)      map_temp%TQU(:,1)=0        
      where( map_temp%TQU(:,1) .ne. 0)          map_temp%TQU(:,1)=1        
!      where( map_temp%TQU(:,1) .le. threshold)  map_temp%TQU(:,1)=0.        
!      where( map_temp%TQU(:,1) .gt. threshold)  map_temp%TQU(:,1)=1.        

      do pp=0,mask%npix-1
         mask%tqu(pp,1) = map_temp%tqu(pp,1)      
      end do

      !-------------------------------
      !Read in other maps and multiply
      !-------------------------------
      do mm=2,nmaps
         if(verbose) print *,CODE,': reading ',trim(infile(mm))      
         call HealpixMap_read(map_temp,trim(infile(mm)))      
         where( map_temp%TQU(:,1) .eq. undef)      map_temp%TQU(:,1)=0        
         where( map_temp%TQU(:,1) .ne. 0)          map_temp%TQU(:,1)=1        
!         where( map_temp%TQU(:,1) .le. threshold)  map_temp%TQU(:,1)=0.        
!         where( map_temp%TQU(:,1) .gt. threshold)  map_temp%TQU(:,1)=1.
         do pp=0,mask%npix-1
            mask%TQU(pp,1) = mask%TQU(pp,1) * map_temp%TQU(pp,1) 
         end do
      enddo
      
      if(verbose) print *,CODE,': writing ',trim(outfile)      
      call HealpixMap_write(mask,trim(outfile),units='')

      

    end program combine_mask
