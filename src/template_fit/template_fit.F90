program lincom_map

      use IniFile
      use mpi_stop
      use HealpixObj
      use HealpixVis
      
      character(LEN=Ini_max_string_len) InputFile
      logical bad

      type(HealpixMap)   :: map,mask,weight,templatemap
!      type(HealpixMap)   :: map_temp

      real undef      
      integer mm!,nn,cc
      
      !Parameters
      integer ntemplates
      logical verbose
      character(LEN=100) :: outfile,infile,weightfile,maskfile
      character(256), allocatable :: templatefile(:)
      real(dp) template_threshold
      
      character(len=*), parameter :: CODE = "TEMPLATE_FIT"

      real(dp) a,b
      integer pp

      undef= -1.6375e30
      !----------------------
      !Read in parameter file
      !----------------------
      InputFile = GetParam(1)
      call Ini_Open(InputFile, 1, bad, .false.)
      if (bad) then
         print *,CODE,': Template fitting'
         print *,'---------------------'
         print *,'Parameter file format'
         print *,'---------------------'
         print *,'verbose    = T or F (default = F).'
         print *,'ntemplates = Number of templates to fit (must equal 1 for the moment).'
         print *,'infile     = Input Healpix map filename.'
         print *,'weightfile = Input Healpix weight filename (inverse variance).'
         print *,'maskfile   = Input Healpix mask filename (for excluding pixels from the fit).'
         print *,'outfile    = output (template subtracted) Healpix filename.'
         print *,'templatefile(n) = nth template.'
         print *,'---------------------'
         call DoStop(trim(concat(CODE,': Error opening parameter file')))
      endif
      Ini_fail_on_not_found = .false.
      if (InputFile == '') call DoStop('No parameter input file')      

      verbose    = Ini_Read_logical('verbose',.false.)
      ntemplates = Ini_Read_Int('ntemplates',1)
      infile     = Ini_Read_String('infile')
      weightfile = Ini_Read_String('weightfile')
      maskfile   = Ini_Read_String('maskfile')
      outfile    = Ini_Read_String('outfile')
      template_threshold = Ini_Read_real('template_threshold',-1.)
      allocate(templatefile(ntemplates))
      do mm = 1,ntemplates
         templatefile(mm) = Ini_Read_String_array('templatefile',mm)
      end do

      if(verbose) print *,CODE,': Running'

      !Read in signal map
      if(verbose) print *,CODE,': Reading ',trim(infile)      
      call HealpixMap_read(map,trim(infile))      

      !Read in weight map
      if(verbose) print *,CODE,': Reading ',trim(weightfile)      
      call HealpixMap_read(weight,trim(weightfile))      

      !Read in mask map
      if(verbose) print *,CODE,': Reading ',trim(maskfile)      
      call HealpixMap_read(mask,trim(maskfile))      
      !Add missing pixels to mask
      where( map%TQU(:,1) .eq. undef)  mask%TQU(:,1)= 0.
      where( map%TQU(:,1) .eq. 0.)     mask%TQU(:,1)= 0.


      !Read in first template
      if(verbose) print *,CODE,': Reading ',trim(templatefile(1))      
      call HealpixMap_read(templatemap,trim(templatefile(1)))      

      if( template_threshold .ne. -1) then
         where( templatemap%TQU(:,1) .lt. template_threshold)     mask%TQU(:,1)= 0.
      endif

      templatemap%TQU(:,1) = templatemap%TQU(:,1)

      a = 0.
      b = 0.
      do pp=0,map%npix-1
         if (mask%tqu(pp,1) .gt. 0.) then
            a = a + templatemap%tqu(pp,1)*weight%tqu(pp,1)*map%tqu(pp,1)
            b = b + templatemap%tqu(pp,1)*weight%tqu(pp,1)*templatemap%tqu(pp,1)
         endif
      end do
      if(verbose) print *,CODE,': Template coefficient = ',a/b 

      
!      !Read in other maps and coadd
!      do mm=2,nmaps
!        if(verbose) print *,CODE,': reading ',trim(infile(mm))      
!        call HealpixMap_read(map_temp,trim(infile(mm)))      
!
!        do cc = 1,map_temp%nmaps
!          where( map_temp%TQU(:,cc) .eq. undef)  map_temp%TQU(:,cc)= 0.
!        end do
!        map_temp%TQU      = map_temp%TQU*a(mm) 
!        map_temp%TQU(:,1) = map_temp%TQU(:,1) + c(mm) 
!        
!        do nn=1,map_temp%nmaps
!           if (nn .gt. map%nmaps) cycle
!           map%tqu(:,nn) = map%tqu(:,nn) + map_temp%tqu(:,nn)
!        end do
!      enddo

!      if(verbose) print *,CODE,': writing ',trim(outfile)      
!      call HealpixMap_write(map,trim(outfile),units='')

      

    end program lincom_map
