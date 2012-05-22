program map_to_dermap

      use IniFile
      use mpi_stop
      use HealpixObj
      use HealpixVis
      use alm_tools
      
      character(LEN=Ini_max_string_len) InputFile
      logical bad

      type(HealpixMap)   :: map
      type(HealpixMap)   :: map_out, map_out_split
      Type(HealpixInfo)  :: HH
      Type(HealpixAlm)   :: alm
            

      real undef      
      
      !Parameters
      integer nder
      logical verbose,polmap,splitmap
      character(LEN=256) :: outfile , infile, outfile_split
      real*8 fwhm_arcmin, fwhm_deg
      
      character(len=*), parameter :: CODE = "MAP_TO_DERMAP"
      integer nmapin, nmapout, nlmax, i, nsmax , nmmax, mm
      character(4) :: pout_str


#ifdef MPIPIX
      integer mpi_division_method
      call mpi_init(i)
      mpi_division_method = 3
#endif

      undef= -1.6375e30
      !----------------------
      !Read in parameter file
      !----------------------
      InputFile = GetParam(1)
      call Ini_Open(InputFile, 1, bad, .false.)
      if (bad) then
         print *,CODE,': Take a Healpix map and calculate derivatives of map'
         print *,CODE,': outputting in synfast format'
         print *,'---------------------'
         print *,'Parameter file format'
         print *,'---------------------'
         print *,'verbose     = T or F (default = F).'
         print *,'nder        = Number of derivatives to take (1 or 2, default 1).'
         print *,'infile      = Healpix input map filename.'
         print *,'polmap      = T or F - whether infile is polarized (default = F)'
         print *,'nlmax       = lmax for SHTs.'        
         print *,'fwhm_arcmin = Beam FWHM in arcmin (optional).'
         print *,'outfile     = Healpix output map filename (or file root for splitmap=T).'
         print *,'splitmap    = T or F - whether to split the map into N healpix maps (default = F)'         
         print *,'---------------------'
         call DoStop(trim(concat(CODE,': Error opening parameter file')))
      endif
      Ini_fail_on_not_found = .false.
      if (InputFile == '') call DoStop('No parameter input file')      

      verbose     = Ini_Read_logical('verbose',.false.)
      polmap      = Ini_Read_logical('polmap',.false.)
      splitmap    = Ini_Read_logical('splitmap',.false.)
      nder        = Ini_Read_Int('nder',1)
      nlmax       = Ini_Read_Int('nlmax',-1)
      fwhm_arcmin = Ini_Read_real('fwhm_arcmin',-1.)      
      outfile     = Ini_Read_String('outfile')
      infile      = Ini_Read_String('infile')

      if(verbose) print *,CODE,': running'

      !------------
      ! Read in map
      !------------
      if(verbose) print *,CODE,': reading ',trim(infile)      
      call HealpixMap_read(map,trim(infile))
      if (nlmax .eqv. -1) then
         nlmax = 3 * map%nside - 1
      end if

      call HealpixInit(HH,map%nside, nlmax, polmap, w8dir='')

      !---------------
      ! Get alm of map
      !---------------
      if(verbose) print *,trim(CODE),': map2alm, nlmax = '//trim(inttostr(nlmax))
      call HealpixMap2alm(HH,map,alm,nlmax,dopol=polmap)

      !-------------------
      ! Optional smoothing
      !-------------------
      if (fwhm_arcmin .ne. -1.) then         
         if(verbose) print *,trim(CODE),': Smoothing alm'
         fwhm_deg = fwhm_arcmin/60.
         call HealpixAlm_smooth(alm,fwhm_deg)
      endif


      !-----------------------
      ! Get derivatives of map
      !-----------------------
      if (polmap) then
         nmapin = 3
      else 
         nmapin = 1
      endif
      nmapout = nmapin * 3 * nder
      call HealpixMap_init(map_out,nside=map%nside,nmaps=nmapout)
      
      
      ! Interface taken fom Healpix team synfast code.
       if(verbose) print *,trim(CODE),': Calculating map and derivatives'
       nsmax = map%nside
       nmmax = nlmax
       select case ( (nmapin-1)*10 + nder*100)
       case(100)
          call alm2map_der(nsmax,nlmax,nmmax,alm%TEB,map_out%TQU(:,1),map_out%TQU(:,2:3))
       case(200)
          call alm2map_der(nsmax,nlmax,nmmax,alm%TEB,map_out%TQU(:,1),map_out%TQU(:,2:3),map_out%TQU(:,4:6))
       case(120)
          call alm2map_der(nsmax,nlmax,nmmax,alm%TEB,map_out%TQU(:,1:3),map_out%TQU(:,4:9))
       case(220)
          call alm2map_der(nsmax,nlmax,nmmax,alm%TEB,map_out%TQU(:,1:3),map_out%TQU(:,4:9),map_out%TQU(:,10:18))
       end select

      !----------
      ! Write map
      !----------
      if (splitmap .eqv. .false.) then
         if(verbose) print *,CODE,': writing ',trim(outfile)      
         call HealpixMap_write(map_out,trim(outfile),units='')
      else
         call HealpixMap_init(map_out_split,nside=map%nside,nmaps=1)
         do mm = 1, nmapout
            write(pout_str,'(i0.4)') mm - 1
            outfile_split = trim(outfile) // '_' // pout_str// '.fits'
            if(verbose) print *,CODE,': writing ',trim(outfile_split)      
            map_out_split%tqu(:,1) = map_out%tqu(:,mm)
            call HealpixMap_write(map_out_split,trim(outfile_split),units='')         
         end do
      end if


#ifdef MPIPIX
      call HealpixFree(HH)
      call mpi_finalize(i)
#endif


      

    end program map_to_dermap
