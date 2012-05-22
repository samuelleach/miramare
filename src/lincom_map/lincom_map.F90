program lincom_map

      use IniFile
      use mpi_stop
      use HealpixObj
      use HealpixVis
      
      character(LEN=Ini_max_string_len) InputFile
      logical bad

      type(HealpixMap)   :: map
      type(HealpixMap)   :: degraded_map
      type(HealpixMap)   :: map_temp
      type(HealpixMap)   :: covmat

      real undef      
      integer pp,mm,nn,cc,index
      
      !Parameters
      integer nmaps,nside_out
      logical verbose, remove_temperature
      character(LEN=100) :: outfile,covmat_infile,variance_outfile
      real(dp), allocatable :: a(:)
      real(dp), allocatable :: c(:)
      character(256), allocatable :: infile(:)

      !Data
      real(dp), dimension(:,:), allocatable :: covmat_pix
      real(dp) variance_pix
      
      
      character(len=*), parameter :: CODE = "LINCOM_MAP"

      undef= -1.6375e30
      !----------------------
      !Read in parameter file
      !----------------------
      InputFile = GetParam(1)
      call Ini_Open(InputFile, 1, bad, .false.)
      if (bad) then
         print *,CODE,': Linear combination of Healpix maps: m_tot = a_i x m_i + c_i'
         print *,'---------------------'
         print *,'Parameter file format'
         print *,'---------------------'
         print *,'verbose   = T or F (default = F).'
         print *,'nmaps     = number of maps to combine.'
         print *,'outfile   = output Healpix mask filename.'
         print *,'infile(n) = nth Healpix map filename.'
         print *,'a(n)      = nth multiplicative coefficient (default = 1.).'
         print *,'c(n)      = nth temperature offset (default = 0.).'
         print *,'nside     = nside of output maps (Optional, default = nside of input maps)'
         print *,'covmat_infile      = Covariance matrix C of input maps (Optional)'
         print *,'variance_outfile   = Variance of linear combination ( a^t C a ) (Optional)'        
         print *,'remove_temperature = Set temperature map to zero (default = false)'        
         print *,'---------------------'
         call DoStop(trim(concat(CODE,': Error opening parameter file')))
      endif
      Ini_fail_on_not_found = .false.
      if (InputFile == '') call DoStop('No parameter input file')      

      verbose          = Ini_Read_logical('verbose',.false.)
      nmaps            = Ini_Read_Int('nmaps',1)
      nside_out        = Ini_Read_Int('nside',-1)
      outfile          = Ini_Read_String('outfile')
      covmat_infile    = Ini_Read_String('covmat_infile')
      variance_outfile = Ini_Read_String('variance_outfile')
      allocate(infile(nmaps),a(nmaps),c(nmaps))
      do mm = 1, nmaps
         infile(mm) = Ini_Read_String_array('infile',mm)
         a(mm)      = Ini_Read_Real_array('a',mm,1.)
         c(mm)      = Ini_Read_Real_array('c',mm,0.)
      end do
      remove_temperature = Ini_Read_logical('remove_temperature',.false.)

      if(verbose) print *,CODE,': running'
      do mm = 1, nmaps
         if(verbose) print *,CODE,': a, c = ',a(mm),c(mm)
      end do

      !-------------------------------
      !Read in first map and get nside
      !-------------------------------
      if(verbose) print *,CODE,': reading ',trim(infile(1))      
      call HealpixMap_read(map,trim(infile(1)))      
      if (nside_out .eq. -1) then
         nside_out = map%nside
      endif


      do cc = 1,map%nmaps
         where( map%TQU(:,cc) .eq. undef)  map%TQU(:,cc)= 0.
      end do
      map%TQU      = map%TQU*a(1) 
      !--------------------------
      ! Add offsets only to T map
      !--------------------------
      do pp = 0, map%npix-1         
         map%TQU(pp,1) = map%TQU(pp,1) + c(1) 
      end do

      !----------------------------
      !Read in other maps and coadd
      !----------------------------
      do mm = 2, nmaps
        if(verbose) print *,CODE,': reading ',trim(infile(mm))      
        call HealpixMap_read(map_temp,trim(infile(mm)))

        if (map_temp%ordering .ne. map%ordering) then
        !Force the ordering of the maps to be the same         
           if (map%ordering .eq. ord_ring) then
              if(verbose) print *,CODE,': Forcing to ring order map ',trim(infile(mm))      
              call HealpixMap_forcering(map_temp)
           else
              if(verbose) print *,CODE,': Forcing to nest order map ',trim(infile(mm))      
              call HealpixMap_forcenest(map_temp)
           endif                   
        end if

        do cc = 1, map_temp%nmaps
           where( map_temp%TQU(:,cc) .eq. undef)  map_temp%TQU(:,cc)= 0.
        end do


        map_temp%TQU      = map_temp%TQU*a(mm) 
        !--------------------------
        ! Add offsets only to T map
        !--------------------------
        do pp = 0, map%npix-1
           map_temp%TQU(pp,1) = map_temp%TQU(pp,1) + c(mm) 
        end do

        do nn = 1, map_temp%nmaps
           if (nn .gt. map%nmaps) cycle
           do pp = 0, map%npix-1
              map%tqu(pp,nn) = map%tqu(pp,nn) + map_temp%tqu(pp,nn)
           end do
        end do
      enddo

      !---------------------------------------
      ! Optionally set temperature map to zero
      !---------------------------------------
      if(remove_temperature) then
         map%tqu(:,1) = 0.
      end if
      

      !-------------------------
      ! U/degrade map and output
      !-------------------------
      if (nside_out .ne. map%nside) then
         if(verbose) print *,trim(CODE),': U/degrading map from nside ',&
              map%nside,' to ',nside_out
         call HealpixMap_udgrade(map, degraded_map, nside_out)
         if(verbose) print *,CODE,': writing ',trim(outfile)      
         call HealpixMap_write(degraded_map,trim(outfile),units='')      
      else
         if(verbose) print *,CODE,': writing ',trim(outfile)      
         call HealpixMap_write(map,trim(outfile),units='')      
      endif

      !---------------------------------------------------------
      ! Calculation of variance of linear combination (optional)
      !---------------------------------------------------------
      if (covmat_infile .ne. '') then
         if(verbose) print *,CODE,': Performing calculation of covariance '

         if(verbose) print *,CODE,': Reading ',trim(covmat_infile)
         call HealpixMap_read(covmat,trim(covmat_infile))

         allocate(covmat_pix(nmaps,nmaps))

         do pp = 0, covmat%npix - 1            
            index = 0
            do mm = 1, nmaps
               do nn = mm, nmaps
                  index = index + 1
                  covmat_pix(mm,nn) = covmat%tqu(pp,index)
                  covmat_pix(nn,mm) = covmat%tqu(pp,index)
               end do
            end do

            variance_pix  =  dot_product(matmul(a,covmat_pix),a)
            map%tqu(pp,1) = variance_pix
         end do
         map%nmaps = 1 ! Handles intensity only case for now - no polarization

         
         if (nside_out .ne. map%nside) then
            if(verbose) print *,trim(CODE),': Inverting variance map and accumulating from nside ',&
                 map%nside,' to ',nside_out

            !U/degrade inverse variance and convert back to a variance
            do pp = 0, covmat%npix - 1            
               map%tqu(pp,1) = 1./map%tqu(pp,1)            
            end do
            call HealpixMap_udgrade(map, degraded_map, nside_out)
            do pp = 0, degraded_map%npix - 1            
               degraded_map%tqu(pp,1) = 1./degraded_map%tqu(pp,1)*(float(nside_out)/float(map%nside))**2
            end do
            if(verbose) print *,CODE,': writing variance map ',trim(variance_outfile)      
            call HealpixMap_write(degraded_map,trim(variance_outfile),units='')                  
         else
            if(verbose) print *,CODE,': writing variance map ',trim(variance_outfile)      
            call HealpixMap_write(map,trim(variance_outfile),units='')      
         endif
            
      endif


    end program lincom_map
