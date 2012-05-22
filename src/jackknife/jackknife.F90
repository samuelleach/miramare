program jackknife

      use IniFile
      use mpi_stop
      use HealpixObj
      use HealpixVis
      use HealpixTools
      
      character(LEN=Ini_max_string_len) InputFile
      logical bad

      type(HealpixMap)   :: map
      type(HealpixMap)   :: mask,degraded_mask
      type(HealpixMap)   :: degraded_map
      type(HealpixMap)   :: map_temp
      type(HealpixMap)   :: covmat,covmat_tmp,degraded_covmat
      type(HealpixMap)   :: chi2

      real undef      
      integer pp,mm,ii,nn,cc,index,nmap_covmat,nmap_chi2
      
      !Parameters
      integer nmaps,nside_out,nnoisesim,noisesim_randseed
      logical verbose, pessimistic_degrade
      character(LEN=100) :: inmaskfile, fileroot_out
      real(dp), allocatable :: a(:)
      real(dp), allocatable :: c(:)
      integer, allocatable :: covmat_firstcol(:)
      character(256), allocatable :: infile(:),covmat_infile(:),covmat_format(:)

      !Data
      real(dp), dimension(:,:), allocatable :: covmat_pix
      
      !Output      
      character(LEN=100) :: outfile,fileextension_out

      character(len=*), parameter :: CODE = "JACKKNIFE"

      logical nested

      undef= -1.6375e30
      !----------------------
      !Read in parameter file
      !----------------------
      InputFile = GetParam(1)
      call Ini_Open(InputFile, 1, bad, .false.)
      if (bad) then
         print *,CODE,': Linear combination of Healpix maps and covariances: m_tot = a_i x m_i + c_i'
         print *,CODE,': followed by noise simulations given by the total covariance.'
         print *,CODE,': Typically set a(1) = 1. and a(2) = -1.'
         print *,'---------------------'
         print *,'Parameter file format'
         print *,'---------------------'
         print *,'verbose   = T or F (default = F).'
         print *,'nmaps     = number of maps to combine.'
         print *,'inmaskfile= input Healpix mask filename (optional)'
         print *,'infile(n) = nth Healpix map filename.'
         print *,'a(n)      = nth multiplicative coefficient (default = 1.).'
         print *,'c(n)      = nth temperature offset (default = 0.).' 
         print *,'nside     = nside of output maps (Optional, default = nside of input maps)'
         print *,'pessimistic_degrade = Handling of empty pixels in mask degrade (Optional, default = T)'
         print *,'covmat_infile(n)  = Covariance matrix C of input maps'
         print *,'covmat_format(n)  = Format of covariance matrix (i, iqu_uncorr, iqu_corr, rms_to_iqu)'
         print *,'covmat_firstcol(n)= Index of covmat first column'
         print *,'nnoisesim         = number of noise sims to make (default = 0)'
         print *,'noisesim_randseed = random number for first noise simulations (default = 123)'
         print *,'fileroot_out      = root name for output files.'
         print *,'---------------------'
         call DoStop(trim(concat(CODE,': Error opening parameter file')))
      endif
      Ini_fail_on_not_found = .false.
      if (InputFile == '') call DoStop('No parameter input file')      

      verbose          = Ini_Read_logical('verbose',.false.)
      pessimistic_degrade = Ini_Read_logical('pessimistic_degrade',.true.)
      nmaps            = Ini_Read_Int('nmaps',1)
      nside_out        = Ini_Read_Int('nside',-1)
      fileroot_out     = Ini_Read_String('fileroot_out')
      allocate(infile(nmaps),a(nmaps),c(nmaps),covmat_infile(nmaps),covmat_firstcol(nmaps))
      allocate(covmat_format(nmaps))
      do mm = 1, nmaps
         infile(mm)          = Ini_Read_String_array('infile',mm)
         covmat_infile(mm)   = Ini_Read_String_array('covmat_infile',mm)
         covmat_format(mm)   = Ini_Read_String_array('covmat_format',mm)
         a(mm)               = Ini_Read_Real_array('a',mm,1.)
         c(mm)               = Ini_Read_Real_array('c',mm,0.)
         covmat_firstcol(mm) = Ini_Read_Int_array('covmat_firstcol',mm,1)
      end do
      nnoisesim            = Ini_Read_Int('nnoisesim',0)
      noisesim_randseed    = Ini_Read_Int('noisesim_randseed',123)

      if(verbose) print *,CODE,': running'
      do mm = 1, nmaps
         if(verbose) print *,CODE,': a, c = ',a(mm),c(mm)
      end do

      
      
      !-------------------------------
      !Read in first map and get nside
      !-------------------------------
      if(verbose) print *,CODE,': reading map ',trim(infile(1))      
      call HealpixMap_read(map,trim(infile(1)))      
      if (nside_out .eq. -1) then
         nside_out = map%nside
      endif
      do cc = 1,map%nmaps
         where( map%TQU(:,cc) .eq. undef)  map%TQU(:,cc)= 0.
      end do
      map%TQU      = map%TQU*a(1) 

      nested = .true.
      if (map%ordering .eq. ord_ring) nested=.false.

      !--------------------------
      ! Add offsets only to T map
      !--------------------------
      do pp = 0, map%npix-1         
         map%TQU(pp,1) = map%TQU(pp,1) + c(1) 
      end do

      !------------
      ! Set up mask
      !------------
      call HealpixMap_init(mask,nside2npix(map%nside),1,nested=nested) 
      mask%TQU(:,1) = 1.
      where( map%TQU(:,1) .eq. 0.)  mask%TQU(:,1)= 0.

      !----------------------------
      !Read in other maps and coadd
      !----------------------------
      do mm = 2, nmaps
        if(verbose) print *,CODE,': reading map ',trim(infile(mm))      
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
        !-------------
        !Build up mask
        !-------------
        if (sum(map_temp%TQU(:,1)) .gt. 0. ) then  !(Handle case where T map is all zeroes)
           where( map_temp%TQU(:,1) .eq. 0.)  mask%TQU(:,1)= 0.
        endif

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
              !-----------------------------
              !Zero under mask (inefficient)
              !-----------------------------
              if (mask%tqu(pp,1) .lt. 1.) then
                 map%tqu(pp,nn) = 0.
              endif

           end do
        end do
      enddo


      !------------------------
      ! Make filename extension
      !------------------------
      if (nside_out .lt. map%nside) then
         fileextension_out =  '_'//trim(inttostr(nside_out))// '.fits'
      else
         fileextension_out =  '_'//trim(inttostr(map%nside))// '.fits'
      endif


      !-----------------------------------
      ! U/degrade map and mask, and output
      !-----------------------------------
      if (nside_out .ne. map%nside) then
         if(verbose) print *,trim(CODE),': U/degrading mask from nside ',&
              map%nside,' to ',nside_out
         call HealpixMap_udgrade(mask, degraded_mask, nside_out,pessimistic=pessimistic_degrade)
         where( degraded_mask%TQU(:,1) .lt. 1.)  degraded_mask%TQU(:,1) = 0.

         if(verbose) print *,trim(CODE),': U/degrading map from nside ',&
              map%nside,' to ',nside_out
         call HealpixMap_udgrade(map, degraded_map, nside_out)

         !--------
         !Trim map
         !--------
         do nn = 1, degraded_map%nmaps
            where( degraded_mask%TQU(:,1) .lt. 1.)  degraded_map%TQU(:,1) = 0.
         end do
         
         !-------
         !Map out
         !-------
         outfile = trim(fileroot_out)//trim('lincom')//trim(fileextension_out)
         if(verbose) print *,CODE,': writing ',trim(outfile)      
         call HealpixMap_write(degraded_map,trim(outfile),units='')      

         !--------
         !Mask out
         !--------
         outfile = trim(fileroot_out)//trim('mask')//trim(fileextension_out)
         if(verbose) print *,CODE,': writing ',trim(outfile)      
         call HealpixMap_write(degraded_mask,trim(outfile),units='')      

      else
         !-------
         !Map out
         !-------
         outfile = trim(fileroot_out)//trim('lincom')//trim(fileextension_out)
         if(verbose) print *,CODE,': writing ',trim(outfile)      
         call HealpixMap_write(map,trim(outfile),units='')      

         !--------
         !Mask out
         !--------
         outfile = trim(fileroot_out)//trim('mask')//trim(fileextension_out)
         if(verbose) print *,CODE,': writing ',trim(outfile)      
         call HealpixMap_write(mask,trim(outfile),units='')      
      endif


      !--------------------
      !Calculate covariance
      !--------------------
      !Assume six element covariance
      call HealpixMap_init(covmat,nside2npix(map%nside),6,nested=nested) 
      
      
      !------------------------
      !Read in covariance files
      !------------------------
      if (covmat_infile(1) .ne. '') then
         if(verbose) print *,CODE,': Performing calculation of covariance '         
         do mm = 1, nmaps 
            if(verbose) print *,CODE,': Reading covmat ',trim(covmat_infile(mm))
            call HealpixMap_read(covmat_tmp,trim(covmat_infile(mm)))

            select case(covmat_format(mm))
            case("i")
               covmat%tqu(:,1) = covmat%tqu(:,1) + a(mm)**2 * covmat_tmp%tqu(:,covmat_firstcol(mm))
               nmap_covmat     = 1
               nmap_chi2       = 1
            case("rms_to_iqu")
               covmat%tqu(:,1) = covmat%tqu(:,1) + a(mm)**2 * covmat_tmp%tqu(:,covmat_firstcol(mm))**2
               covmat%tqu(:,4) = covmat%tqu(:,4) + a(mm)**2 * covmat_tmp%tqu(:,covmat_firstcol(mm))**2
               covmat%tqu(:,6) = covmat%tqu(:,6) + a(mm)**2 * covmat_tmp%tqu(:,covmat_firstcol(mm))**2
               nmap_covmat     = 1
               nmap_chi2       = 4
            case("iqu_uncorr")
               covmat%tqu(:,1) = covmat%tqu(:,1) + a(mm)**2 * covmat_tmp%tqu(:,covmat_firstcol(mm)+0) 
               covmat%tqu(:,4) = covmat%tqu(:,4) + a(mm)**2 * covmat_tmp%tqu(:,covmat_firstcol(mm)+1) 
               covmat%tqu(:,6) = covmat%tqu(:,6) + a(mm)**2 * covmat_tmp%tqu(:,covmat_firstcol(mm)+2) 
               nmap_covmat     = 3
               nmap_chi2       = 4
            case("iqu_corr")
               covmat%tqu(:,1) = covmat%tqu(:,1) + a(mm)**2 * covmat_tmp%tqu(:,covmat_firstcol(mm)+0) 
               covmat%tqu(:,2) = covmat%tqu(:,2) + a(mm)**2 * covmat_tmp%tqu(:,covmat_firstcol(mm)+1) 
               covmat%tqu(:,3) = covmat%tqu(:,3) + a(mm)**2 * covmat_tmp%tqu(:,covmat_firstcol(mm)+2) 
               covmat%tqu(:,4) = covmat%tqu(:,4) + a(mm)**2 * covmat_tmp%tqu(:,covmat_firstcol(mm)+3) 
               covmat%tqu(:,5) = covmat%tqu(:,5) + a(mm)**2 * covmat_tmp%tqu(:,covmat_firstcol(mm)+4) 
               covmat%tqu(:,6) = covmat%tqu(:,6) + a(mm)**2 * covmat_tmp%tqu(:,covmat_firstcol(mm)+5) 
               nmap_covmat     = 6
               nmap_chi2       = 4               
            case default
               stop 'Choose covmat_format = i, iqu_corr or iqu_uncorr'
            end select
         end do


      !------------------------------
      !Invert total covariance matrix
      !------------------------------
         if(verbose) print *,trim(CODE),': Inverting N to get N^-1 for accumulation'
         call HealpixMap_invert(covmat,covmat_firstcol=1)

      !---------------
      !Accumulate N^-1
      !---------------
         if(verbose) print *,trim(CODE),': Accumulating N^-1 from nside '//&
              trim(inttostr(map%nside))//' to '//trim(inttostr(nside_out))
         call HealpixMap_udgrade(covmat, degraded_covmat, nside_out)
         do mm = 1,nmap_covmat
            degraded_covmat%TQU(:,mm) = degraded_covmat%TQU(:,mm)*(float(map%nside)/float(nside_out))**2
         end do

      !-------------------------
      ! Calculate (dm).N^-1.(dm)
      !-------------------------
         if(verbose) print *,trim(CODE),': Calculating chi^2 map '
         call HealpixMap_init(chi2,degraded_map%npix,nmap_chi2,nested=nested) 

         select case(covmat_format(1)) !Based on covmat format of first map. Attention
         case("i")
            chi2%tqu(:,1) = degraded_map%tqu(:,1)*degraded_map%tqu(:,1)*degraded_covmat%tqu(:,1)
         case("iqu_uncorr")
            chi2%tqu(:,1) = degraded_map%tqu(:,1)*degraded_map%tqu(:,1)*degraded_covmat%tqu(:,1)
            chi2%tqu(:,2) = degraded_map%tqu(:,2)*degraded_map%tqu(:,4)*degraded_covmat%tqu(:,4)
            chi2%tqu(:,3) = degraded_map%tqu(:,3)*degraded_map%tqu(:,6)*degraded_covmat%tqu(:,6)
            chi2%tqu(:,4) = chi2%tqu(:,1) + chi2%tqu(:,2) + chi2%tqu(:,3)
         case("iqu_corr")
            call healpixmap_gausslike(degraded_map,degraded_covmat,chi2)
         case("rms_to_iqu")
            call healpixmap_gausslike(degraded_map,degraded_covmat,chi2)
         case default
            stop 'Choose covmat_format = i, iqu_corr or iqu_uncorr'
         end select
      end if

      !---------------
      ! Write Chi2 map
      !---------------
      outfile = trim(fileroot_out)//trim('chi2')//trim(fileextension_out)
      if(verbose) print *,CODE,': writing chi2 map ',trim(outfile)      
      call HealpixMap_write(chi2,trim(outfile),units='')      


      !---------------------------------------------------------
      ! Calculation of variance of linear combination (optional)
      !---------------------------------------------------------
!       if (covmat_infile .ne. '') then
!          if(verbose) print *,CODE,': Performing calculation of covariance '

!          if(verbose) print *,CODE,': Reading ',trim(covmat_infile)
!          call HealpixMap_read(covmat,trim(covmat_infile))

!          allocate(covmat_pix(nmaps,nmaps))

!          do pp = 0, covmat%npix - 1            
!             index = 0
!             do mm = 1, nmaps
!                do nn = mm, nmaps
!                   index = index + 1
!                   covmat_pix(mm,nn) = covmat%tqu(pp,index)
!                   covmat_pix(nn,mm) = covmat%tqu(pp,index)
!                end do
!             end do

!             variance_pix  =  dot_product(matmul(a,covmat_pix),a)
!             map%tqu(pp,1) = variance_pix
!          end do
!          map%nmaps = 1 ! Handles intensity only case for now - no polarization

         
!          if (nside_out .ne. map%nside) then
!             if(verbose) print *,trim(CODE),': Inverting variance map and accumulating from nside ',&
!                  map%nside,' to ',nside_out

!             !U/degrade inverse variance and convert back to a variance
!             do pp = 0, covmat%npix - 1            
!                map%tqu(pp,1) = 1./map%tqu(pp,1)            
!             end do
!             call HealpixMap_udgrade(map, degraded_map, nside_out)
!             do pp = 0, degraded_map%npix - 1            
!                degraded_map%tqu(pp,1) = 1./degraded_map%tqu(pp,1)*(float(nside_out)/float(map%nside))**2
!             end do
!             if(verbose) print *,CODE,': writing variance map ',trim(variance_outfile)      
!             call HealpixMap_write(degraded_map,trim(variance_outfile),units='')                  
!          else
!             if(verbose) print *,CODE,': writing variance map ',trim(variance_outfile)      
!             call HealpixMap_write(map,trim(variance_outfile),units='')      
!          endif
!            
!      endif


    end program jackknife
