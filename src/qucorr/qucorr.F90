program qucorr

      use IniFile
      use mpi_stop
      use HealpixObj
      use HealpixVis
      
      character(LEN=Ini_max_string_len) InputFile
      logical bad

      type(HealpixMap)   :: amap_q
      type(HealpixMap)   :: amap_u
      type(HealpixMap)   :: projectionmap
      type(HealpixMap)   :: correctionmap
      type(HealpixMap)   :: correctionvar
      type(HealpixMap)   :: leakage_correctionvar
      type(HealpixMap)   :: afactor_correctionvar
      type(HealpixMap)   :: leakagemap
      type(HealpixMap)   :: leakagevar
      type(HealpixMap)   :: tempmap

      real undef      
      integer pp,mm,ss,nn,cc,index
      
      !Parameters
      integer nhorn
      logical verbose
      character(LEN=100) :: leakagemapfile,leakagevarfile,outfilestring
      character(LEN=100) :: correctionmapfile,projectionmapfile
      character(LEN=100) :: correctionvarfile,leakage_correctionvarfile,afactor_correctionvarfile
      real(dp), allocatable :: afactor(:)
      real(dp), allocatable :: dafactor(:)
      character(256), allocatable :: amap_qfile(:),amap_ufile(:)


      integer nside_out

      !Data
!      real(dp), dimension(:,:), allocatable :: covmat_pix
!      real(dp) variance_pix
      
      
      character(len=*), parameter :: CODE = "QUCORR"

      undef= -1.6375e30
      !----------------------
      !Read in parameter file
      !----------------------
      InputFile = GetParam(1)
      call Ini_Open(InputFile, 1, bad, .false.)
      if (bad) then
         print *,CODE,': Make qu correction maps and propagate uncertainty'
         print *,'---------------------'
         print *,'Parameter file format'
         print *,'---------------------'
         print *,'verbose         = T or F (default = F).'
         print *,'nhorn           = number of horns to combine'
         print *,'correctionmap   = output correction map filename'
         print *,'correctionvar   = output correction map (total) variance filename'
         print *,'leakage_correctionvar   = output correction map variance filename (leakage only part)'
         print *,'afactor_correctionvar   = output correction map variance filename (afactor only part)'
         print *,'projectionmap   = output leakage projection map filename'
         print *,'amap_q(n)   = nth Healpix Q amap filename'
         print *,'amap_u(n)   = nth Healpix U amap filename'         
         print *,'afactor(n)  = nth bandpass a factor'
         print *,'dafactor(n) = nth bandpass a factor RMS error'
         print *,'leakagemap   = input leakage map filename'
         print *,'leakagevar   = input leakage map variance filename'
         print *,'---------------------'
         call DoStop(trim(concat(CODE,': Error opening parameter file')))
      endif
      Ini_fail_on_not_found = .false.
      if (InputFile == '') call DoStop('No parameter input file')      

      verbose          = Ini_Read_logical('verbose',.false.)
      nhorn            = Ini_Read_Int('nhorn',1)
      leakagemapfile   = Ini_Read_String('leakagemap')
      leakagevarfile   = Ini_Read_String('leakagevar')
      projectionmapfile= Ini_Read_String('projectionmap')
      correctionmapfile= Ini_Read_String('correctionmap')
      correctionvarfile= Ini_Read_String('correctionvar')
      leakage_correctionvarfile= Ini_Read_String('leakage_correctionvar')
      afactor_correctionvarfile= Ini_Read_String('afactor_correctionvar')
      allocate(amap_qfile(nhorn),amap_ufile(nhorn),afactor(nhorn),dafactor(nhorn))
      do mm = 1, nhorn
         amap_qfile(mm)       = Ini_Read_String_array('amap_q',mm)
         amap_ufile(mm)       = Ini_Read_String_array('amap_u',mm)
         afactor(mm)      = Ini_Read_Real_array('afactor',mm,0.)
         dafactor(mm)     = Ini_Read_Real_array('dafactor',mm,0.)
      end do

      if(verbose) print *,CODE,': running'
      do mm = 1, nhorn
         if(verbose) print *,CODE,': a, da = ',afactor(mm),dafactor(mm)
      end do

      !---------------------------------
      !Read in leakage map and get nside
      !---------------------------------
      if(verbose) print *,CODE,': reading ',trim(leakagemapfile)      
      call HealpixMap_read(leakagemap,trim(leakagemapfile))      
!      if (nside_out .eq. -1) then
      nside_out = leakagemap%nside
!      endif
      if(verbose) print *,CODE,': reading ',trim(leakagevarfile)      
      call HealpixMap_read(leakagevar,trim(leakagevarfile))      


      !--------------
      ! Read in Amaps
      !--------------
      call HealpixMap_init(amap_q,leakagemap%npix,nmaps=nhorn,nested=.true.)
      call HealpixMap_init(amap_u,leakagemap%npix,nmaps=nhorn,nested=.true.)
      do mm = 1, nhorn
         if(verbose) print *,CODE,': reading ',trim(amap_qfile(mm))      
         call healpixMap_read(tempmap,trim(amap_qfile(mm)),&
              ordering=leakagemap%ordering,nside=leakagemap%nside)
         amap_q%tqu(:,mm) = tempmap%tqu(:,1)

         if(verbose) print *,CODE,': reading ',trim(amap_ufile(mm))      
         call healpixMap_read(tempmap,trim(amap_ufile(mm)),&
              ordering=leakagemap%ordering,nside=leakagemap%nside)
         amap_u%tqu(:,mm) = tempmap%tqu(:,1)
      end do


      !-------------------------
      ! Calculate projection map
      !-------------------------
      call HealpixMap_init(projectionmap,leakagemap%npix,nmaps=3,nested=.true.)
      do mm = 1, nhorn
         projectionmap%tqu(:,2) = projectionmap%tqu(:,2) - afactor(mm)*amap_q%tqu(:,mm)
         projectionmap%tqu(:,3) = projectionmap%tqu(:,3) - afactor(mm)*amap_u%tqu(:,mm)         
      end do
      if(verbose) print *,CODE,': writing ',trim(projectionmapfile)      
      call HealpixMap_write(projectionmap,trim(projectionmapfile),units='')


      !-------------------------
      ! Calculate correction map
      !-------------------------
      call HealpixMap_init(correctionmap,leakagemap%npix,nmaps=3,nested=.true.)
      correctionmap%tqu(:,2) = leakagemap%tqu(:,1) * projectionmap%tqu(:,2)
      correctionmap%tqu(:,3) = leakagemap%tqu(:,1) * projectionmap%tqu(:,3)
      if(verbose) print *,CODE,': writing ',trim(correctionmapfile)      
      call HealpixMap_write(correctionmap,trim(correctionmapfile),units='')


      !------------------------------
      ! Calculate correction variance (leakage, afactor and sum)
      !------------------------------
      call HealpixMap_init(correctionvar,leakagemap%npix,nmaps=3,nested=.true.)
      call HealpixMap_init(leakage_correctionvar,leakagemap%npix,nmaps=3,nested=.true.)
      call HealpixMap_init(afactor_correctionvar,leakagemap%npix,nmaps=3,nested=.true.)

      !leakage uncertainty
      leakage_correctionvar%tqu(:,2) = projectionmap%tqu(:,2)**2 * leakagevar%tqu(:,1)
      leakage_correctionvar%tqu(:,3) = projectionmap%tqu(:,3)**2 * leakagevar%tqu(:,1)
      if(verbose) print *,CODE,': writing ',trim(leakage_correctionvarfile)      
      call HealpixMap_write(leakage_correctionvar,trim(leakage_correctionvarfile),units='')

      !afactor uncertainty
      afactor_correctionvar%tqu(:,2) = dafactor(1)**2 * amap_q%tqu(:,1)**2 * leakagemap%tqu(:,1)**2
      afactor_correctionvar%tqu(:,3) = dafactor(1)**2 * amap_u%tqu(:,1)**2 * leakagemap%tqu(:,1)**2
      do mm = 2, nhorn
         afactor_correctionvar%tqu(:,2) = afactor_correctionvar%tqu(:,2) + &
              dafactor(mm)**2 * amap_q%tqu(:,mm)**2 * leakagemap%tqu(:,1)**2
         afactor_correctionvar%tqu(:,3) = afactor_correctionvar%tqu(:,3) + &
              dafactor(mm)**2 * amap_u%tqu(:,mm)**2 * leakagemap%tqu(:,1)**2
      end do
      if(verbose) print *,CODE,': writing ',trim(afactor_correctionvarfile)      
      call HealpixMap_write(afactor_correctionvar,trim(afactor_correctionvarfile),units='')

      !sum
      correctionvar%tqu(:,2) = afactor_correctionvar%tqu(:,2) + leakage_correctionvar%tqu(:,2)
      correctionvar%tqu(:,3) = afactor_correctionvar%tqu(:,3) + leakage_correctionvar%tqu(:,3)
      if(verbose) print *,CODE,': writing ',trim(correctionvarfile)      
      call HealpixMap_write(correctionvar,trim(correctionvarfile),units='')



!      correctionvar%tqu(:,2) = projectionmap%tqu(:,2)**2 * leakagevar%tqu(:,1)
!      correctionvar%tqu(:,3) = projectionmap%tqu(:,3)**2 * leakagevar%tqu(:,1)
!      do mm = 1, nhorn
!         correctionvar%tqu(:,2) = correctionvar%tqu(:,2) + &
!              dafactor(mm)**2 * amap_q%tqu(:,mm)**2 * leakagemap%tqu(:,1)**2
!         correctionvar%tqu(:,3) = correctionvar%tqu(:,3) + &
!              dafactor(mm)**2 * amap_u%tqu(:,mm)**2 * leakagemap%tqu(:,1)**2
!      end do
!      if(verbose) print *,CODE,': writing ',trim(correctionvarfile)      
!      call HealpixMap_write(correctionvar,trim(correctionvarfile),units='')



 !      !----------------------------
!       !Read in other maps and coadd
!       !----------------------------
!       do mm = 2, nmaps
!         if(verbose) print *,CODE,': reading ',trim(infile(mm))      
!         call HealpixMap_read(map_temp,trim(infile(mm)))

!         if (map_temp%ordering .ne. map%ordering) then
!         !Force the ordering of the maps to be the same         
!            if (map%ordering .eq. ord_ring) then
!               if(verbose) print *,CODE,': Forcing to ring order map ',trim(infile(mm))      
!               call HealpixMap_forcering(map_temp)
!            else
!               if(verbose) print *,CODE,': Forcing to nest order map ',trim(infile(mm))      
!               call HealpixMap_forcenest(map_temp)
!            endif                   
!         end if

!         do cc = 1, map_temp%nmaps
!            where( map_temp%TQU(:,cc) .eq. undef)  map_temp%TQU(:,cc)= 0.
!         end do


!         map_temp%TQU      = map_temp%TQU*a(mm) 
!         !--------------------------
!         ! Add offsets only to T map
!         !--------------------------
!         do pp = 0, map%npix-1
!            map_temp%TQU(pp,1) = map_temp%TQU(pp,1) + c(mm) 
!         end do

!         do nn = 1, map_temp%nmaps
!            if (nn .gt. map%nmaps) cycle
!            do pp = 0, map%npix-1
!               map%tqu(pp,nn) = map%tqu(pp,nn) + map_temp%tqu(pp,nn)
!            end do
!         end do
!       enddo

!       !---------------------------------------
!       ! Optionally set temperature map to zero
!       !---------------------------------------
!       if(remove_temperature) then
!          map%tqu(:,1) = 0.
!       end if
      

!       !-------------------------
!       ! U/degrade map and output
!       !-------------------------
!       if (nside_out .ne. map%nside) then
!          if(verbose) print *,trim(CODE),': U/degrading map from nside ',&
!               map%nside,' to ',nside_out
!          call HealpixMap_udgrade(map, degraded_map, nside_out)
!          if(verbose) print *,CODE,': writing ',trim(outfile)      
!          call HealpixMap_write(degraded_map,trim(outfile),units='')      
!       else
!          if(verbose) print *,CODE,': writing ',trim(outfile)      
!          call HealpixMap_write(map,trim(outfile),units='')      
!       endif

!       !---------------------------------------------------------
!       ! Calculation of variance of linear combination (optional)
!       !---------------------------------------------------------
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
            
!       endif


    end program qucorr
