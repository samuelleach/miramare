program degrade_map

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
      type(HealpixMap)   :: degraded_map
      type(HealpixMap)   :: degraded_hitmap

      real undef      
      integer pp,mm,index
      real(dp)           :: a(3,3)!,ainv(3,3)



      !Parameters
      logical verbose
      character(LEN=100) :: infile
      character(LEN=100) :: outfile
      character(LEN=100) :: outweightfile
      character(LEN=100) :: weightfile
      integer nside,weightformat
      logical hitmap_is_iqu_format

      character(len=*), parameter :: CODE = "DEGRADE_MAP"
      
      
      undef= -1.63750e30
      !----------------------
      !Read in parameter file
      !----------------------
      InputFile = GetParam(1)
      call Ini_Open(InputFile, 1, bad, .false.)
      if (bad) then
         print *,CODE,': Degrading of healpix maps via N^-1 or hitmap weighting.'
         print *,'---------------------'
         print *,'Parameter file format'
         print *,'---------------------'
         print *,'verbose = T or F (default = F).'
         print *,'infile = Healpix map filename.'
         print *,'outfile = Output Healpix map filename.'
         print *,'outweightfile = Output degraded weight filename (optional).'
         print *,'nside = Output Healpix map nside value.'
         print *,'weightformat = format of weighting file (default = 1).'
         print *,' weightformat = 1 : Weighting via hitmap or 3 col IQU weight in weightfile.'
         print *,' weightformat = 2 : Weighting via hitmap in 2nd or 4th column of an'
         print *,'                    I or IQU map respectively.'
         print *,' weightformat = 3 : Weighting via (3 x 3) madam N^-1 map.'
         print *,' weightformat = 4 : Degrading of WMAP IQU(N,NQQ,NQU,NUU).'
         print *,'weightfile = Healpix weight map filename.'
         print *,'---------------------'
         call DoStop(trim(concat(CODE,': Error opening parameter file')))
      endif
      Ini_fail_on_not_found = .false.
      if (InputFile == '') call DoStop('No parameter input file')      

      verbose = Ini_Read_logical('verbose',.false.)
      infile = Ini_Read_String('infile')
      outfile = Ini_Read_String('outfile')
      outweightfile = Ini_Read_String('outweightfile')
      weightfile = Ini_Read_String('weightfile')
      nside = Ini_Read_Int('nside',-1)
      weightformat = Ini_Read_Int('weightformat',1)

      if(verbose) print *,trim(CODE),': running'

      if(verbose) print *,trim(CODE),': reading ',trim(infile)
      call HealpixMap_read(map,trim(infile))


      select case(weightformat)
      !-----------------------------------------------------------------------------------------
      case(1)
         if(verbose) print *,trim(CODE),': weightformat = 1 : Weighting via hitmap or 3 col IQU weight in weightfile.'
      !-----------------------------------------------------------------------------------------

         if (map%nmaps .eq. 2 .or. map%nmaps .eq. 4) then
            print *,trim(CODE), ': Warning: expected input map to have 1 or 3 columns.'
            stop
         end if

         if(verbose) print *,trim(CODE),': reading ',trim(weightfile)
         call HealpixMap_read(hitmap,trim(weightfile))
         if (hitmap%nside .ne. map%nside) then
            call DoStop(trim(concat(CODE,': signal and weight map nside is not equal.')))            
         end if

         if (map%nmaps .eq. hitmap%nmaps) then
            hitmap_is_iqu_format = .true.
         else		 
            hitmap_is_iqu_format = .false.
         endif

         !Apply weighting to map
         do mm=1, map%nmaps
            if(hitmap_is_iqu_format) then
               index = mm
	    else
               index = 1
            endif
            map%TQU(:,mm)= map%TQU(:,mm)*hitmap%tqu(:,index)
         end do

         if(verbose) print *,trim(CODE),': Degrading map from nside ',&
              hitmap%nside,' to ',nside
         call HealpixMap_udgrade(map, degraded_map, nside)
         if(verbose) print *,trim(CODE),': Accumulating hitmap from nside ',&
              hitmap%nside,' to ',nside
         call HealpixMap_udgrade(hitmap, degraded_hitmap, nside)

         do mm=1, map%nmaps
            if(hitmap_is_iqu_format) then
               index = mm
	    else
               index = 1
            endif
            degraded_map%tqu(:,mm)=degraded_map%tqu(:,mm)/degraded_hitmap%tqu(:,index)
         end do

         degraded_hitmap%TQU(:,:) = degraded_hitmap%TQU(:,:)*&
              (float(hitmap%nside)/float(nside))**2

         if(outweightfile .ne. '') then
            if(verbose) print *,trim(CODE),': writing ',trim(outweightfile)
            call HealpixMap_write(degraded_hitmap,trim(outweightfile),units='')
         end if
      !-----------------------------------------------------------------------------------------
      case(2)
      !-----------------------------------------------------------------------------------------
         if(verbose) print *,trim(CODE),': weightformat = 2 : Weighting via hitmap in 2nd or 4th column of an'
         if(verbose) print *,trim(CODE),':                    I or IQU map respectively.'

         if (map%nmaps .ne. 2 .and. map%nmaps .ne. 4) then
            print *,map%nmaps
            print *,trim(CODE), ': Warning: expected input map to have 2 or 4 columns.'
            stop
         end if

         !Apply weighting to map
         do mm = 1, map%nmaps-1
            map%TQU(:,mm)= map%TQU(:,mm)*map%tqu(:,map%nmaps)
         end do

         if(verbose) print *,trim(CODE),': degrading map from nside ',&
              map%nside,' to ',nside
         call HealpixMap_udgrade(map, degraded_map, nside)

         do mm =1, map%nmaps-1
            degraded_map%tqu(:,mm)=degraded_map%tqu(:,mm)/degraded_map%tqu(:,map%nmaps)
         end do

         degraded_map%TQU(:,map%nmaps) = degraded_map%TQU(:,map%nmaps)*&
              (float(map%nside)/float(nside))**2

         if(outweightfile .ne. '') then
            if(verbose) print *,trim(CODE),&
                 ': WARNING : not expecting value for outweightfile in parameter file.'
            stop
         end if
      !-----------------------------------------------------------------------------------------
      case(3)
      !-----------------------------------------------------------------------------------------
         if(verbose) print *,trim(CODE),&
              ': weightformat = 3 : Weighting via (3 x 3) madam N^-1 map.'

         if(verbose) print *,trim(CODE),': Reading ',trim(weightfile)
         call HealpixMap_read(hitmap,trim(weightfile))
         if (hitmap%nside .ne. map%nside) then
            call DoStop(trim(concat(CODE,': Signal and weight map nside is not equal.')))            
         end if

         if (hitmap%nmaps .ne. 6) then
            call DoStop(trim(concat(CODE,': Expecting a 6 column Healpix map.')))            
            stop
         end if

         if (hitmap%ordering .ne. map%ordering) then
            select case(map%ordering)
            case(ord_nest)
               print *,trim(CODE),': Warning: reordering hitmap to nested.'
               call HealpixMap_ForceNest(hitmap)
            case(ord_ring)
               print *,trim(CODE),': Warning: reordering hitmap to ring format.'
               call HealpixMap_ForceRing(hitmap)
            end select
         end if

         if(verbose) print *,trim(CODE),': Multiplying N^-1 x map.'
         do pp = 0, map%npix-1 
            a(1,1) = hitmap%tqu(pp,1)
            a(1,2) = hitmap%tqu(pp,2)
            a(1,3) = hitmap%tqu(pp,3)
            a(2,1) = hitmap%tqu(pp,2)
            a(2,2) = hitmap%tqu(pp,4)
            a(2,3) = hitmap%tqu(pp,5)
            a(3,1) = hitmap%tqu(pp,3)
            a(3,2) = hitmap%tqu(pp,5)
            a(3,3) = hitmap%tqu(pp,6)
            map%tqu(pp,1:3)= matmul(a,map%tqu(pp,1:3))
         end do

         if(verbose) print *,trim(CODE),': Degrading map from nside ',&
              map%nside,' to ',nside
         call HealpixMap_udgrade(map, degraded_map, nside)

         if(verbose) print *,trim(CODE),': Accumulating pixel matrix from nside ',&
              hitmap%nside,' to ',nside
         call HealpixMap_udgrade(hitmap, degraded_hitmap, nside)


         !NEED TO CHECK THIS PART OF THE ALGORITHM
         if(verbose) print *,trim(CODE),': Rescaling signal map'
         degraded_map%tqu(:,1)= degraded_map%tqu(:,1)/degraded_hitmap%tqu(:,1)
         degraded_map%tqu(:,2)= degraded_map%tqu(:,2)/degraded_hitmap%tqu(:,4)
         degraded_map%tqu(:,3)= degraded_map%tqu(:,3)/degraded_hitmap%tqu(:,6)

         degraded_hitmap%TQU(:,:) = degraded_hitmap%TQU(:,:)*&
              (float(hitmap%nside)/float(nside))**2

         if(outweightfile .ne. '') then
            if(verbose) print *,trim(CODE),': writing ',trim(outweightfile)
            call HealpixMap_write(degraded_hitmap,trim(outweightfile),units='')
         end if

      !-----------------------------------------------------------------------------------------
      case(4)
      !-----------------------------------------------------------------------------------------
         if(verbose) print *,trim(CODE),&
              ': weightformat = 4 : Degrading of WMAP IQU(N,NQQ,NQU,NUU).'

         if(verbose) print *,trim(CODE),': reading hitcounts from ',trim(infile)
         call HealpixMap_read(hitmap,trim(infile),extno=1)

         if (hitmap%nmaps .ne. 4) then
            print *,map%nmaps
            call DoStop(trim(concat(CODE,': Expecting 4 columns in extension 1.')))            
            stop
         end if

         if(verbose) print *,trim(CODE),': Multiplying N^-1 x map.'
         do pp = 0, map%npix-1 
            a(1,1) = hitmap%tqu(pp,1)
            a(1,2) = 0.
            a(1,3) = 0.
            a(2,1) = 0.
            a(2,2) = hitmap%tqu(pp,2)
            a(2,3) = hitmap%tqu(pp,3)
            a(3,1) = 0.
            a(3,2) = hitmap%tqu(pp,3)
            a(3,3) = hitmap%tqu(pp,4)
            map%tqu(pp,1:3)= matmul(a,map%tqu(pp,1:3))
         end do

         if(verbose) print *,trim(CODE),': Degrading map from nside ',&
              map%nside,' to ',nside
         call HealpixMap_udgrade(map, degraded_map, nside)

         if(verbose) print *,trim(CODE),': Accumulating pixel matrix from nside ',&
              hitmap%nside,' to ',nside
         call HealpixMap_udgrade(hitmap, degraded_hitmap, nside)

         !NEED TO CHECK THIS PART OF THE ALGORITHM
         if(verbose) print *,trim(CODE),': Rescaling signal map'
         degraded_map%tqu(:,1)= degraded_map%tqu(:,1)/degraded_hitmap%tqu(:,1)
         degraded_map%tqu(:,2)= degraded_map%tqu(:,2)/degraded_hitmap%tqu(:,2)
         degraded_map%tqu(:,3)= degraded_map%tqu(:,3)/degraded_hitmap%tqu(:,4)

         !Fourth column of WMAP extension 0 is Nobs
         degraded_map%TQU(:,4) = degraded_map%TQU(:,4)*&
              (float(hitmap%nside)/float(nside))**2

         degraded_hitmap%TQU(:,:) = degraded_hitmap%TQU(:,:)*&
              (float(hitmap%nside)/float(nside))**2

         if(outweightfile .ne. '') then
            if(verbose) print *,trim(CODE),': writing ',trim(outweightfile)
            call HealpixMap_write(degraded_hitmap,trim(outweightfile),units='')
         end if
      case default
         print *,trim(CODE), 'Not implemented this weightformat.'
         stop
      end select

      if(verbose) print *,trim(CODE),': writing ',trim(outfile)
      call HealpixMap_write(degraded_map,trim(outfile),units='')

      

    end program degrade_map
