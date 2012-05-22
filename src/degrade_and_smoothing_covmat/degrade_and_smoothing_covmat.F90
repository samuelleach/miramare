program degrade_and_smoothing_covmat

      use IniFile
      use mpi_stop
      use HealpixObj
      use HealpixVis
      use HealpixTools
      use pix_tools
      use beamtools
      use alm_tools
      
      character(LEN=Ini_max_string_len) InputFile
      logical bad,want_pol

      type(HealpixMap)   :: map,crop_mask
      Type(HealpixAlm)   :: alm
      type(HealpixMap)   :: map_tmp,map_tmp2
      Type(HealpixInfo)  :: HH
      type(HealpixMap)   :: map_out,crop_mask_out

      real undef, medianval      
      real(dp), pointer  :: beam(:)

      !----------
      !Parameters
      !----------
      logical verbose,infile_is_covmat,outfile_is_covmat
      character(LEN=256) :: infile
      character(LEN=256) :: outfile
      character(LEN=256) :: outfile_median
      character(LEN=256) :: outfile_rms
      character(LEN=256) :: variance_beam_file
      character(LEN=256) :: ordering_out
      logical discard_offdiagonal
      real(dp) rescale_factor,factor
      integer nlmax,nlmax_beam,nlmax_parfile,nside_out,nside_in
      integer covmat_firstcol
      real(dp) missing_pixel_value, variance_in
      integer extno,fits_hdu

      character(len=*), parameter :: CODE = "DEGRADE_AND_SMOOTHING_COVMAT"
      integer mm,nelement_in,nelement_out,offset,ordering_in,pp
      logical do_smoothing,do_accumulation,covmat_is_threebythree,crop_accumulated_map
      logical crop_map
      logical is_nested,inhomogeneous_noise
      integer ord_in,ord_out
      
#ifdef MPIPIX
      integer mpi_division_method
      integer i
      call mpi_init(i)
      mpi_division_method = 3
#endif
      
      undef                = -1.63750e30
      crop_accumulated_map = .false.
      crop_map             = .false.
      is_nested            = .false.
      !----------------------
      !Read in parameter file
      !----------------------
      InputFile = GetParam(1)
      call Ini_Open(InputFile, 1, bad, .false.)
      if (bad) then
         print *,CODE,': Accumulating and smoothing of N or N^-1 maps.'
         print *,'---------------------'
         print *,'Parameter file format'
         print *,'---------------------'
         print *,'verbose  = T or F (default = F).'
         print *,'infile   = Healpix map filename.'
         print *,'fits_hdu = FITS HDU (default = 2)'
         print *,'outfile  = Output Healpix map filename.'
         print *,'nlmax   = Maximum l for smoothing. (default = 3*nside_in+1)'
         print *,'outfile_median = Output Healpix map filename (map containing median value of outfile).'
         print *,'outfile_rms    = Output Healpix map filename containing rms noise level (only for diagonal noise).'
         print *,'infile_is_covmat  = T or F (default = F => Input is N^-1)'
         print *,'outfile_is_covmat = T or F (default = F => Output is N^-1)'
         print *,'covmat_firstcol   = First column of covmat data (default = 1)'
         print *,'nside_out         = Output nside'
         print *,'discard_offdiagonal = T or F (whether to drop the off-diagonal terms'
         print *,'                        of the output object.'
         print *,'rescale_factor      = Multiplicative factor for N (e.g. sigma^2) (default = 1.) '
         print *,'variance_beam_file  = Name of variance B_l file (output from make_window).'
         print *,'missing_pixel_value = Value for indicating missing pixels in map.'
         print *,'ordering_out = nest or ring (default = ordering of input map)'
         print *,' Instead of infile use:'
         print *,'nside_in    = nside of input homogeneous variance map.'
         print *,'variance_in = Variance of input homogeneous variance map.'
         print *,'---------------------'
         call DoStop(trim(concat(CODE,': Error opening parameter file')))
      endif
      Ini_fail_on_not_found = .false.
      if (InputFile == '') call DoStop('No parameter input file')      

      verbose             = Ini_Read_logical('verbose',.false.)
      infile              = Ini_Read_String('infile')
      fits_hdu            = Ini_Read_int('fits_hdu',2)
      outfile             = Ini_Read_String('outfile')
      nlmax_parfile       = Ini_Read_int('nlmax',8000)
      outfile_median      = Ini_Read_String('outfile_median')
      outfile_rms         = Ini_Read_String('outfile_rms')
      nside_out           = Ini_Read_int('nside_out',-1)
      covmat_firstcol     = Ini_Read_int('covmat_firstcol',1)
      variance_beam_file  = Ini_Read_String('variance_beam_file')
      rescale_factor      = Ini_Read_Real('rescale_factor',1.)
      infile_is_covmat    = Ini_Read_logical('infile_is_covmat',.false.)
      outfile_is_covmat   = Ini_Read_logical('outfile_is_covmat',.false.)
      discard_offdiagonal = Ini_Read_logical('discard_offdiagonal',.false.)
      missing_pixel_value = Ini_Read_Real('missing_pixel_value',undef)
      ordering_out        = Ini_Read_String('ordering_out')
      nside_in            = Ini_Read_int('nside_in',-1)
      variance_in         = Ini_Read_Real('variance_in',-1.)

      if(verbose) print *,trim(CODE),': Running'

      !-------------------------------
      ! Parameter defaulting for nlmax
      !-------------------------------
      if (nlmax_parfile .eq. 8000) then
         nlmax = 3*nside_out + 1
         if(verbose) print *,trim(CODE),': Setting default nlmax to '//trim(inttostr(nlmax))
      else
         nlmax = nlmax_parfile
         if(verbose) print *,trim(CODE),': Setting default nlmax to value in parfile '//trim(inttostr(nlmax))
      endif

      !---------------------
      !Read in variance beam
      !---------------------
      if (variance_beam_file .ne. '') then
         if(verbose) print *,trim(CODE),': Reading ',trim(variance_beam_file)
         call ReadBeam(beam, variance_beam_file, nlmax_beam)
         do_smoothing = .true.
         if(nlmax_beam .lt. nlmax) then
            nlmax = nlmax_beam
            print *,trim(CODE),': using nlmax from signal beam file: '//trim(inttostr(nlmax))
         endif
         if(nlmax_parfile .lt. nlmax) then
            nlmax = nlmax_parfile
            print *,trim(CODE),': using nlmax from parameter file: '//trim(inttostr(nlmax))
         endif
      else
         do_smoothing = .false.
         if(verbose) print *,trim(CODE),': No variance beam information supplied so no smoothing will be applied.'         
      endif

      !-------------------------------------------------------
      ! Read in covmat file or set up homogeneous variance map
      !-------------------------------------------------------
      if(infile .ne. '') then
         extno = fits_hdu - 2 
         if(verbose) print *,trim(CODE)//': Reading '//trim(infile)
         call HealpixMap_read(map,trim(infile),extno=extno)
         if (any(map%TQU(:,1) .eq. missing_pixel_value)) then
            if(verbose) print *,trim(CODE),': There are missing pixel in '//trim(infile)
            crop_accumulated_map = .true.
            where( abs(map%TQU(:,:)/missing_pixel_value -1.) .lt. 1.e-5 )  map%TQU(:,:) = 0.
         endif
         inhomogeneous_noise = .true.
         ord_in              = map%ordering

         !Set undef pixels to zero
         where( abs(map%TQU(:,:)/undef -1.) .lt. 1.e-5 )  map%TQU(:,:) = 0.
      else
         if(verbose) print *,trim(CODE),': Setting up homogeneous variance map'
         call HealpixMap_init(map_out,nside2npix(nside_out),nmaps = 1,nested=is_nested)
         if (do_smoothing) then
            factor = beam(0)
         else
            factor = 1.
         endif
         map_out%tqu(:,1)    = variance_in * (float(nside_out)/float(nside_in))**2*factor
         inhomogeneous_noise = .false.
      endif

      if(inhomogeneous_noise) then

!      if (map%nside .gt. nside_out) then
      if (map%nside .ne. nside_out) then
         do_accumulation = .true.
      else
         do_accumulation = .false.
      end if
      want_pol = .false.


      select case(ordering_out)
      case("")
         ord_out = ord_in
      case("nest")
         ord_out = ord_nest
      case("ring")
         ord_out = ord_ring
      case default
         stop 'Choose ordering_out = nest or ring'
      end select


      !------------------------------------
      ! Evaluate number of colums of covmat
      !------------------------------------
      offset      = covmat_firstcol - 1
      nelement_in = map%nmaps - offset
      select case(nelement_in)
      case(6:7)
         covmat_is_threebythree = .true.         
         if(discard_offdiagonal) then
            nelement_out = 3
         else
            nelement_out = 6
         endif
      case(3)
         covmat_is_threebythree = .false.
         nelement_out           = nelement_in
      case(2)
         covmat_is_threebythree = .false.
         nelement_out           = 1 ! Drop second column
      case(1)
         covmat_is_threebythree = .false.
         nelement_out           = 1 
      end select

      !------------------------
      !Invert covariance matrix
      !------------------------
      if(infile_is_covmat .and. do_accumulation) then
         if(verbose) print *,trim(CODE),': Inverting N to get N^-1 for accumulation'
         call HealpixMap_invert(map,covmat_firstcol=covmat_firstcol)
      endif
      
      !--------------------------------------------------------
      !ACCUMULATE N^-1, and optionally crop away missing pixels
      !--------------------------------------------------------
      if(do_accumulation) then
         if(crop_accumulated_map) then
            if(verbose) print *,trim(CODE),': Allocating memory for cropping mask.'
            if (map%ordering .eq. ord_nest) is_nested = .true.
            call HealpixMap_init(crop_mask,map%npix,nmaps = 1,nested=is_nested)
            crop_mask%tqu(:,1) = 1.
            do pp = 0,map%npix-1
               if( map%TQU(pp,covmat_firstcol) .eq. 0. )  crop_mask%TQU(pp,1) = 0.
            end do
            
            if(verbose) print *,trim(CODE),': Degrading mask from nside '//&
                 trim(inttostr(map%nside))//' to '//trim(inttostr(nside_out))
            call HealpixMap_udgrade(crop_mask, crop_mask_out, nside_out)
            where( crop_mask_out%TQU(:,1) .lt. 1. ) crop_mask_out%TQU(:,1) = 0.
         endif
         
         if(verbose) print *,trim(CODE),': Accumulating N^-1 from nside '//&
              trim(inttostr(map%nside))//' to '//trim(inttostr(nside_out))
         call HealpixMap_udgrade(map, map_tmp2, nside_out)
         do mm = covmat_firstcol,covmat_firstcol+nelement_in-1
            map_tmp2%TQU(:,mm) = map_tmp2%TQU(:,mm)*(float(map%nside)/float(nside_out))**2
         end do
!         call HealpixMap_write(map_tmp2,'!testdegrade1.fits',units='') !SL containt zero weight elements
         
         if(crop_accumulated_map) then
            if(verbose) print *,trim(CODE),': Cropping degraded map.'
            do mm = covmat_firstcol, covmat_firstcol+nelement_in-1
               do pp = 0,crop_mask_out%npix - 1
                  if(crop_mask_out%TQU(pp,1)  .lt. 1.) map_tmp2%TQU(pp,mm) = 0.
               end do
            enddo
         endif
         map = map_tmp2
!         call HealpixMap_write(map,'!testdegrade.fits',units='')

      else
         if (map%ordering .eq. ord_nest) is_nested = .true.
         call HealpixMap_init(crop_mask_out,map%npix,nmaps = 1,nested=is_nested)
         crop_mask_out%tqu(:,1) = 1.
         do pp = 0,map%npix - 1
            if( map%TQU(pp,covmat_firstcol) .eq. 0. )  then
               crop_mask_out%TQU(pp,1) = 0.
               crop_map = .true.
            endif
         end do
      endif


      !--------------------------------------
      !Invert matrix to get covariance matrix
      !--------------------------------------
      if(do_accumulation .or.  .not. infile_is_covmat) then
         if(verbose) print *,trim(CODE),': Inverting N^-1 to get covariance matrix'
         call HealpixMap_invert(map,covmat_firstcol=covmat_firstcol)
      endif
!      call HealpixMap_write(map,'!testinvert1.fits',units='')
 
      call HealpixMap_init(map_tmp,map%npix,nmaps=1)
      if(verbose) print *,trim(CODE),': Initialising MPI healpix'
      call HealpixInit(HH,map%nside, 2*nlmax,want_pol, w8dir='')
      if (HH%MpiID .eq. 0) then
      !----------------------------------------------------------------
      !Do the variance smoothing of each entry in the covariance matrix
      !----------------------------------------------------------------
         if (do_smoothing) then
            ordering_in = map%ordering
            if (ordering_in .eq. ord_nest) then
               if(verbose) print *,trim(CODE),': Forcing nest to ring ordering.'
               call HealpixMap_forceRing(map)
               call HealpixMap_forceRing(crop_mask_out)
            end if

            do mm = 1, nelement_in
               if(verbose) print *,trim(CODE),': map2alm, map '//trim(inttostr(mm))//&
                    ' out of '//trim(inttostr(nelement_in))
               do pp = 0,map_tmp%npix-1
                  map_tmp%tqu(pp,1) = map%tqu(pp,offset+mm)
               end do
               call HealpixMap2alm(HH,map_tmp,alm,nlmax)
               
               if(verbose) print *,trim(CODE),': Smoothing alm '//trim(inttostr(mm))
               call HealpixAlm_smooth_beam(alm,beam)
               
               if(verbose) print *,trim(CODE),': alm2map       '//trim(inttostr(mm))
               call HealpixAlm2Map(HH,alm, map_tmp, map%npix)

               do pp=0,crop_mask_out%npix - 1
                  if ( crop_mask_out%TQU(pp,1) .lt. 1. ) map_tmp%TQU(pp,1) = 0.
               enddo
               do pp=0,map%npix - 1
                  map%tqu(pp,offset+mm) = map_tmp%tqu(pp,1)
               end do
            end do
            if (ord_out .ne. ord_ring) then
               if(verbose) print *,trim(CODE),': Forcing ring back to nest ordering.'
               call HealpixMap_forceNest(map)
            end if
         end if

         !-----------------------------------------------------------
         !Discard off-diagonal terms of covariance matrix, if desired
         !-----------------------------------------------------------
         call HealpixMap_init(map_out,nside2npix(nside_out),nmaps = nelement_out) !Buggy ?
         map_out%ordering = map%ordering
         if(discard_offdiagonal .and. nelement_in .ge. 6) then
            if(verbose) print *,trim(CODE),': Discarding off diagonal terms.'
               do pp=0,map_out%npix - 1
                  map_out%tqu(pp,1) = map%tqu(pp,offset+1)
                  map_out%tqu(pp,2) = map%tqu(pp,offset+4)
                  map_out%tqu(pp,3) = map%tqu(pp,offset+6)
               end do
         else
            do mm = 1,nelement_out
               do pp=0,map_out%npix - 1
                  map_out%tqu(pp,mm) = map%tqu(pp,offset+mm)
               end do
            end do
         endif
      endif
      endif

      if (HH%MpiID .eq. 0) then
         !-------------
         ! Rescale maps
         !-------------
         if(verbose) print *,trim(CODE),': rescaling covariance map by factor = ',rescale_factor
         do mm = 1, nelement_out
            do pp = 0,map_out%npix - 1
               map_out%tqu(pp,mm) = map_out%tqu(pp,mm) * rescale_factor 
            end do
         end do

         !-----------
         !Reorder map
         !-----------
         select case(ord_out)
         case(ord_ring)
            if(verbose) print *,trim(CODE)//': Forcing output map to ring ordering.'
            call HealpixMap_ForceRing(map_out)
            if(crop_accumulated_map .or. crop_map) then
               call HealpixMap_ForceRing(crop_mask_out)
            end if
         case(ord_nest)
            if(verbose) print *,trim(CODE)//': Forcing output map to nest ordering.'
            call HealpixMap_ForceNest(map_out)
            call HealpixMap_ForceNest(crop_mask_out)
         end select

         !--------------------------
         !Optionally output RMS maps
         !--------------------------
         if (outfile_rms .ne. '' .and. nelement_out .le. 3) then
            do mm=1, map_out%nmaps
               do pp = 0,map_out%npix - 1
                  map_out%tqu(pp,mm) = sqrt(map_out%tqu(pp,mm))
               end do               
               if(crop_accumulated_map .or. crop_map) then                  
                  do pp = 0,crop_mask_out%npix - 1
                     if ( crop_mask_out%TQU(pp,1) .lt. 1. ) map_out%TQU(pp,mm) = 0.
                  end do
               end if
            end do
            if(verbose) print *,trim(CODE)//': Writing rms noise map to ',trim(outfile_rms)
            call HealpixMap_write(map_out,trim(outfile_rms),units='')            
            do mm=1, map_out%nmaps
               do pp = 0,map_out%npix - 1
                  map_out%tqu(pp,mm) = map_out%tqu(pp,mm)**2
               end do
            end do
         endif
         
         !---------------------------------------
         !Do final inversion of matrix if desired
         !---------------------------------------
         if(.not. outfile_is_covmat) then
            if(verbose) print *,trim(CODE),': Inverting covariance matrix to get N^-1.'
            call HealpixMap_invert(map_out)
         endif

         !-----------------
         !Write map to disk
         !-----------------
         if(verbose) print *,trim(CODE),': Writing map ',trim(outfile)
         call HealpixMap_write(map_out,trim(outfile),units='')

         !------------------------
         !Write median map to disk
         !------------------------
         if (outfile_median .ne. '') then
            do mm=1, map_out%nmaps
               medianval = HealpixMap_getmedian(map_out,mapindex = mm)
               if(verbose) print *,trim(CODE)//': Map '//trim(inttostr(mm))//&
                    ' median value = ',medianval
               do pp = 0,map_out%npix - 1
                  map_out%tqu(pp,mm) = medianval
                end do
               if(crop_accumulated_map) then
                  do pp = 0,crop_mask_out%npix - 1
                     if ( crop_mask_out%TQU(pp,1) .lt. 1. ) then
                        map_out%TQU(pp,mm) = 0.
                     end if
                  end do
               endif
            end do
            if(verbose) print *,trim(CODE)//': Writing median map ',trim(outfile_median)
            call HealpixMap_write(map_out,trim(outfile_median),units='')
         endif

      endif
      
#ifdef MPIPIX
      call HealpixFree(HH)
      call mpi_finalize(i)
#endif
      
    end program degrade_and_smoothing_covmat

