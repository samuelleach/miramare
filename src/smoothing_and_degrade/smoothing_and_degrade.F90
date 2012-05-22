program smoothing_and_degrade

      use IniFile
      use mpi_stop
      use HealpixObj
      use HealpixVis
      use pix_tools
      use beamtools
      use alm_tools
      use mask_maker
      
      character(LEN=Ini_max_string_len) InputFile
      logical bad,want_pol

      type(HealpixMap)   :: map,map_smooth
      type(HealpixMap)   :: mask,inmask
      Type(HealpixAlm)   :: alm
      type(HealpixMap)   :: map_tmp
      Type(HealpixInfo)  :: HH,HH2
      type(HealpixMap)   :: map_out

      real undef      
      integer npix,pp,i,ff
      integer*8 mask_npix
      real(dp), pointer     :: beam(:)
      real(dp), allocatable :: gaussian_beam(:,:)

      !Parameters
      logical verbose, remove_monopole,removedipole,polmap,gapfill
      character(LEN=256) :: infile
      character(LEN=256) :: inmaskfile
      character(LEN=256) :: outfile
      character(LEN=256) :: beam_file
      character(LEN=256) :: ordering_out
      real(dp) fwhm_arcmin,gal_cut,target_mean_value,gapfill_fwhm
      integer nlmax,nlmax_beam,gapfill_niter,gapfill_lmax
      integer nside,ord_in,ord_out
      real(dp) rescale_factor, offset,missval,output_missval

      character(len=*), parameter :: CODE = "SMOOTHING_AND_DEGRADE"
      logical treat_missing_pixels,do_smoothing
      real(kind=DP), dimension(0:3,1:2) :: mono_dip
      real(kind=DP) zbounds(2)
      integer(I4B) lowlreg,mm,nmapsout
      real(kind=DP) mean_value
      
#ifdef MPIPIX
      integer mpi_division_method
      call mpi_init(i)
      mpi_division_method = 3
#endif

      
      undef = -1.63750e30
      !----------------------
      !Read in parameter file
      !----------------------
      InputFile = GetParam(1)
      call Ini_Open(InputFile, 1, bad, .false.)
      if (bad) then
         print *,CODE,': Smoothing (+ degrading) of Healpix temperature and polarisation maps.'
         print *,'---------------------'
         print *,'Parameter file format'
         print *,'---------------------'
         print *,'verbose = T or F (default = F).'
         print *,'infile  = Healpix map filename.'
         print *,'inmaskfile  = (Optional) Healpix mask filename for masking data before smoothing.'
         print *,'polmap  = T or F (default = F)'
         print *,'outfile = Output Healpix map filename.'
         print *,'nlmax   = Maximum l for smoothing. (default = 3*nside_in+1)'
         print *,'nside   = Healpix nside for output map.'
         print *,'rescale_factor = Multiplicative factor for map (default = 1.)'
         print *,'offset         = Offset to add to rescaled map (default = 0.)'
         print *,'Then use either:'
         print *,'fwhm_arcmin = Beam FWHM in arcmin.'
         print *,'Or:'
         print *,'beam_file   = Name of B_l file. If both beam_file and fwhm_arcmin are left blank'
         print *,'              then no smoothing is done.'
         print *,'remove_monopole = T or F (default = F)'
         print *,'removedipole    = T or F (default = F, works only if remove_monopole = T)'
         print *,'gal_cut         = Galactic latitude (deg) below which pixels are ignored for removing monopole.'
         print *,'ordering_out    = nest or ring (default = ordering of input map)'
         print *,'target_mean_value = Desired mean value of map, before rescaling (default, parameter not used)'
         print *,'gapfill = Do gap filling via iterative smoothing (default = false)'
         print *,'gapfill_fwhm = FWHM for iterative smoothing'
         print *,'gapfill_niter = Number of iterations iterative smoothing (default = 5)'
         print *,'gapfill_lmax = lmax iterative smoothing (default = 5)'
         print *,'missval = Missing pixel value for setting to zero (default = -1.63750e30)'
         print *,'output_missval = Missing pixel value in output maps (default = -1.63750e30)'
         print *,'---------------------'
         call DoStop(trim(concat(CODE,': Error opening parameter file')))
      endif
      Ini_fail_on_not_found = .false.
      if (InputFile == '') call DoStop('No parameter input file')      

      verbose         = Ini_Read_logical('verbose',.false.)
      infile          = Ini_Read_String('infile')
      inmaskfile      = Ini_Read_String('inmaskfile')
      outfile         = Ini_Read_String('outfile')
      polmap          = Ini_Read_logical('polmap',.false.)
      nlmax           = Ini_Read_int('nlmax',-1)
      nside           = Ini_Read_int('nside',-1)
      rescale_factor  = Ini_Read_Real('rescale_factor',1.)
      offset          = Ini_Read_Real('offset',0.)
      fwhm_arcmin     = Ini_Read_Real('fwhm_arcmin',-1.)
      beam_file       = Ini_Read_String('beam_file')
      remove_monopole = Ini_Read_logical('remove_monopole',.false.)
      removedipole    = Ini_Read_logical('remove_dipole',.false.)
      gal_cut         = Ini_Read_Real('gal_cut',0.)
      ordering_out    = Ini_Read_String('ordering_out')
      target_mean_value = Ini_Read_Real('target_mean_value',-1.)
      gapfill       = Ini_Read_logical('gapfill',.false.)
      gapfill_fwhm  = Ini_Read_Real('gapfill_fwhm',60.)
      gapfill_niter = Ini_Read_int('gapfill_niter',5)
      gapfill_lmax  = Ini_Read_int('gapfill_lmax',1024)
      missval       = Ini_Read_Real('missval',undef)
      output_missval= Ini_Read_Real('output_missval',undef)


      if(verbose) print *,trim(CODE),': Running'

      !-------------
      !Read in a map
      !-------------
      if(verbose) print *,trim(CODE),': Reading ',trim(infile)
      call HealpixMap_read(map,trim(infile))
      ord_in = map%ordering

      !---------------------------------
      !Read in mask for applying to data
      !---------------------------------
      if (inmaskfile .ne. '') then
         if(verbose) print *,trim(CODE),': Reading ',trim(inmaskfile)
         call HealpixMap_read(inmask,trim(inmaskfile))
         if(inmask%nside .ne. map%nside) then
            print *,trim(CODE),': ERROR: inmaskfile and infile have different nside'
         end if
         if(inmask%ordering .ne. map%ordering) then
            print *,trim(CODE),': ERROR: inmaskfile and infile have different Healpix ordering'
         end if

         do pp = 0, map%npix-1
            if(inmask%tqu(pp,1) .lt. 1) then
!               map%tqu(pp,:) = fmissval
               map%tqu(pp,:) = undef
!               print *,fmissval
            endif
         end do
      end if

      !--------------------------------------
      ! Get the number of maps for outputting
      !--------------------------------------
      select case(polmap)
      case(.false.)
         nmapsout = 1
         want_pol = .false.
      case(.true.)
         nmapsout = 3
         want_pol = .true.
      end select

      !-------------------------
      !Get the mask from the map
      !-------------------------
      call HealpixMap_init(mask,map%npix,nmaps=1)
      call get_mask_from_map(map,mask,mask_npix)
!      print *, 'mask_npix ',mask_npix

      !-------------------------------
      ! Parameter defaulting for nlmax
      !-------------------------------
      if (nlmax .eq. -1) then
         nlmax = 3*map%nside + 1
         if(verbose) print *,trim(CODE),': Setting nlmax to '//trim(inttostr(nlmax))
      endif

      !---------------------
      !Read in the beam file
      !---------------------
      if (beam_file .ne. '') then
         print *,trim(CODE),': Reading beam file '//trim(beam_file)
         call ReadBeam(beam, beam_file, nlmax_beam)
         if (fwhm_arcmin .ne. -1.) then
            print *,trim(CODE),': WARNING: Choose fwhm_arcmin or beam_file parameter.'
         endif
         if(nlmax_beam .lt. nlmax) then
            print *,trim(CODE),': using nlmax from signal beam file: '//trim(inttostr(nlmax_beam))
            nlmax = nlmax_beam
         endif
      endif

      !-----------------------
      !or make a Gaussian beam
      !-----------------------
      if (fwhm_arcmin .ne. -1.) then
         allocate(gaussian_beam(0:nlmax,3))
         call gaussbeam(fwhm_arcmin,nlmax,gaussian_beam)
         allocate(beam(0:nlmax))
         beam(:) = gaussian_beam(:,1)
         deallocate(gaussian_beam)
      endif



      !---------------------------------
      ! Gap fill via iterative smoothing
      !---------------------------------
      if(verbose) print *,trim(CODE),': Initialise healpix'
      call HealpixInit(HH,map%nside, nlmax,want_pol, w8dir='')

      if (gapfill) then
         !check to see if mask has holes

         !Set undef pixels to zero
!         where( abs(map%TQU(:,:)/undef -1.) .lt. 1.e-5 )  map%TQU(:,:) = 0.
         do pp = 0, map%npix -1
            if ( abs(map%TQU(pp,1)/undef -1.) .lt. 1.e-5 )  then
               map%TQU(pp,:) = 0.
            endif
         enddo

         call HealpixMap_init(map_smooth,map%npix,nmaps=nmapsout)

         if (HH%MpiID .eq. 0) then
         do ff = 1,gapfill_niter
            print *,trim(CODE),': Gapfill iteration, gapfill_lmax ',ff,gapfill_lmax

            call HealpixMap2alm(HH,map,alm,gapfill_lmax,dopol=want_pol)
            call HealpixAlm_smooth(alm,gapfill_fwhm)
            call HealpixAlm2Map(HH,alm, map_smooth, map%npix) 

            do pp = 0, mask%npix-1
               if( mask%TQU(pp,1) .lt. 1. ) then
                  map%TQU(pp,1) = map_smooth%TQU(pp,1)
                  if(want_pol) then
                     map%TQU(pp,2) = map_smooth%TQU(pp,2)
                     map%TQU(pp,3) = map_smooth%TQU(pp,3)
                  endif
               endif
            end do
         end do
         endif

         !Update the mask.
         call get_mask_from_map(map,mask,mask_npix)
!         print *, 'mask_npix ',mask_npix

      endif

      
      !----------------------------------------------------
      !For handling the case where no smoothing is required
      !----------------------------------------------------
      if (fwhm_arcmin .eq. -1 .and. beam_file .eq. '') then
         do_smoothing = .false.         
         if(verbose) print *,trim(CODE),': No beam information supplied so no smoothing will be applied.'         
      else
         do_smoothing = .true.
      endif      
      
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


      !----------------------------
      !Set undefined pixels to zero
      !----------------------------
!      if (any( abs(map%TQU(:,1)/undef -1.) .lt. 1.e-5)) then 
      if (any( abs(map%TQU(:,1)/missval -1.) .lt. 1.e-5)) then 
         treat_missing_pixels = .true.
         do mm = 1,map%nmaps
            where( abs(map%TQU(:,mm)/missval -1.) .lt. 1.e-5 )  map%TQU(:,mm) = 0.
         end do
      else        
         treat_missing_pixels = .false.
      endif

      !------------------------
      !Remove monopole from map
      !------------------------
      if(remove_monopole) then
         if(verbose) print *,trim(CODE),': Removing monopole'
         call HealpixMap_init(map_tmp,map%npix,nmaps=1)
         do pp = 0, map%npix-1
            map_tmp%tqu(pp,1) = map%tqu(pp,1)
         end do
         zbounds(1) = sin(gal_cut*3.14159/180.)
         zbounds(2) = -zbounds(1) 
         if (removedipole) then
            lowlreg    = 2 ! Remove monopole and dipole
         else
            lowlreg    = 1 ! Remove only the monopole
         endif
         call remove_dipole(map_tmp%nside,map_tmp%tqu(0:,1),map_tmp%ordering,lowlreg,mono_dip(0:,1),&
              zbounds,fmissval = undef, mask = mask%tqu(0:,1))
         if(verbose) print *,trim(CODE),': Monopole value = ',mono_dip(0,1)
         mono_dip(0,1) = -mono_dip(0,1)
         if (removedipole) then
            if(verbose) print *,trim(CODE),': Dipole X value = ',mono_dip(1,1)
            if(verbose) print *,trim(CODE),': Dipole Y value = ',mono_dip(2,1)
            if(verbose) print *,trim(CODE),': Dipole Z value = ',mono_dip(3,1)
            mono_dip(1,1) = -mono_dip(1,1)
            mono_dip(2,1) = -mono_dip(2,1)
            mono_dip(3,1) = -mono_dip(3,1)
         endif
         !------------------------------
         !Subtract monopole (and dipole)
         !------------------------------
         call add_dipole(map%nside,map%tqu(0:,1),map%ordering,lowlreg,&
              mono_dip(0:,1),fmissval=undef)

      end if

      !----------------------------------
      ! Correct the mean value of the map
      !----------------------------------
      if(target_mean_value .ne. -1) then
         if(verbose) print *,trim(CODE),': Adjusting the mean of the map'
         if(removedipole .eqv. .true.) then 
            if(verbose) print *,trim(CODE),': set remove_dipole = F for the case where target_mean_value is set'
            stop
         endif
         npix = 0
         do pp = 0, map%npix-1
            if (map%tqu(pp,1) .ne. undef) then
               mean_value = mean_value + map%tqu(pp,1)
               npix    = npix + 1
            endif
         end do
         mean_value = mean_value/dble(npix)
         if(verbose) print *,trim(CODE),': (Rescaled) mean of the map, offset to add = ',&
              mean_value*rescale_factor,(target_mean_value - mean_value)*rescale_factor
         do pp = 0, map%npix-1
            if (map%tqu(pp,1) .ne. undef) then
               map%tqu(pp,1) = map%tqu(pp,1) + target_mean_value - mean_value
            endif
         end do
      end if


!      if(verbose) print *,trim(CODE),': Initialise healpix'
!      call HealpixInit(HH,map%nside, nlmax,want_pol, w8dir='')
!      call HealpixInit(HH,map%nside, nlmax,want_pol, w8dir='/home/leach/healpix/Healpix_2.13a/data/')

      if (HH%MpiID .eq. 0) then
      !--------------
      !Smooth the map
      !--------------
      if (do_smoothing) then
         if(verbose) print *,trim(CODE),': map2alm, nlmax = '//trim(inttostr(nlmax))
         call HealpixMap2alm(HH,map,alm,nlmax,dopol=want_pol)

         if(verbose) print *,trim(CODE),': Smoothing alm'
         call HealpixAlm_smooth_beam(alm,beam)
         
         if(verbose) print *,trim(CODE),': alm2map '
         call HealpixAlm2Map(HH,alm, map, map%npix) 
      endif
!      call HealpixMap_write(map,'!testsmooth.fits',units='')
      
      !---------------
      !U/Degrade the map
      !---------------
!      if (nside .ne. map%nside) then 
!         npix = nside2npix(nside)
!         if(verbose) print *,trim(CODE),': (U)/Degrading map to nside '//trim(inttostr(nside))
!         call HealpixMap_init(map_out,npix,nmaps=map%nmaps)
!         call HealpixMap_udgrade(map,map_out,nside)
!      else
!         call HealpixMap_init(map_out,map%npix,nmaps=map%nmaps)
!         map_out = map
!      end if
      if (nside .lt. map%nside) then 
         if(verbose) print *,trim(CODE),': Degrading map to nside '//trim(inttostr(nside))
         npix = nside2npix(nside)
         call HealpixMap_init(map_out,npix,nmaps=map%nmaps)
         call HealpixMap_udgrade(map,map_out,nside)
      else if(nside .gt. map%nside) then 
         if(verbose) print *,trim(CODE),': Resythesising map at nside '//trim(inttostr(nside))
         npix = nside2npix(nside)
!SL         call HealpixMap_init(map_out,npix,nmaps=map%nmaps)
         call HealpixMap_init(map_out,npix,nmaps=nmapsout)
         if(verbose) print *,trim(CODE),': map2alm, nlmax = '//trim(inttostr(nlmax))
         call HealpixInit(HH2,nside, nlmax,want_pol, w8dir='')
         call HealpixAlm2Map(HH2,alm, map_out, npix) 
         call HealpixFree(HH2)
      else
         call HealpixMap_init(map_out,map%npix,nmaps=map%nmaps)
         map_out = map
      end if
!      call HealpixMap_write(map,'!testdegrade.fits',units='')

      !------------------------------
      !Rescale the map and add offset
      !------------------------------
      do pp = 0, map_out%npix-1
         map_out%tqu(pp,1) = rescale_factor*map_out%tqu(pp,1) + offset
         if (want_pol) then
            map_out%tqu(pp,2) = rescale_factor*map_out%tqu(pp,2)
            map_out%tqu(pp,3) = rescale_factor*map_out%tqu(pp,3)
         end if
      end do

      !-------------------------------------------
      !Mask away pixels that are close to the mask
      !-------------------------------------------
      if (treat_missing_pixels) then 
         if(verbose) print *,trim(CODE),': Dealing with missing pixels'
         if(do_smoothing) then
            if(verbose) print *,trim(CODE),': Smoothing mask'
            call HealpixMap2alm(HH,mask,alm,nlmax)
            call HealpixAlm_smooth_beam(alm,beam)
            call HealpixAlm2Map(HH,alm, mask, mask%npix) 
         end if
         if(nside .lt. map%nside) then
            if(verbose) print *,trim(CODE),': Degrading mask'
            npix = nside2npix(nside)
            call HealpixMap_init(map_tmp,npix,nmaps=1)
            call HealpixMap_udgrade(mask,map_tmp,nside)
            mask = map_tmp          
         end if
!         call HealpixMap_write(mask,'!testmask1.fits',units='')

         if(verbose) print *,trim(CODE),': Trimming mask and cropping data'
         do pp = 0, mask%npix-1
            if (1.003*mask%TQU(pp,1) .gt. 0.99999) mask%TQU(pp,1) = 1.
            if (mask%TQU(pp,1) .lt. 1. )           mask%TQU(pp,1) = 0.
         end do
!         call HealpixMap_write(mask,'!testmask2.fits',units='')

         do pp = 0, mask%npix-1
            if( mask%TQU(pp,1) .lt. 1. )     map_out%TQU(pp,1) = output_missval
         end do

         if (want_pol) then
            do pp = 0, mask%npix-1
               if( mask%TQU(pp,1) .lt. 1. ) then
                  map_out%TQU(pp,2) = output_missval
                  map_out%TQU(pp,3) = output_missval
               endif
            end do
         endif
         
      endif      

      !-----------
      !Reorder map
      !-----------
      select case(ord_out)
      case(ord_ring)
         if(verbose) print *,trim(CODE)//': Forcing output map to ring ordering.'
         call HealpixMap_ForceRing(map_out)
      case(ord_nest)
         if(verbose) print *,trim(CODE)//': Forcing output map to nest ordering.'
         call HealpixMap_ForceNest(map_out)
      end select
      
      
      if(verbose) print *,trim(CODE)//': Writing map ',trim(outfile)
      call HealpixMap_write(map_out,trim(outfile),units='')
   endif

#ifdef MPIPIX
   call HealpixFree(HH)
   call mpi_finalize(i)
#endif

   
 end program smoothing_and_degrade
