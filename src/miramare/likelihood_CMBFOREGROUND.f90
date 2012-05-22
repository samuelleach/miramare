module Likelihoodcmbforeground
  use InstrumentData
  use MatrixUtils
  use ParamDef
  use HealpixObj
  use foregroundscalings
  use foreground_simulations
  use noise_simulations
  use runparameters
  use CMBData
  use AMLutils
  use omp_lib
  implicit none

  !-----------
  !GLOBAL DATA
  !-----------
  Type(HealpixMap),    allocatable, target :: MapDataArray(:)        ! Data
  Type(HealpixMap),    allocatable, target :: MapWeightArray(:)      ! Weight (N^-1)
  Type(HealpixMap),    allocatable, target :: MapDataWeightArray(:)  ! Weight x data
  Type(RegionalModel), allocatable, target :: Region(:)              ! Foreground model
  integer*8,           allocatable         :: pixlist_by_region(:)   ! Pixel list

  type(ParamSet) LastParams
  real LastLnLike

  type(HealpixMask) ::  Hmask! Defines region for simulations and single region analysis

 
contains

 subroutine Likelihood_init
   !----------------------------------------------------
   !Initialise global data.  
   !Reads in or simulates Q/U Healpix maps and templates
   !Read in Q/U Healpix weights
   !----------------------------------------------------

   implicit none
   
   integer stokes_index,ff,ss,tt,cc
   character(LEN=512) filename
   type(HealpixMap) hitmap, Htemplate
   real factor, dummy
   logical all_simulations_cached, writetocache
   integer comp_list(1),rr
   integer*8 counter,pp,mm
   integer nweightmaps
   real(dp), allocatable :: calibfactor(:)
   real(dp), allocatable :: unity(:)
   real(dp) Amat(rp%nchannel,1)
   character(32) :: component_list(1)
   real(dp)      :: fg_param(1)
   real(dp)      :: reference_frequency_list(1)

   if (rp%Feedback > 0) print *,'-----------------------'
   if (rp%Feedback > 0) print *,'Initialising likelihood'
   if (rp%Feedback > 0) print *,'-----------------------'

   !--------------------
   !Allocate global data
   !--------------------
   allocate(MapDataArray(rp%nchannel))
   allocate(MapWeightArray(rp%nchannel))
   allocate(unity(rp%nchannel))
   select case(rp%weight_format)
   case(5)
      nweightmaps    = 6
   case default
      if (rp%want_pol) then
         nweightmaps = 3
      else
         nweightmaps = 1
      end if
   end select

   do ff = 1,rp%nchannel 
      unity(ff) = 1.
      call HealpixMap_Nullify(MapWeightArray(ff))
!      call HealpixMap_Init(MapWeightArray(ff), nside2npix(rp%nside_in), pol = rp%want_pol)
      call HealpixMap_Init(MapWeightArray(ff), nside2npix(rp%nside_in), nweightmaps)
      !Getting ready for reduced memory mode.
!!      call HealpixMap_Init(MapWeightArray(ff), rp%npixtot, pol = rp%want_pol)
      call HealpixMap_Nullify(MapDataArray(ff))
      call HealpixMap_Init(MapDataArray(ff), nside2npix(rp%nside_in), pol = rp%want_pol)
      !Getting ready for reduced memory mode.
!!      call HealpixMap_Init(MapDataArray(ff), rp%npixtot, pol = rp%want_pol)
   end do

   !Getting ready for reduced memory mode.
   allocate(pixlist_by_region(rp%npixtot))
   counter = 0
   do rr = 1, rp%nregion
      do pp = 1, region(rr)%npix
         counter                    = counter + 1
         pixlist_by_region(counter) = region(rr)%pix(pp)
      end do
   enddo

   !-------------------------------------------
   !Check for cached simulated data and weights
   !-------------------------------------------
   all_simulations_cached = .false.
   if(rp%sim_signal) then
      all_simulations_cached = .true.
      if (rp%feedback > 0) print *,' Checking for cached simulations in '//trim(rp%cachedir)
      do ff = 1, rp%nchannel 
         filename = trim(rp%cachedir)//'/'//trim(rp%rootname)//'_signal_'//&
              trim(inttostr(ff))//'.fits'
         if (.not. FileExists(filename)) then
            if (rp%feedback > 0) print *,' Cannot find cached simulation '//trim(filename)
            all_simulations_cached = .false.
         endif
         
         filename = trim(rp%cachedir)//'/'//trim(rp%rootname)//'_weight_'//&
              trim(inttostr(ff))//'.fits'
         if (.not. FileExists(filename)) then
            if (rp%feedback > 0) print *,' Cannot find cached weight '//trim(filename)
            all_simulations_cached = .false.
         endif
      end do
   end if

   if(all_simulations_cached .and. rp%do_compsep_on_simulated_component .eq. -1) then
      !--------------------------
      !Read in cached simulations
      !--------------------------
      do ff = 1,rp%nchannel 
         filename = trim(rp%cachedir)//'/'//&
              trim(rp%rootname)//'_signal_'//trim(inttostr(ff))//'.fits'
         if (rp%feedback > 0) print *, ' Reading simulated data ',trim(filename)
         call HealpixMap_read(MapDataArray(ff),filename)
         !Getting ready for reduced memory mode ??
         !!         call HealpixMap_readpixlist(MapDataArray(ff),filename,pixlist_by_region,rp%npixtot)
         filename = trim(rp%cachedir)//'/'//&
              trim(rp%rootname)//'_weight_'//trim(inttostr(ff))//'.fits'
         if (rp%feedback > 0) print *, ' Reading simulated weight ',trim(filename)
         call HealpixMap_read(MapWeightArray(ff),filename)
         !Getting ready for reduced memory mode ??
      end do
   else
      !-----------------------------------
      !Read in data or perform simulations
      !-----------------------------------

      !-------------------------------
      !First get recalibration factors
      !-------------------------------
      allocate(calibfactor(rp%nchannel))
      if(rp%sim_caliberror) then
         if (rp%Feedback>0) print *,' Simulating calibration errors.'
         call InitRandom(rp%caliberror_rand_seed)
         dummy = gaussian1()         

         do ff = 1, rp%nchannel 
            calibfactor(ff) = rp%sim_calib(ff) + rp%sim_calib_rms(ff)*Gaussian1()
            if (rp%Feedback>1) print *,' Channel '//trim(inttostr(ff))//&
                 ', Simulated calibration factor = ',real(calibfactor(ff))
         end do
         
         !-------------------------------------------
         !Write simulated calibration factors to disk
         !-------------------------------------------
         filename = trim(rp%cachedir)//'/'//trim(rp%rootname)//'_calibfactor.txt'
         if (rp%feedback > 0) print *, ' Writing simulated calibration factors to ',trim(filename)
         call CreateTxtFile(filename,tmp_file_unit)
         do ff = 1, rp%nchannel 
            write (tmp_file_unit,*) trim('calib(')//trim(inttostr(ff))//') =',calibfactor(ff)
         end do
         close(tmp_file_unit)
      else
         do ff = 1, rp%nchannel 
            calibfactor(ff) = rp%calib(ff) 
            if (rp%Feedback>1) print *,' Channel '//trim(inttostr(ff))//', Calibration factor = ',calibfactor(ff)
         end do
      endif
      
      !------------
      !Read in data
      !------------
      if (.not. rp%sim_signal) then
         do ff = 1, rp%nchannel 
            filename = trim(rp%datafile(ff))      

            if (filename .ne. '') then
               if (rp%feedback > 0) print *, ' Reading data ',trim(filename)
               call HealpixMap_read(MapDataArray(ff),filename)
               !Getting ready for reduced memory mode.
               !!            call HealpixMap_readpixlist(MapDataArray(ff),filename,pixlist_by_region,rp%npixtot)
               do pp = 0, MapDataArray(ff)%npix-1
                  do ss = 1, rp%nstokes
                     MapDataArray(ff)%TQU(pp,rp%stokes_list(ss)) = &
                          MapDataArray(ff)%TQU(pp,rp%stokes_list(ss)) * calibfactor(ff)
                  end do
               end do
            else
               if (rp%feedback > 0) print *, ' WARNING: datafile is not set for channel ',ff
               if (rp%feedback > 0) print *, ' WARNING: Setting data to zero for channel ',ff
               do pp = 0, MapDataArray(ff)%npix-1
                  do ss = 1, rp%nstokes
                     MapDataArray(ff)%TQU(pp,rp%stokes_list(ss)) = 0.
                  end do
               end do
            endif
         end do
      endif
      
      !---------------
      !Read in weights
      !---------------
      !----------------------------------------------------------------------------------------
      !    weight_format = 1 => Read IQU weight maps (units 1/data^2)
      !    weight_format = 2 => Read hit maps + I and QU pixel RMS.
      !    weight_format = 3 => Read "WMAP style" hits in column 4 of data + I and QU pixel RMS
      !    weight_format = 4 => Homogeneous noise => I and QU pixel RMS.
      !    weight_format = 5 => "Madam" style 3x3 IQU inverse noise covariance.
      !----------------------------------------------------------------------------------------
      select case(rp%weight_format)
      case(1)
         do ff = 1,rp%nchannel 
            filename = trim(rp%weightfile(ff))
            if (rp%feedback > 0) print *, ' Reading weight ',trim(filename)
            call HealpixMap_read(MapWeightArray(ff),filename)
            !Getting ready for reduced memory mode.
!!            call HealpixMap_readpixlist(MapWeightArray(ff),filename,pixlist_by_region,rp%npixtot)
         end do
      case(2)
         call HealpixMap_Nullify(hitmap)
!         call HealpixMap_Init(hitmap, nside2npix(rp%nside_in))
         !-------------------------------------------------
         !Read hitmaps from column 1 of a separate FITS map
         !-------------------------------------------------
         do ff = 1, rp%nchannel 
            filename = trim(rp%weightfile(ff))
            factor   = 1.
            if (rp%feedback > 0) print *, ' Reading hitmap ',trim(filename)
            call HealpixMap_read(hitmap,filename)
            !Getting ready for reduced memory mode.
!!            call HealpixMap_readpixlist(hitmap,filename,pixlist_by_region,rp%npixtot)
            !Convert hit maps to weights
            do ss = 1, rp%nstokes
               stokes_index = rp%stokes_list(ss)
               if (stokes_index .eq. 1) then
                  factor = 1./rp%stokesi_noise_sig0(ff)**2  
                  MapWeightArray(ff)%TQU(Hmask%good,stokes_index)=&
                       hitmap%tqu(Hmask%good,1)*factor
               else
                  factor = 1./rp%stokesqu_noise_sig0(ff)**2  
                  MapWeightArray(ff)%TQU(Hmask%good,stokes_index)=&
                       hitmap%tqu(Hmask%good,1)*factor
               end if
            end do
         end do
         call HealpixMap_Free(Hitmap)
      case(3) !WMAP style (hitmap in 4 column of IQU map)
         !Read hitmaps from column 4 of FITS data map      
         do ff = 1, rp%nchannel
            filename = trim(rp%datafile(ff))
            if (rp%feedback > 0) print *, ' Reading ',trim(filename)
            call HealpixMap_read(hitmap,filename)
            !Getting ready for reduced memory mode.
!!           HealpixMap_readpixlist(hitmap,filename,pixlist_by_region,rp%npixtot)
            !Convert hit maps to weights
            do ss = 1,rp%nstokes
               stokes_index = rp%stokes_list(ss)
               if (stokes_index .eq. 1) then
                  factor = 1./rp%stokesi_noise_sig0(ff)**2  
                  MapWeightArray(ff)%TQU(Hmask%good,stokes_index)=&
                       hitmap%tqu(Hmask%good,4)*factor
               else
                  factor = 1./rp%stokesqu_noise_sig0(ff)**2  
                  MapWeightArray(ff)%TQU(Hmask%good,stokes_index)=&
                       hitmap%tqu(Hmask%good,4)*factor
               end if
            end do
         end do
         call HealpixMap_Free(Hitmap)
      case(4) !Homogeneous noise case
         do ff = 1,rp%nchannel
            do ss = 1, rp%nstokes
               stokes_index = rp%stokes_list(ss)
               if (stokes_index .eq. 1) then
                  factor = 1./rp%stokesi_noise_sig0(ff)**2
               else
                  factor = 1./rp%stokesqu_noise_sig0(ff)**2
               endif
               MapWeightArray(ff)%TQU(Hmask%good,stokes_index) = factor
            enddo
         enddo
      case(5)
         do ff = 1,rp%nchannel 
            filename = trim(rp%weightfile(ff))
            if (rp%feedback > 0) print *, ' Reading 3x3 inverse noise covariance  ',trim(filename)
            call HealpixMap_read(MapWeightArray(ff),filename)
            !Getting ready for reduced memory mode.
!!            call HealpixMap_readpixlist(MapWeightArray(ff),filename,pixlist_by_region,rp%npixtot)
         end do
      end select

      !-----------------------
      !Recalibrate the weights
      !-----------------------
      if (.not. rp%sim_caliberror) then
         if (rp%Feedback>0) print *,' Recalibrating weights using calib(i) parameters'
!         do ff = 1, rp%nchannel 
!            MapWeightArray(ff)%TQU(Hmask%good,:) = MapWeightArray(ff)%TQU(Hmask%good,:)/&
!                 calibfactor(ff)**2
!         end do
         do ff = 1, rp%nchannel 
            do pp = 1, Hmask%npix 
               do ss = 1, rp%nstokes
                  MapWeightArray(ff)%TQU(Hmask%good(pp),rp%stokes_list(ss)) = &
                       MapWeightArray(ff)%TQU(Hmask%good(pp),rp%stokes_list(ss))/&
                       calibfactor(ff)**2
               end do
            end do
         end do
      else
         if (rp%Feedback>0) print *,' sim_caliberror = T => No recalibration of weights using calib(i) parameters'
      endif  !??? revisit this
         

      !------------------------------------------------------------------------------
      !Add templates e.g. Subtract free-free (WMAP analysis) or add dust (simulation)
      !------------------------------------------------------------------------------
      if(rp%ntemplate .gt. 0) then
         do tt = 1, rp%ntemplate
            if (rp%Feedback>0) print *,' Reading template from ',trim(rp%templatefile(tt))         
            call HealpixMap_read(Htemplate,trim(rp%templatefile(tt)))
            !Getting ready for reduced memory mode.
            !!           HealpixMap_readpixlist(Htemplate,filename,pixlist_by_region,rp%npixtot)
            do ff = 1 ,rp%nchannel 
               do pp = 0, MapDataArray(ff)%npix-1
                  do ss = 1, rp%nstokes
                     MapDataArray(ff)%TQU(pp,rp%stokes_list(ss)) = &
                          MapDataArray(ff)%TQU(pp,rp%stokes_list(ss)) + &
                          Htemplate%TQU(pp,rp%stokes_list(ss))*rp%template_coefficient(tt,ff)
                  end do
               end do
            end do
         end do
         call HealpixMap_Free(Htemplate)
      endif


      !----------------
      !Simulate the CMB
      !----------------
      if(rp%sim_signal) then
         if(rp%do_compsep_on_simulated_component .lt. 1 .and. rp%do_compsep_on_simulated_component .gt. -2) then
            writetocache = .true.
            if (rp%Feedback>0) print *,' Getting CMB simulation'
            call get_cmb_simulation(rp%Inst%channel_nu0_ghz,MapDataArray,rp%nchannel,rp,&
                 weight = calibfactor,writetocache=writetocache)
            !Getting ready for reduced memory mode.
            !         call get_cmb_simulation(rp%Inst%channel_nu0_ghz,MapDataArray,rp%nchannel,rp,&
            !              weight=calibfactor,writetocache=.true.,pixlist=pixlist_by_region)
         endif
      endif
      !----------------------------------------------------------
      !Add simulated noise using the level set by the weight file
      !----------------------------------------------------------
      if(rp%sim_noise) then
         writetocache = .true.
         if(rp%do_compsep_on_simulated_component .gt. -1) writetocache = .false.
         call get_noise_simulation(MapWeightArray,MapDataArray,rp%nchannel,rp,&
              writetocache=writetocache)
            !Getting ready for reduced memory mode.
      end if

      !------------------------
      !Add simulated components
      !------------------------
      if(rp%sim_component) then
         do cc = 1, rp%sim_ncomponent
            comp_list(1) = cc
            if((rp%do_compsep_on_simulated_component .eq. -1) .or.&
                 (rp%do_compsep_on_simulated_component .eq. cc)) then
               call get_fg_simulation(comp_list,1,rp%Inst%channel_nu0_ghz,& 
                    MapDataArray,rp%nchannel,rp,weight=calibfactor)
               !Getting ready for reduced memory mode.
            endif
         end do
      endif

      !------------------
      !Add offset to data
      !------------------
      do ff = 1 , rp%nchannel
         do pp = 0, MapDataArray(ff)%npix-1
            do mm = 1, MapDataArray(ff)%nmaps               
               MapDataArray(ff)%TQU(pp,mm) = MapDataArray(ff)%TQU(pp,mm) + rp%offset(ff)
            end do
         end do
      end do
      
      !--------------------------------------------
      !Cache the simulated data and weights to disk
      !--------------------------------------------
      if(rp%sim_noise .or. rp%sim_signal .or. rp%sim_component) then
         if(rp%do_compsep_on_simulated_component .eq. -1) then
            do ff = 1, rp%nchannel 
!               filename = trim(rp%cachedir)//'/'//trim(rp%rootname)//'_signal_'//trim(inttostr(ff))//'.fits'
               filename = "!"//trim(rp%cachedir)//'/'//trim(rp%rootname)//'_signal_'//trim(inttostr(ff))//'.fits'
               if (rp%feedback > 0) print *, ' Writing ',trim(filename)
               call HealpixMap_write(MapDataArray(ff),filename)
!               filename = trim(rp%cachedir)//'/'//trim(rp%rootname)//'_weight_'//trim(inttostr(ff))//'.fits'
               filename = "!"//trim(rp%cachedir)//'/'//trim(rp%rootname)//'_weight_'//trim(inttostr(ff))//'.fits'
               if (rp%feedback > 0) print *, ' Writing ',trim(filename)
               call HealpixMap_write(MapWeightArray(ff),filename)
            end do
         endif
      endif

   endif

   !------------------------------------
   !Remove raw offsets (average of data)
   !------------------------------------
   if (rp%remove_raw_offsets) then
      do ff = 1, rp%nchannel 
         do ss = 1, rp%nstokes
            rp%Inst%channel_raw_offset(ff,ss) = &
                 sum(MapDataArray(ff)%TQU(Hmask%good,rp%stokes_list(ss)))
            rp%Inst%channel_raw_offset(ff,ss) = &
                 rp%Inst%channel_raw_offset(ff,ss)/real(Hmask%npix)
            if (rp%Feedback>0) print *,' Raw offset of channel,stokes '//trim(inttostr(ff))//' '//&
                 trim(inttostr(ss))//' = ',rp%Inst%channel_raw_offset(ff,ss)
            MapDataArray(ff)%TQU(:,rp%stokes_list(ss)) = &
                 MapDataArray(ff)%TQU(:,rp%stokes_list(ss))-rp%Inst%channel_raw_offset(ff,ss)
         end do
      end do
   endif

   !------------------------
   ! Perform unit conversion
   !------------------------
   if (.not. rp%sim_signal .and. rp%convert_units_cmb_to_rj) then
      if (rp%Feedback>0) print *,' Performing unit conversion from CMB to RJ units.'
      !--------------------------
      ! Get mixing matrix for CMB
      !--------------------------
      component_list(1) = 'cmb'
      fg_param(1)       = 1.
      Amat = get_mixingmatrix(component_list,reference_frequency_list,1,fg_param,&
           rp%Inst%bandpass_ghz,rp%inst%bandpass_weight,&
           rp%nchannel,rp%Inst%nabscissa)                  
      do ff = 1,rp%nchannel
         if (rp%Feedback>0) print *,'  Channel '//trim(inttostr(ff))//' x ',Amat(ff,1)
         do pp = 0, MapDataArray(ff)%npix -1
            do ss = 1, rp%nstokes
               MapDataArray(ff)%TQU(pp,rp%stokes_list(ss))   = MapDataArray(ff)%TQU(pp,rp%stokes_list(ss))*Amat(ff,1)
               MapWeightArray(ff)%TQU(pp,rp%stokes_list(ss)) = MapWeightArray(ff)%TQU(pp,rp%stokes_list(ss))/Amat(ff,1)**2
            end do
         end do
!         MapDataArray(ff)%TQU   = MapDataArray(ff)%TQU*Amat(ff,1)
!         MapWeightArray(ff)%TQU = MapWeightArray(ff)%TQU/Amat(ff,1)**2
      end do
   endif

   
   !----------------------------------------------
   ! Multiply the weight maps by a constant factor
   !----------------------------------------------
   do ff = 1,rp%nchannel
      do ss = 1, rp%nstokes
         do pp = 0, MapWeightArray(ff)%npix-1
            MapWeightArray(ff)%TQU(pp,rp%stokes_list(ss)) = &
                 MapWeightArray(ff)%TQU(pp,rp%stokes_list(ss))*rp%weight_factor(ff)
         end do
      end do
   end do


   if (rp%Feedback>0) print *,'--------------------------------'
   if (rp%Feedback>0) print *,'Finished initialising likelihood'
   if (rp%Feedback>0) print *,'--------------------------------'

 end subroutine Likelihood_init

 subroutine Multiply_Ninv_data
   !-------------------------------------
   !One time multiplcation of Ninv x data
   !-------------------------------------
   integer ff,cc,ss
   integer*8 pp

   allocate(MapDataWeightArray(rp%nchannel))
   do ff = 1, rp%nchannel
      call HealpixMap_Nullify(MapDataWeightArray(ff))
      call HealpixMap_Init(MapDataWeightArray(ff), nside2npix(rp%nside_in), pol = rp%want_pol)
      !Getting ready for reduced memory mode.
      !call HealpixMap_Init(MapDataWeightArray(ff), rp%npixtot, pol = rp%want_pol)
      select case(MapWeightArray(ff)%nmaps)
      case(6)
!         continue
!         !INSERT MATRIX VECTOR MULTIPLICATION HERE.
         select case(rp%stokes)
            case('IQU')
               do pp = 0, nside2npix(rp%nside_in) - 1
                  MapDataWeightArray(ff)%TQU(pp,1) = &
                       MapDataArray(ff)%TQU(pp,1)*MapWeightArray(ff)%TQU(pp,1)+&
                       MapDataArray(ff)%TQU(pp,2)*MapWeightArray(ff)%TQU(pp,2)+&
                       MapDataArray(ff)%TQU(pp,3)*MapWeightArray(ff)%TQU(pp,3)
                  MapDataWeightArray(ff)%TQU(pp,2) = &
                       MapDataArray(ff)%TQU(pp,1)*MapWeightArray(ff)%TQU(pp,2)+&
                       MapDataArray(ff)%TQU(pp,2)*MapWeightArray(ff)%TQU(pp,4)+&
                       MapDataArray(ff)%TQU(pp,3)*MapWeightArray(ff)%TQU(pp,5)
                  MapDataWeightArray(ff)%TQU(pp,3) = &
                       MapDataArray(ff)%TQU(pp,1)*MapWeightArray(ff)%TQU(pp,3)+&
                       MapDataArray(ff)%TQU(pp,2)*MapWeightArray(ff)%TQU(pp,5)+&
                       MapDataArray(ff)%TQU(pp,3)*MapWeightArray(ff)%TQU(pp,6)
               end do
            case('QU')
               do pp = 0, nside2npix(rp%nside_in) - 1
                  MapDataWeightArray(ff)%TQU(pp,2) = &
                       MapDataArray(ff)%TQU(pp,2)*MapWeightArray(ff)%TQU(pp,4)+&
                       MapDataArray(ff)%TQU(pp,3)*MapWeightArray(ff)%TQU(pp,5)
                  MapDataWeightArray(ff)%TQU(pp,3) = &
                       MapDataArray(ff)%TQU(pp,2)*MapWeightArray(ff)%TQU(pp,5)+&
                       MapDataArray(ff)%TQU(pp,3)*MapWeightArray(ff)%TQU(pp,6)
               end do
            case default
               print *,'Fix matrix multiplication'
            end select

      case default
         do ss = 1, rp%nstokes
            do pp = 0, nside2npix(rp%nside_in) - 1
               MapDataWeightArray(ff)%TQU(pp,rp%stokes_list(ss)) = &
                    MapDataArray(ff)%TQU(pp,rp%stokes_list(ss))*&
                    MapWeightArray(ff)%TQU(pp,rp%stokes_list(ss))
            end do
         end do
      end select
   end do
   do cc = 1, rp%nchannel
      call HealpixMap_free(MapDataArray(cc))             
   enddo
   deallocate(MapDataArray)
 end subroutine Multiply_Ninv_data

 function cmbforeground_LnLike4(Params) result (LnLike)
   !--------------------------------------
   !Multi-region spectral index likelihood
   !--------------------------------------

   implicit none
   type(ParamSet) Params
   real LnLike
   integer ii,ff,pp,ss,rr,firstparam_index,lastparam_index,cc
   integer*8 pix_index1,pix_index2

   real(dp) tmap(rp%nchannel)
   real(dp) invdiagNoise(rp%nchannel,rp%nchannel)
   real(dp) tmpmat(rp%nchannel,rp%ncomponent)
   real(dp) tmpmat_reduced(rp%nchannel,rp%npolcomponent)
   real(dp), allocatable, dimension(:,:) :: A, NinvA   !(rp%nchannel,rp%npolcomponent) 
   real(dp), allocatable, dimension(:,:) :: At
   real(dp), allocatable, dimension(:)   :: chi2
   real(dp), allocatable, dimension(:)   :: smap
   real(dp), allocatable, dimension(:,:) :: nmat
   real(dp), allocatable, dimension(:)   :: cmap
   real(dp) sum_chi2,factor

   sum_chi2          = 0.d0
   LnLike            = 0.d0
   invdiagNoise(:,:) = 0.d0
   do rr = 1, rp%nregion !Loop over regions
      lastparam_index  = region(rr)%lastparam_index
      firstparam_index = region(rr)%firstparam_index

      if (all(abs(Params%P(firstparam_index:lastparam_index) - &
           LastParams%P(firstparam_index:lastparam_index)) < 1d-10)) then
         !---------------------------------------------------------------------------
         !.. the parameters, and hence the likelihood of this region have not changed
         !---------------------------------------------------------------------------
         continue
      else
         !--------------------
         !Set up mixing matrix
         !--------------------
         tmpmat = get_mixingmatrix(region(rr)%component_list,&
              rp%component_reference_frequency_list,region(rr)%ncomponent, &
              Params%P(firstparam_index:lastparam_index),&
              rp%Inst%bandpass_ghz,rp%inst%bandpass_weight,rp%nchannel,&
              rp%Inst%nabscissa)         
         
         do ss = 1,rp%nstokes !Loop over Stokes parameters
            !--------------------------------------------------------------------------
            ! Set to zero the Q and U mixing matrix elements for unpolarized components
            !--------------------------------------------------------------------------
            if (rp%stokes_list(ss) .gt. 1 .and. rp%npolcomponent .lt. rp%ncomponent) then
               do cc=1,rp%npolcomponent
                  tmpmat_reduced(:,cc) = tmpmat(:,rp%polcomponent_list(cc))
               end do
               allocate(A(rp%nchannel,rp%npolcomponent))
               allocate(NinvA(rp%nchannel,rp%npolcomponent))
               allocate(At(rp%npolcomponent,rp%nchannel))
               allocate(cmap(rp%npolcomponent))
               allocate(smap(rp%npolcomponent))
               allocate(nmat(rp%npolcomponent,rp%npolcomponent))
               A(:,:) = tmpmat_reduced(:,:)
               At     = transpose(A)
            else
               allocate(A(rp%nchannel,rp%ncomponent))
               allocate(NinvA(rp%nchannel,rp%ncomponent))
               allocate(At(rp%ncomponent,rp%nchannel))
               allocate(cmap(rp%ncomponent))
               allocate(smap(rp%ncomponent))
               allocate(nmat(rp%ncomponent,rp%ncomponent))
               A(:,:) = tmpmat(:,:)
               At     = transpose(A)
            endif

            pix_index1 = 1
            pix_index2 = region(rr)%npix
            allocate(chi2(region(rr)%npix))
            !Getting ready for reduced memory mode.
!            pix_index1 = region(rr)%pix_index1   
!            pix_index2 = region(rr)%pix_index1 + region(rr)%npix - 1

!$OMP parallel do private(pp,tmap,NinvA,smap,nmat,cmap,ff) shared(rr,ss,chi2,MapWeightArray,MapDataWeightArray,A,At,pix_index1,pix_index2,rp)
            do pp = pix_index1, pix_index2, rp%modulo_npix !Loop over pixels
               !-------------------
               !Set up noise matrix
               !-------------------
               do ff = 1, rp%nchannel     
                  tmap(ff)            = MapDataWeightArray(ff)%TQU(region(rr)%pix(pp),rp%stokes_list(ss))
                  !Getting ready for reduced memory mode.
                  !invdiagNoise(ff,ff) = MapWeightArray(ff)%TQU(pp,rp%stokes_list(ss))
                  !tmap(ff)            = MapDataWeightArray(ff)%TQU(pp,rp%stokes_list(ss))
                  NinvA(ff,:)          = A(ff,:)*MapWeightArray(ff)%TQU(region(rr)%pix(pp),rp%stokes_list(ss))
               end do
               
               !-----------------
               !Matrix operations
               !-----------------
               smap     = matmul(At,tmap)
               nmat     = matmul(At,NinvA) 
               call Matrix_Inverse(nmat)
               cmap     = matmul(nmat,smap)
               chi2(pp) = dot_product(smap,cmap)

!               call matrix_mulvec(At,tmap,smap)
!               call matrix_mult(At,NinvA,nmat)
!               call Matrix_Inverse(nmat)
!               call matrix_mulvec(nmat,smap,cmap)
!               chi2(pp) = Matrix_vecdot(smap,cmap)
               
            end do !Loop over pixels
            deallocate(A,At,smap,cmap,nmat,NinvA)

            region(rr)%chi2 = 0.0d0
            do pp = pix_index1, pix_index2, rp%modulo_npix
               region(rr)%chi2 = region(rr)%chi2 + chi2(pp)
            enddo
            deallocate(chi2)

         end do !Loop over Stokes

         if(use_prior) then
            factor = real(rp%nstokes)
            do ii = firstparam_index,lastparam_index
               if (abs(scales%pvar(ii)) > 1d-10) then        
                  region(rr)%chi2 = region(rr)%chi2 - &
                       factor*(params%p(ii)-scales%pmean(ii))**2/scales%pvar(ii) - &
                       factor*dlog(scales%pvar(ii))
               endif
            end do
         end if

      end if
      sum_chi2 = sum_chi2 + region(rr)%chi2
   end do !Loop over region

   LnLike = -0.5d0*sum_chi2 + rp%lnlike_offset
   if (rp%Feedback>2) then
      print *, 'Parameters : ',Params%P(1:rp%nregion*rp%nparam_component)
      print *, 'Ln(like) = '  ,lnlike
!      print *, 'Ln(like), offset = ',  lnlike,-0.5d0*sum_chi2,rp%lnlike_offset
!      if (lnlike 
   end if

 end function cmbforeground_LnLike4

 function cmbforeground_LnLike5(Params) result (LnLike)
   !----------------------------------------------------------------------------
   !Multi-region spectral index likelihood, with assumption of homogeneous noise
   !----------------------------------------------------------------------------

   implicit none
   type(ParamSet) Params
   real LnLike
   integer ii,ff,pp,ss,rr,firstparam_index,lastparam_index,cc
   integer*8 pix_index1,pix_index2

   real(dp) tmap(rp%nchannel)
   real(dp) invdiagNoise(rp%nchannel,rp%nchannel)
   real(dp) tmpmat(rp%nchannel,rp%ncomponent)
   real(dp) tmpmat_reduced(rp%nchannel,rp%npolcomponent)
   real(dp), allocatable, dimension(:,:) :: A, NinvA   !(rp%nchannel,rp%npolcomponent) 
   real(dp), allocatable, dimension(:,:) :: At
   real(dp), allocatable, dimension(:)   :: smap,chi2
   real(dp), allocatable, dimension(:,:) :: nmat
   real(dp), allocatable, dimension(:)   :: cmap
   real(dp) sum_chi2,factor

   sum_chi2          = 0.d0
   LnLike            = 0.d0
   invdiagNoise(:,:) = 0.d0
   do rr = 1, rp%nregion !Loop over regions
      lastparam_index  = region(rr)%lastparam_index
      firstparam_index = region(rr)%firstparam_index

      if (all(abs(Params%P(firstparam_index:lastparam_index) - &
           LastParams%P(firstparam_index:lastparam_index)) < 1d-10)) then
         !---------------------------------------------------------------------------
         !.. the parameters, and hence the likelihood of this region have not changed
         !---------------------------------------------------------------------------
         continue
      else
         region(rr)%chi2 = 0.0d0
         !--------------------
         !Set up mixing matrix
         !--------------------
         tmpmat = get_mixingmatrix(region(rr)%component_list,&
              rp%component_reference_frequency_list,region(rr)%ncomponent, &
              Params%P(firstparam_index:lastparam_index),&
              rp%Inst%bandpass_ghz,rp%inst%bandpass_weight,rp%nchannel,&
              rp%Inst%nabscissa)         
         
         do ss = 1,rp%nstokes !Loop over Stokes parameters
            !--------------------------------------------------------------------------
            ! Set to zero the Q and U mixing matrix elements for unpolarized components
            !--------------------------------------------------------------------------
            if (rp%stokes_list(ss) .gt. 1 .and. rp%npolcomponent .lt. rp%ncomponent) then
               do cc=1,rp%npolcomponent
                  tmpmat_reduced(:,cc) = tmpmat(:,rp%polcomponent_list(cc))
               end do
               allocate(A(rp%nchannel,rp%npolcomponent))
               allocate(NinvA(rp%nchannel,rp%npolcomponent))
               allocate(At(rp%npolcomponent,rp%nchannel))
               allocate(cmap(rp%npolcomponent))
               allocate(smap(rp%npolcomponent))
               allocate(nmat(rp%npolcomponent,rp%npolcomponent))
               A(:,:) = tmpmat_reduced(:,:)
               At     = transpose(A)
            else
               allocate(A(rp%nchannel,rp%ncomponent))
               allocate(NinvA(rp%nchannel,rp%ncomponent))
               allocate(At(rp%ncomponent,rp%nchannel))
               allocate(cmap(rp%ncomponent))
               allocate(smap(rp%ncomponent))
               allocate(nmat(rp%ncomponent,rp%ncomponent))
               A(:,:) = tmpmat(:,:)
               At     = transpose(A)
            endif

            pix_index1 = 1
            pix_index2 = region(rr)%npix
            allocate(chi2(region(rr)%npix))
            !Getting ready for reduced memory mode.
!            pix_index1 = region(rr)%pix_index1   
!            pix_index2 = region(rr)%pix_index1 + region(rr)%npix - 1

            !-------------------
            !Set up noise matrix
            !-------------------
            do ff = 1, rp%nchannel     
               !Getting ready for reduced memory mode.
               !invdiagNoise(ff,ff) = MapWeightArray(ff)%TQU(pp,rp%stokes_list(ss))
               NinvA(ff,:)          = A(ff,:)*MapWeightArray(ff)%TQU(region(rr)%pix(1),rp%stokes_list(ss))
            end do
            
            !--------------------------------------------------
            ! Calculate once for all pixels (homegeneous noise)
            !--------------------------------------------------
            nmat     = matmul(At,NinvA)
            call Matrix_Inverse(nmat)

            do pp = pix_index1, pix_index2, rp%modulo_npix
               
               do ff = 1, rp%nchannel     
                  tmap(ff)            = MapDataWeightArray(ff)%TQU(region(rr)%pix(pp),rp%stokes_list(ss))
                  !Getting ready for reduced memory mode.
                  !tmap(ff)            = MapDataWeightArray(ff)%TQU(pp,rp%stokes_list(ss))
               end do
               
               !-----------------
               !Matrix operations
               !-----------------
               smap     = matmul(At,tmap)
               cmap     = matmul(nmat,smap)
               chi2(pp) = dot_product(smap,cmap)
            end do !Loop over pixels

            deallocate(A,At,smap,cmap,nmat,NinvA)
            do pp = pix_index1, pix_index2, rp%modulo_npix
               region(rr)%chi2 = region(rr)%chi2 + chi2(pp)
            enddo
            deallocate(chi2)
         end do !Loop over Stokes

         if(use_prior) then
            factor = real(rp%nstokes)
            do ii = firstparam_index,lastparam_index
               if (abs(scales%pvar(ii)) > 1d-10) then        
                  region(rr)%chi2 = region(rr)%chi2 - &
                       factor*(params%p(ii)-scales%pmean(ii))**2/scales%pvar(ii) - &
                       factor*dlog(scales%pvar(ii))
               endif
            end do
         end if

      end if
      sum_chi2 = sum_chi2 + region(rr)%chi2
   end do !Loop over region

   LnLike = -0.5d0*sum_chi2 + rp%lnlike_offset
   if (rp%Feedback>2) then
      print *, 'Parameters : ',Params%P(1:rp%nregion*rp%nparam_component)
      print *, 'Ln(like) = '  ,lnlike
   end if

 end function cmbforeground_LnLike5


 function cmbforeground_LnLike11(Params) result (LnLike)
   !-------------------------------------------------------------------
   !Multi-region spectral index likelihood, allowing for IQU covariance
   !-------------------------------------------------------------------

   !        
   !        ( A_11 A_12 A_13 )         (d_I1)                (N_11 0   0  )     
   !  A^t = ( A_21 A_22 A_23 ), tmap = (d_Q1), invdiagNoise =( 0  N_22 0  )
   !        ( A_31 A_32 A_33 )         (d_U1)                ( 0   0  N_33)
   !                                   (d_I2)
   !                                   (d_Q2)
   !                                   (d_U2)
   !                                   ( .. )
   !
   !   Each block of A^t is a diagonal matrix.
   !
   ! Need to set up the A matrix in a different way.
   ! A^t x tmap - time savings due to zeroes


   !   d = A s + n
   
   !         (s_I1)
   !         (s_Q1)
   !         (s_U1)
   !         (s_I2)
   !         (s_Q2)
   !         (s_U2)
   !         (s_I3)
   !         (s_Q3)
   !         (s_U4)
   

   
   implicit none
   type(ParamSet) Params
   real LnLike
   integer ff,pp,ss,rr,firstparam_index,lastparam_index,cc
   integer*8 pix_index1,pix_index2

   real(dp) tmap(rp%nchannel*rp%nstokes)
   real(dp) invdiagNoise(rp%nchannel*rp%nstokes,rp%nchannel*rp%nstokes)
   real(dp) pmat(rp%nchannel,rp%ncomponent)
   real(dp), allocatable, dimension(:,:) :: A   !(rp%nchannel,rp%npolcomponent) 
   real(dp), allocatable, dimension(:)   :: smap
   real(dp), allocatable, dimension(:,:) :: nmat
   real(dp), allocatable, dimension(:)   :: cmap
   real(dp) sum_chi2

   sum_chi2          = 0.d0
   LnLike            = 0.d0
   invdiagNoise(:,:) = 0.d0
   do rr = 1, rp%nregion !Loop over regions
      lastparam_index  = region(rr)%lastparam_index
      firstparam_index = region(rr)%firstparam_index

      if (all(abs(Params%P(firstparam_index:lastparam_index) - &
           LastParams%P(firstparam_index:lastparam_index)) < 1d-10)) then
         !---------------------------------------------------------------------------
         !.. the parameters, and hence the likelihood of this region have not changed
         !---------------------------------------------------------------------------
         continue
      else
         region(rr)%chi2 = 0.0d0
         !--------------------
         !Set up mixing matrix
         !--------------------
         pmat = get_mixingmatrix(region(rr)%component_list,&
              rp%component_reference_frequency_list,region(rr)%ncomponent, &
              Params%P(firstparam_index:lastparam_index),&
              rp%Inst%bandpass_ghz,rp%inst%bandpass_weight,rp%nchannel,&
              rp%Inst%nabscissa)         
         allocate(A(rp%nchannel*rp%nstokes,rp%ncomponent*rp%nstokes))
         A(:,:) = 0.
         do cc = 1, rp%ncomponent
            do ff = 1, rp%nchannel
               do ss = 1, rp%nstokes
                  A((ff-1)*rp%nstokes+ss,(cc-1)*rp%nstokes+ss) = pmat(ff,cc)
               end do
            end do
         end do

         allocate(cmap(rp%ncomponent*rp%nstokes))
         allocate(smap(rp%ncomponent*rp%nstokes))
         allocate(nmat(rp%ncomponent*rp%nstokes,rp%ncomponent*rp%nstokes))

         pix_index1        = 1
         pix_index2        = region(rr)%npix
         invDiagNoise(:,:) = 0.d0
         do pp = pix_index1, pix_index2, rp%modulo_npix !Loop over pixels
            !-------------------
            !Set up noise matrix
            !-------------------
            do ff = 1, rp%nchannel  !*rp%nstokes, rp%nstokes
               do ss = 1, rp%nstokes
                  tmap((ff-1)*rp%nstokes + ss)  = &
                       MapDataWeightArray(ff)%TQU(region(rr)%pix(pp), rp%stokes_list(ss))
                  !Diagonal part for now:
                  invdiagNoise((ff-1)*rp%nstokes + ss, (ff-1)*rp%nstokes + ss) = &
                       MapWeightArray(ff)%TQU(region(rr)%pix(pp),rp%stokes_list(ss)*2)
               end do
            end do
            !-----------------
            !Matrix operations
            !-----------------
            smap = matmul(transpose(A),tmap)
            nmat = matmul(transpose(A),matmul(invdiagNoise,A))
!            print *,tmap ! all zeroes
            call Matrix_Inverse(nmat)
            cmap = matmul(nmat,smap)
            region(rr)%chi2 = region(rr)%chi2 + dot_product(smap,cmap)
         end do
         deallocate(A,smap,cmap,nmat)
      end if
      sum_chi2 = sum_chi2 + region(rr)%chi2
   end do !Loop over region

   LnLike = -0.5d0*sum_chi2 + rp%lnlike_offset
   if (rp%Feedback>2) then
      print *, 'Parameters : ',Params%P(1:rp%nregion*rp%nparam_component)
      print *, 'Ln(like) = '  ,lnlike
   end if

 end function cmbforeground_LnLike11

 function cmbforeground_LnLike7(Params) result (LnLike)
   !---------------------------------------------------------------
   !Single region spectral index likelihood with calibration errors
   !---------------------------------------------------------------
   
   implicit none
   type(ParamSet) Params
   real LnLike
   integer ff,pp,ss,index1,index2,ii

   real(dp) tmap(rp%nchannel)
   real(dp) invdiagNoise(rp%nchannel,rp%nchannel)
   real(dp) NinvOmegaA(rp%nchannel,rp%nchannel)
   real(dp) A(rp%nchannel,rp%ncomponent)
   real(dp) OmegaA(rp%nchannel,rp%ncomponent)
   real(dp) OmegaAt(rp%ncomponent,rp%nchannel)
   real(dp) smap(rp%ncomponent)
   real(dp) nmat(rp%ncomponent,rp%ncomponent)
   real(dp) cmap(rp%ncomponent)
   real(dp) sum_chi2,calib_penalty
   real(dp) Omega(rp%nchannel,rp%nchannel)
   real(dp) omegavec(rp%nchannel)
   real(dp) InvSigma(rp%nchannel,rp%nchannel)
   real(dp) omegabar(rp%nchannel)
   real(dp) factor

   !--------------------
   !Set up mixing matrix
   !--------------------
   A = get_mixingmatrix(rp%component_list,&
        rp%component_reference_frequency_list,rp%ncomponent, &
        Params%P(1:rp%nparam_component),&
        rp%Inst%bandpass_ghz,rp%inst%bandpass_weight,rp%nchannel,&
        rp%Inst%nabscissa)
   
   if (rp%calibration_marginalisation) then
      !--------------------------------------
      !Set up calibration data and parameters
      !--------------------------------------
      InvSigma(:,:) = 0.d0
      Omega(:,:)    = 0.d0
      do ff = 1, rp%nchannel
         omegabar(ff)    = 1.
         Omega(ff,ff)    = Params%P(rp%nparam_component+ff)
         omegavec(ff)    = Params%P(rp%nparam_component+ff)
         InvSigma(ff,ff) = 1./rp%calib_rms(ff)**2
      end do
      calib_penalty = dot_product(matmul(omegavec-omegabar,InvSigma),omegavec-omegabar)
      OmegaA        = matmul(Omega,A)
   else      
      calib_penalty = 0.0
      OmegaA        = A
   endif
   OmegaAt          = transpose(OmegaA)

   !----------------
   !Loop over pixels
   !----------------
   LnLike            = 0.d0
   sum_chi2          = 0.d0
   invdiagNoise(:,:) = 0.d0
   do ss = 1, rp%nstokes
      do pp = 1, Hmask%npix 
         !-------------------
         !Set up noise matrix
         !-------------------
         do ff = 1, rp%nchannel     
            tmap(ff)         = MapDataWeightArray(ff)%TQU(Hmask%good(pp),rp%stokes_list(ss))                  
            NinvOmegaA(ff,:) = OmegaA(ff,:)*MapWeightArray(ff)%TQU(Hmask%good(pp),rp%stokes_list(ss))
         end do
         
         !-----------------
         !Matrix operations
         !-----------------
         smap     = matmul(OmegaAt,tmap)
         nmat     = matmul(OmegaAt,NinvOmegaA)
         call Matrix_Inverse(nmat)
         cmap     = matmul(nmat,smap)
         sum_chi2 = sum_chi2 + dot_product(smap,cmap)
      end do
   end do

   !----------------------------
   ! Apply spectral index priors
   !----------------------------
   if(use_prior) then
      factor = real(rp%nstokes) 
      do ii = 1,rp%nparam_component
         if (abs(scales%pvar(ii)) > 1d-10) then        
            sum_chi2 = sum_chi2 - &
                 factor*(params%p(ii)-scales%pmean(ii))**2/scales%pvar(ii) - &
                 factor*dlog(scales%pvar(ii))
         endif
      end do
   end if


   LnLike = -0.5*(sum_chi2 - calib_penalty) + rp%lnlike_offset
   if (rp%Feedback > 2) then
      print *, 'Foreground parameters : ',Params%P(1:rp%nregion*rp%nparam_component)
      if (rp%calibration_marginalisation) then
         index1 = rp%nregion*rp%nparam_component + 1
         index2 = index1 - 1 + rp%nchannel
         print *, 'Calibrations parameters : ',Params%P(index1:index2)
      endif
      print *, 'Ln(like) = ',lnlike
   end if

 end function cmbforeground_LnLike7


 function cmbforeground_LnLike8(Params) result (LnLike)
   !-------------------------------------------------------------------------------------
   !Single region spectral index likelihood with calibration errors and homogeneous noise
   !-------------------------------------------------------------------------------------
   
   implicit none
   type(ParamSet) Params
   real LnLike
   integer ff,pp,ss,index1,index2,ii

   real(dp) tmap(rp%nchannel)
   real(dp) invdiagNoise(rp%nchannel,rp%nchannel)
   real(dp) NinvOmegaA(rp%nchannel,rp%nchannel)
   real(dp) A(rp%nchannel,rp%ncomponent)
   real(dp) OmegaA(rp%nchannel,rp%ncomponent)
   real(dp) OmegaAt(rp%ncomponent,rp%nchannel)
   real(dp) smap(rp%ncomponent)
   real(dp) nmat(rp%ncomponent,rp%ncomponent)
   real(dp) cmap(rp%ncomponent)
   real(dp) sum_chi2,calib_penalty
   real(dp) Omega(rp%nchannel,rp%nchannel)
   real(dp) omegavec(rp%nchannel)
   real(dp) InvSigma(rp%nchannel,rp%nchannel)
   real(dp) omegabar(rp%nchannel)
   real(dp) factor

   !--------------------
   !Set up mixing matrix
   !--------------------
   A = get_mixingmatrix(rp%component_list,&
        rp%component_reference_frequency_list,rp%ncomponent, &
        Params%P(1:rp%nparam_component),&
        rp%Inst%bandpass_ghz,rp%inst%bandpass_weight,rp%nchannel,&
        rp%Inst%nabscissa)
   
   if (rp%calibration_marginalisation) then
      !--------------------------------------
      !Set up calibration data and parameters
      !--------------------------------------
      InvSigma(:,:) = 0.d0
      Omega(:,:)    = 0.d0
      do ff = 1, rp%nchannel
         omegabar(ff)    = 1.
         Omega(ff,ff)    = Params%P(rp%nparam_component+ff)
         omegavec(ff)    = Params%P(rp%nparam_component+ff)
         InvSigma(ff,ff) = 1./rp%calib_rms(ff)**2
      end do
      calib_penalty = dot_product(matmul(omegavec-omegabar,InvSigma),omegavec-omegabar)
      OmegaA        = matmul(Omega,A)
   else      
      calib_penalty = 0.0
      OmegaA        = A
   endif
   OmegaAt          = transpose(OmegaA)

   !----------------
   !Loop over pixels
   !----------------
   LnLike            = 0.d0
   sum_chi2          = 0.d0
   invdiagNoise(:,:) = 0.d0
   do ss = 1, rp%nstokes
      !-------------------
      !Set up noise matrix
      !-------------------
      do ff = 1, rp%nchannel     
         NinvOmegaA(ff,:) = OmegaA(ff,:)*MapWeightArray(ff)%TQU(Hmask%good(1),rp%stokes_list(ss))
      end do
      nmat     = matmul(OmegaAt,NinvOmegaA)
      call Matrix_Inverse(nmat)

      do pp = 1, Hmask%npix 
         !-----------
         !Set up data
         !-----------
         do ff = 1, rp%nchannel     
            tmap(ff)         = MapDataWeightArray(ff)%TQU(Hmask%good(pp),rp%stokes_list(ss))                  
         end do
         
         !-----------------
         !Matrix operations
         !-----------------
         smap     = matmul(OmegaAt,tmap)
         cmap     = matmul(nmat,smap)
         sum_chi2 = sum_chi2 + dot_product(smap,cmap)
      end do
   end do

   !----------------------------
   ! Apply spectral index priors
   !----------------------------
   if(use_prior) then
      factor = real(rp%nstokes) 
      do ii = 1,rp%nparam_component
         if (abs(scales%pvar(ii)) > 1d-10) then        
            sum_chi2 = sum_chi2 - &
                 factor*(params%p(ii)-scales%pmean(ii))**2/scales%pvar(ii) - &
                 factor*dlog(scales%pvar(ii))
         endif
      end do
   end if


   LnLike = -0.5*(sum_chi2 - calib_penalty) + rp%lnlike_offset
   if (rp%Feedback > 2) then
      print *, 'Foreground parameters : ',Params%P(1:rp%nregion*rp%nparam_component)
      if (rp%calibration_marginalisation) then
         index1 = rp%nregion*rp%nparam_component + 1
         index2 = index1 - 1 + rp%nchannel
         print *, 'Calibrations parameters : ',Params%P(index1:index2)
      endif
      print *, 'Ln(like) = ',lnlike
   end if

 end function cmbforeground_LnLike8


 function cmbforeground_LnLike10(Params) result (LnLike)
   !------------------------------------------------------------------------------------------
   !Single region spectral index likelihood with calibration errors and offset marginalisation
   !------------------------------------------------------------------------------------------

   implicit none
   type(ParamSet) Params
   real LnLike
   integer ff,pp,ss,stokes_index,status,nfrq

   real(dp) invdiagNoise(rp%nchannel,rp%nchannel)
   real(dp) pmat(rp%nchannel,rp%ncomponent)
   real(dp) tmpmat(rp%nchannel,rp%ncomponent)
   real(dp) invNtmpmat(rp%nchannel,rp%ncomponent)
   real(dp) sum_chi2,calib_penalty
   real(dp) Omega(rp%nchannel,rp%nchannel)
   real(dp) omegavec(rp%nchannel)
   real(dp) InvSigma(rp%nchannel,rp%nchannel)
   real(dp) omegabar(rp%nchannel)

   !Offset stuff - precomputable    
   real(dp) invuNu(rp%nchannel,rp%nchannel)! (U^t N^-1 U)^-1   
   real(dp) compMap(rp%nchannel)! (U^t N^-1 U)^-1 U^t N^-1 d : compressed map (projected onto the offset space)   
   real(dp) const2 ! = 1d0 :\gamma^2 for stabilizing the matrix inversion
   !Offset stuff - on the fly
   real(dp) U(rp%nchannel)!,V(rp%nchannel)
   integer nsing  ! A number of singular modes of M^{-1} and A^t M^{-1} A
   real(dp), dimension(rp%ncomponent,rp%nchannel) ::  singVecttmp ! dimensions (nsing, ncomp)
   real(dp), dimension(:,:), allocatable :: singVect ! dimensions (nsing, ncomp)
   real(dp) QQterm(rp%nchannel,rp%nchannel)
   real(dp) XXterm(rp%ncomponent,rp%ncomponent)
   real(dp) XQterm(rp%ncomponent,rp%nchannel)
   real(dp), dimension(:), allocatable :: qVect ! dimensions (nsing+nchan)
   real(dp), dimension(:), allocatable :: pVect ! dimensions (nsing+nchan)
   real(dp) tmpVect(rp%nchannel)
   real(dp) aMd(rp%ncomponent)
   real(dp) invaNa(rp%ncomponent,rp%ncomponent)
   real(dp) aN(rp%ncomponent,rp%nchannel) !??
   real(dp) invRMat(rp%nchannel*2,rp%nchannel*2)
   real(dp), dimension(:,:), allocatable :: rKernel ! dimensions (nsing+nchan,nsing+nchan)
   real(dp)  tmpMap(rp%ncomponent)! 

   nfrq   = rp%nchannel
   const2 = 1.d0

   !--------------------
   !Set up mixing matrix
   !--------------------
   pmat = get_mixingmatrix(rp%component_list,&
        rp%component_reference_frequency_list,rp%ncomponent, &
        Params%P(1:rp%nparam_component),&
        rp%Inst%bandpass_ghz,rp%inst%bandpass_weight,rp%nchannel,&
        rp%Inst%nabscissa)

   if (rp%calibration_marginalisation) then
      !Set up calibration data and parameters
      InvSigma(:,:) = 0.0d0
      Omega(:,:)    = 0.0d0
      do ff=1, nfrq
         omegabar(ff)    = 1d0
         Omega(ff,ff)    = Params%P(rp%nparam_component+ff)
         omegavec(ff)    = Params%P(rp%nparam_component+ff)
         InvSigma(ff,ff) = 1d0/rp%calib_rms(ff)**2
      end do
      calib_penalty = dot_product(matmul(omegavec-omegabar,InvSigma),omegavec-omegabar)
      tmpmat        = matmul(Omega,pmat)
   else
      calib_penalty = 0.0d0
      tmpmat        = pmat
   endif

   LnLike   = 0.d0
   sum_chi2 = 0.d0
   do ss=1,rp%nstokes
     stokes_index=ss
      call offsetmarg_computations(tmpmat,nfrq,Hmask%npix, &
           rp%ncomponent,stokes_index,SingVecttmp,U,nsing,compMap,invuNu)
      allocate(singvect(rp%ncomponent,nsing))
      singvect(:,:)= SingVecttmp(1:rp%ncomponent,1:nsing)
      QQterm(:,:) = 0.d0
      XXterm(:,:) = 0.d0
      XQterm(:,:) = 0.d0
      allocate(qVect(nsing+nfrq))
      allocate(pVect(nsing+nfrq))
      qVect(:)          = 0.d0
      pVect(:)          = 0.d0      
      invdiagNoise(:,:) = 0.d0
      !Loop over pixels
      do pp = 1, Hmask%npix 
         ! Compute A^t M^-1 d for this pix and set up noise matrix

         do ff = 1, nfrq     
            tmpVect(ff)         =(MapDataArray(ff)%TQU(Hmask%good(pp),rp%stokes_list(ss))-compMap(ff))*&
                 MapWeightArray(ff)%TQU(Hmask%good(pp),rp%stokes_list(ss))
            invdiagNoise(ff,ff) = MapWeightArray(ff)%TQU(Hmask%good(pp),rp%stokes_list(ss))
         end do
                  
         aMd        = matmul(transpose(tmPmat),tmpVect)
         ! Compute p-p block of (A^t N^-1 A)^-1
         invNtmpmat = matmul(invdiagNoise,tmPmat)
         invaNa     = matmul(transpose(tmPmat),invNtmpmat)
         call Matrix_Inverse(invaNa)

         ! Precompute (A^t N^-1 A)^-1 A^t M^-1 d
         tmpMap = matmul(invANA, aMd)

         ! U^t (A^t M^-1 d)^t (A^t N^-1 A)^-1 (A^t M^-1 d) U for each pix
         sum_chi2 = sum_chi2 + dot_product(aMd,tmpMap) ! = transpose( aMd) ## fmap[ p, 0:ncomp-1]

         ! (A^t N^-1 U)^t (A^t N^-1 A)^-1 A^t M^-1 d
         ! and X^t (A^t N^-1 A)^-1 A^t M^-1 d
         qVect(1:nfrq)            = qVect(1:nfrq) - &
              matmul(invdiagnoise,matmul(tmPmat,tmpMap))
         qVect(nfrq+1:nfrq+nsing) = qVect(nfrq+1:nfrq+nsing) + &
              matmul(transpose(singVect),tmpMap)

         ! U^t N^-1 A (A^t N^-1 A)^-1 A^t N^-1 U
         aN     = matmul(transpose(tmPmat),invdiagnoise)
         QQterm = QQterm+ matmul(transpose(aN),matmul(invANA,aN))

         ! (A^t N^-1 A)^-1 A^t N^-1 U
         XQterm = XQterm + matmul(invANA,aN)

         ! (A^t N^-1 A)^-1
         XXterm = XXterm +invANA
      end do !Over pixels

      ! Compute the R kernel
      invrMat(:,:)           = 0.d0
      invrMat(1:nfrq,1:nfrq) = invuNu(1:nfrq,1:nfrq)
      do ff=1, nfrq
         invrMat(nfrq+ff,nfrq+ff)=const2
      end do
      call Matrix_Inverse(invRmat)

      allocate(rKernel(nfrq+nsing,nfrq+nsing))
      Rkernel                           = 0.0d0
      Rkernel(1:nfrq,1:nfrq)            = -QQterm 
      Rkernel(1:nfrq,nfrq+1:nfrq+nsing) = -matmul(transpose(XQterm),singvect)!CHECK
      Rkernel(nfrq+1:nfrq+nsing,1:nfrq) =  matmul(transpose(singvect),XQterm)
      Rkernel(nfrq+1:nfrq+nsing,nfrq+1:nfrq+nsing) = &
           matmul(transpose(singvect),matmul(XXterm,singvect))

      Rkernel(1:nfrq+nsing,1:nfrq+nsing) = &
           Rkernel(1:nfrq+nsing,1:nfrq+nsing) + &
           invRmat(1:nfrq+nsing,1:nfrq+nsing)
      call Matrix_Inverse_new(rKernel) !!NEED TO FIX THIS??

      ! And now the relavant SMW correction
      pVect(1:nfrq)            = -qVect(1:nfrq)
      pVect(nfrq+1:nfrq+nsing) =  qVect(nfrq+1:nfrq+nsing)
      sum_chi2                 = sum_chi2 - dot_product(pVect,matmul(rKernel,pVect))

      if (nsing .gt. 0) then
        deallocate(singVect,stat=status)
      endif
      deallocate(pVect,stat=status)
      deallocate(qVect,stat=status)
      deallocate(rKernel,stat=status)
   end do !Loop over stokes parameters

   LnLike = -0.5*(sum_chi2 - calib_penalty) + rp%lnlike_offset
   print *, Lnlike,params%p(1), params%p(2), params%p(3), params%p(4)

 end function cmbforeground_LnLike10

 subroutine offsetmarg_computations(pmat,nchannel, &
      npix,ncomponent,stokes_index,SingVect,U,nsing,compMap,invuNu)

  !-----
  !Input                                    
  !-----
  integer*8, intent(inout) :: npix
  integer,   intent(inout) :: nchannel,ncomponent,stokes_index
  real(dp),  intent(inout) :: pmat(nchannel,ncomponent)

  !------
  !Output                                    
  !------
  real(dp), intent(out) :: singVect(ncomponent,nchannel)  ! dimensions (nsing, ncomp)
  real(dp), intent(out) :: U(nchannel)
  integer,  intent(out) :: nsing
  real(dp), intent(out) :: compMap(nchannel)! (U^t N^-1 U)^-1 U^t N^-1 d : compressed map (projected onto $
  real(dp), intent(out) :: invuNu(nchannel,nchannel)! (U^t N^-1 U)^-1
  
  !Internal
  real(dp), dimension(:,:), allocatable :: sVect 
  real(dp) V(nchannel)
  integer ff,pp,status
  real(dp) uNd(nchannel) ! U^t N^-1 d
  real(dp) invaTa(ncomponent,ncomponent)
  real(dp) projmat(nchannel,nchannel)
  real(dp) W(nchannel,nchannel)
  integer indxsing(nchannel) !Which eigenmodes of projmat are singular
  real(dp) sv_sum
  real(dp), dimension(:,:), allocatable :: wsing ! dimensions (nsing, nfreq)
  
  !U^t N^-1 U kernel and U^t N^-1 d
  invuNu(:,:) = 0.d0
  uNd(:)      = 0.d0
  do ff = 1, nchannel
     do pp = 1, npix
        invuNu(ff,ff) = invuNu(ff,ff)+MapWeightArray(ff)%TQU(Hmask%good(pp),rp%stokes_list(stokes_index))
        uNd(ff)       = uNd(ff)+ MapdataArray(ff)%TQU(Hmask%good(pp),rp%stokes_list(stokes_index))* &
             MapWeightArray(ff)%TQU(Hmask%good(pp),rp%stokes_list(stokes_index))
     end do
  end do

  !(U^t N^-1 U)^-1
  call Matrix_Inverse(invuNu)
  
  !(U^t N^-1 U)^-1 U^t N^-1 d
  compMap = matmul(invunu, uNd)
  
  ! Compute singular vectors (identical for all pixs)
  invaTa = matmul(transpose(pmat),pmat)
  call Matrix_Inverse(invaTa)
  
  projmat= matmul(pmat,matmul(invaTa,transpose(pmat)))
  
  W=projmat
  call Matrix_SVD(W,nchannel,nchannel,U,V)
  
  !Work out which eigen modes are singular
  nsing=0
  do ff=1, nchannel
     if( abs(U(ff)-1d0) .lt. 1d-10) then
        nsing           = nsing+1
        indxsing(nsing) = ff
     endif
  end do
  
  if (nsing .gt. 0) then
     allocate(wsing(nchannel,nsing))
     allocate(sVect(ncomponent,nsing))
     wsing = W(1:nchannel,indxsing(1:nsing))
     sVect = matmul(invAta,matmul(transpose(pmat),wsing))
  endif

  do ff=1,nsing
    sv_sum         = dot_product(sVect(:,ff),sVect(:,ff))
    singVect(:,ff) = sVect(:,ff)/sqrt(sv_sum)
  end do
  
  if (nsing .gt. 0) then
     deallocate(wsing,stat=status)
     deallocate(svect)
  endif

end subroutine offsetmarg_computations


  

end module Likelihoodcmbforeground
