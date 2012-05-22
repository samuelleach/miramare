module RunParameters
!Module for reading and storing the parameters for component separation.
  use IniFile !Reading parameters
  use paramdef !Contains DoStop() subroutine
  use mpi_stop
  use Instrumentdata
  use foregroundscalings
  use AMLutils
  use miramare_utils

  type RunParams !Parameters associated with component separation and I/O.
     character(LEN=80)  :: rootname
     character(LEN=256) :: cachedir !Directory for cached simulations data.
     character(LEN=256) :: plotdir  !Directory for plots from miramare_getdist.
     character(LEN=256) :: out_dir  !Directory for all other output files from miramare_estamp
     character(LEN=256) :: maskfile, regionfile,pixelsofinterestfile    
     character(LEN=256) :: cls_file !Cl of CMB for unlensed simulation
     character(3)       :: stokes   !IQU, QU, I, Q or U    
     integer nregion,feedback,nstokes,nchannel
     integer ntemplate,likelihood_label,nside_in,ordering_in
     logical calibration_marginalisation,offset_marginalisation,remove_raw_offsets
     type(Instrument) Inst
     real, dimension(:), allocatable   :: calib, offset,calib_rms,fwhm_arcmin,sim_calib,sim_calib_rms,sim_offset
     real, dimension(:,:), allocatable :: template_coefficient
     real, dimension(:), allocatable   :: stokesi_noise_sig0, stokesqu_noise_sig0
     character(LEN=256), dimension(:), allocatable :: datafile,weightfile,templatefile
     real*8 lnlike_offset
     logical noise_matrix_correction,want_pol,want_plots
     logical output_whitened_residuals
     logical output_relative_residuals
     logical sim_signal,sim_noise,sim_caliberror,sim_suppress_bmode
     logical sim_component,sim_component_perturbations
     integer npol
     integer weight_format, mcmc_rand_seed, noise_rand_seed
     integer signal_rand_seed, component_rand_seed, caliberror_rand_seed
     integer, dimension(:), allocatable :: stokes_list
     integer modulo_npix, sim_lmax, sim_ncomponent
     character(32), dimension(:), allocatable  :: sim_component_list
     character(256), dimension(:), allocatable :: sim_componentfile
     real(dp), dimension(:), allocatable       :: sim_componentfile_factor
     real(dp), dimension(:,:), allocatable     :: sim_component_param
     real(dp), dimension(:), allocatable       :: sim_component_frequency
     real(dp), dimension(:), allocatable       :: sim_noise_factor
     real(dp), dimension(:), allocatable       :: weight_factor
     character(LEN=256), dimension(5,5) :: sim_component_cls_file                !Cl of component perturbations
     character(LEN=256), dimension(5,5) :: sim_component_paramvariationsmap_file !Map of component perturbations
     integer npixtot

     !Move these parameteters to regional model parameters?
     integer ncomponent,nparam_component,npolcomponent
     character(32), dimension(:), allocatable  :: component_list
     real(dp), dimension(:), allocatable       :: component_reference_frequency_list
     integer, dimension(:), allocatable        :: polcomponent_list
     logical, dimension(:), allocatable        :: component_is_pol
     character(256), dimension(:), allocatable :: tabulated_scaling_component

     integer do_compsep_on_simulated_component, npixel_of_interest
     integer, dimension(:), allocatable       :: pixel_of_interest
     logical homogeneous_noise

     logical convert_units_cmb_to_rj

   end type RunParams
   
   !---------------------
   !GLOBAL run parameters
   !---------------------
   type(RunParams) rp

contains
  
  subroutine get_runparams(parfile,rp)
    character(LEN=Ini_max_string_len) parfile
    Type(RunParams) rp
    logical bad
    integer cc,pp,ff,tt,ss,index,max_params
    character(32) component

    integer nabscissa,ii,numpoint
    real*8 deltaf
    real(dp), pointer, dimension(:,:) :: spectrum


    nabscissa = 5  !Bandpasses 

    call Ini_Open(parfile, 2, bad, .false.)
    if (bad) call DoStop('Error opening parameter file')

    rp%feedback = Ini_Read_Int('feedback',2)
    rp%rootname = Ini_Read_String('file_root')
    rp%rootname = ExtractFileName(rp%rootname)

    if (rp%feedback .gt. 0) then
       print *, '-----------------------'
       print *, 'Reading parameters from ',trim(parfile)
       print *, '-----------------------'
    endif

    rp%out_dir = Ini_Read_String('out_dir')
    if (rp%out_dir /= '') then
       rp%out_dir = CheckTrailingSlash(rp%out_dir)
    else
       rp%out_dir = './'
    endif
    if(.not. directory_exists(rp%out_dir)) then
       print *, 'ERROR: Directory does not exist: out_dir = ',trim(rp%out_dir)
       call mpistop('ERROR: Please mkdir this directory or change out_dir parameter.')
    endif


    if (rp%feedback > 0 ) write (*,*) ' Producing output files in directory '//&
         trim(rp%out_dir)

    rp%cachedir = Ini_Read_String('cachedir')
    if (rp%cachedir .eq. '') rp%cachedir = rp%out_dir
    if(.not. directory_exists(rp%cachedir)) then
       print *, 'ERROR: Directory does not exist: cachedir ',trim(rp%cachedir)
       call mpistop('ERROR: Please mkdir this directory or change cachedir parameter.')
    endif

    rp%plotdir                         = Ini_Read_String('plotdir')
    if (rp%plotdir .eq. '') rp%plotdir = rp%out_dir
    if(.not. directory_exists(rp%plotdir)) then
       print *, 'ERROR: Directory does not exist: ',trim(rp%plotdir)
       call mpistop('ERROR: Please mkdir this directory.')
    endif


    rp%want_plots = Ini_Read_Logical('want_plots',.false.)
    rp%maskfile   = Ini_Read_String('maskfile')
    if (rp%maskfile .eq. '') then
       rp%nside_in = Ini_Read_Int('nside_in',-1)
    endif    
    rp%regionfile                            = Ini_Read_String('regionfile')
    if (rp%regionfile .eq. '') rp%regionfile = rp%maskfile
    rp%stokes  = Ini_Read_String('stokes')
    rp%nstokes = stokes2nstokes(rp%stokes)

    allocate(rp%stokes_list(rp%nstokes))
    do ss = 1, rp%nstokes
       rp%stokes_list(ss) = stokes_healpix_index(rp%stokes,ss)
    end do
    
    rp%ncomponent = Ini_Read_Int('ncomponent',2)
    allocate(rp%component_list(rp%ncomponent))
    allocate(rp%tabulated_scaling_component(rp%ncomponent))
    allocate(rp%component_reference_frequency_list(rp%ncomponent))
    allocate(rp%component_is_pol(rp%ncomponent))
    
    !--------------------------------------------------
    ! Allocate global data in foregroundscalings module
    !--------------------------------------------------
    allocate(tabulated_scaling(rp%ncomponent,nmaxscalingpoints,2))
    allocate(tabulated_scaling_numpoint(rp%ncomponent))

    rp%npolcomponent = 0
    do cc = 1, rp%ncomponent
       rp%component_list(cc)                     = Ini_Read_String_array('component',cc)
       rp%component_reference_frequency_list(cc) = Ini_Read_Real_array('component_reference_frequency',cc,150.)
       rp%component_is_pol(cc)                   = Ini_Read_Logical_array('component_is_pol',cc,.true.)
       if (rp%component_is_pol(cc))  rp%npolcomponent = rp%npolcomponent + 1             
       if (rp%Feedback>0) then 
          if (rp%component_is_pol(cc)) then
             print *,' Component ',cc,' = '//trim(rp%component_list(cc))//' (polarized)'
          else
             print *,' Component ',cc,' = '//trim(rp%component_list(cc))//' (unpolarized)'
          endif
          if (rp%component_list(cc) .ne. 'cmb') then
             print *,' Reference frequency (GHz)  = ',rp%component_reference_frequency_list(cc)
          endif
!          else
!             print *,' CMB component is in thermodynamic units.'
!          endif
          if (rp%component_list(cc) .eq. 'tabulated') then
             rp%tabulated_scaling_component(cc) = Ini_Read_String_array('tabulated_scaling_component',cc)
             print *,' Reading tabulated component: '//trim(rp%tabulated_scaling_component(cc))
             call read_spectrum(rp%tabulated_scaling_component(cc),&
                  spectrum,rp%component_reference_frequency_list(cc),numpoint)
             do pp = 1, numpoint
                tabulated_scaling(cc,pp,1) = spectrum(pp,1)
                tabulated_scaling(cc,pp,2) = spectrum(pp,2)
             end do
             tabulated_scaling_numpoint(cc) = numpoint
             deallocate(spectrum)
          endif
       end if
    end do
    allocate(rp%polcomponent_list(rp%npolcomponent))
    index = 0
    do cc = 1, rp%ncomponent
       if (rp%component_is_pol(cc)) then
          index                       = index + 1
          rp%polcomponent_list(index) = cc
       endif
    end do

    rp%nparam_component            = get_totalforegroundparameters(rp%component_list,rp%ncomponent)
    rp%nchannel                    = Ini_Read_Int('nchannel',3)
    rp%calibration_marginalisation = Ini_Read_Logical('calibration_marginalisation',.true.)
    rp%offset_marginalisation      = Ini_Read_Logical('offset_marginalisation',.false.)
    call Instrument_init(rp%Inst,rp%nchannel,rp%nstokes,nabscissa)
    allocate(rp%calib(rp%nchannel),rp%offset(rp%nchannel),rp%calib_rms(rp%nchannel))
    allocate(rp%datafile(rp%nchannel),rp%fwhm_arcmin(rp%nchannel))
    do ff=1, rp%nchannel
       rp%calib(ff)    = Ini_Read_Real_array('calib',ff,1.)
       rp%offset(ff)   = Ini_Read_Real_array('offset',ff,0.)
       rp%datafile(ff) = Ini_Read_String_array('datafile',ff)
       if(rp%calibration_marginalisation) rp%calib_rms(ff) = Ini_Read_Real_array('calib_rms',ff,0.02)
       rp%Inst%channel_nu0_ghz(ff) = Ini_Read_Real_array('frequency',ff,-1.)
       if (rp%Inst%channel_nu0_ghz(ff) .eq. -1.) then
          print *,'runparameters.f90: ERROR: Set frequency() of channel '//trim(inttostr(ff))//' in parameter file.'
          stop
       endif
       rp%fwhm_arcmin(ff) = Ini_Read_Real_array('fwhm_arcmin',ff,0.)
    end do
    rp%remove_raw_offsets = Ini_Read_Logical('remove_raw_offsets',.false.)
    rp%lnlike_offset      = Ini_Read_Double('lnlike_offset',0.d0)
    
    rp%ntemplate = Ini_Read_Int('ntemplate',0)
    allocate(rp%templatefile(rp%ntemplate))
    allocate(rp%template_coefficient(rp%ntemplate,rp%nchannel))
    if (rp%ntemplate .gt. 0) then
       do tt = 1, rp%ntemplate
          rp%templatefile(tt) = Ini_Read_String_array('templatefile',tt)
          do ff = 1, rp%nchannel
             rp%template_coefficient(tt,ff) =& 
                  Ini_Read_Real_array('template'//trim(inttostr(tt))//'_coefficient',ff,0.)
          enddo
       enddo
    endif

    rp%noise_matrix_correction = Ini_Read_Logical('noise_matrix_correction',.false.)

    do ff = 1, rp%nchannel
       rp%Inst%bandpass_ghz(1,ff) = &
            Ini_Read_Double_array('tophatmin_frequency',ff,rp%Inst%channel_nu0_ghz(ff))
       rp%Inst%bandpass_ghz(nabscissa,ff) = &
            Ini_Read_Double_array('tophatmax_frequency',ff,rp%Inst%channel_nu0_ghz(ff))
       print '(a,I3,X,F8.3,X,F8.3)','  Channel, bandpass lower, upper edge [GHz] ',&
            ff,rp%Inst%bandpass_ghz(1,ff),&
            rp%Inst%bandpass_ghz(nabscissa,ff)
       deltaf = (rp%Inst%bandpass_ghz(nabscissa,ff)-rp%Inst%bandpass_ghz(1,ff))/&
            float(nabscissa-1)
       !-------------------------
       ! Set up tophat bandpasses
       !-------------------------
       do ii = 1, nabscissa
          rp%Inst%bandpass_ghz(ii,ff)    = rp%Inst%bandpass_ghz(1,ff) + &
               float(ii-1)*deltaf
          rp%Inst%bandpass_weight(ii,ff) = 1./float(nabscissa)
       enddo
    enddo
    
    if (rp%nstokes .gt. 1 .or. rp%stokes_list(1) .gt. 1) then
       rp%want_pol = .true.
       rp%npol     = 3
    else
       rp%want_pol = .false.
       rp%npol     = 1
    endif
    
    rp%weight_format = Ini_Read_Int('weight_format',1)
!    weight_format = 1 => Read IQU weight maps (units 1/data^2)
!    weight_format = 2 => Read hit maps + I and QU pixel RMS.
!    weight_format = 3 => Read "WMAP style" hits in column 4 of data + I and QU pixel RMS.
!    weight_format = 4 => Homogeneous noise => I and QU pixel RMS.
!    weight_format = 5 => "Madam" style 3x3 IQU inverse noise covariance.

    if (rp%weight_format .eq. 1 .or. rp%weight_format .eq. 2 .or. rp%weight_format .eq. 5) then
       allocate(rp%weightfile(rp%nchannel))
       do ff = 1, rp%nchannel
          rp%weightfile(ff) = Ini_Read_String_array('weightfile',ff)
       enddo
    endif
    if (rp%weight_format .eq. 2 .or. rp%weight_format .eq. 3&
         .or. rp%weight_format .eq. 4) then
       allocate(rp%stokesi_noise_sig0(rp%nchannel))
       do ff = 1, rp%nchannel
          rp%stokesi_noise_sig0(ff) = Ini_Read_real_array('stokesi_noise_sig0',ff,0.)
          if(rp%stokesi_noise_sig0(ff) .eq. 0.0) then
             print *, ' Warning: stokesi_noise_sig0 = 0.0'
          endif
       enddo
       
       if (rp%want_pol) then
          allocate(rp%stokesqu_noise_sig0(rp%nchannel))
          do ff = 1, rp%nchannel
             rp%stokesqu_noise_sig0(ff) = Ini_Read_real_array('stokesqu_noise_sig0',ff,0.)
             if(rp%stokesqu_noise_sig0(ff) .eq. 0.0) then
                print *, ' Warning: stokesqu_noise_sig0 = 0.0'
                stop
             endif
          enddo
       endif
    endif
    if (rp%weight_format .eq. 4) then
       rp%homogeneous_noise = .true.
    else
       rp%homogeneous_noise = .false.
    end if


    !-------------------------------------------------------
    ! Set up the likelihood label to be used in calclike.f90
    !-------------------------------------------------------
    if((rp%offset_marginalisation .eqv. .false.) .and. (rp%calibration_marginalisation .eqv. .false.)) then
       select case(rp%homogeneous_noise)
       case(.true.)
          rp%likelihood_label = 5
          if (rp%feedback > 0 ) write (*,*) &
               ' Using likelihood function without calibration or offset marginalisation'
          if (rp%feedback > 0 ) write (*,*) &
               ' and assuming homogeneous noise'
       case(.false.)
          rp%likelihood_label = 4
          if (rp%feedback > 0 ) write (*,*) &
               ' Using likelihood function without calibration or offset marginalisation'
       end select
    endif
    if((rp%offset_marginalisation .eqv. .false.) .and. (rp%calibration_marginalisation .eqv. .true.)) then
       select case(rp%homogeneous_noise)
       case(.true.)
          rp%likelihood_label = 8 
          if (rp%feedback > 0 ) write (*,*) ' Using likelihood function with calibration but not offset marginalisation.'
          if (rp%feedback > 0 ) write (*,*) ' Assuming homogeneous noise.'
       case(.false.)
          rp%likelihood_label = 7 
          if (rp%feedback > 0 ) write (*,*) ' Using likelihood function with calibration but not offset marginalisation'
       end select
    endif
    if(rp%offset_marginalisation .eqv. .true.) then
       rp%likelihood_label = 9 
       if (rp%feedback > 0 ) write (*,*) ' Using likelihood function with calibration and offset marginalisation'
       if (rp%feedback > 0 ) write (*,*) ' WARNING: this likelihood function can only be run for single region data'
    endif

    if(rp%weight_format .eq. 5) then
       rp%likelihood_label = 11
    end if

    rp%mcmc_rand_seed       = Ini_Read_Int('mcmc_rand_seed',1)
    rp%noise_rand_seed      = Ini_Read_Int('noise_rand_seed',100)
    rp%signal_rand_seed     = Ini_Read_Int('signal_rand_seed',1000)
    rp%caliberror_rand_seed = Ini_Read_Int('caliberror_rand_seed',10000)


    rp%modulo_npix = Ini_Read_Int('modulo_npix',1)

    rp%sim_noise = Ini_Read_Logical('sim_noise',.false.)
    if (rp%sim_noise) then
       allocate(rp%sim_noise_factor(rp%nchannel))
       do ff = 1, rp%nchannel
          !A linear factor by which to multiply the noise simulation at each channel
          rp%sim_noise_factor(ff) = Ini_Read_Real_array('sim_noise_factor',ff,1.)
       end do
    endif

    allocate(rp%weight_factor(rp%nchannel))
    do ff = 1, rp%nchannel
       !A linear factor by which to multiply the weight (inverse variance) at each channel
       rp%weight_factor(ff) = Ini_Read_Real_array('weight_factor',ff,1.)
    end do


    rp%sim_signal         = Ini_Read_Logical('sim_signal',.false.)
    rp%sim_caliberror     = Ini_Read_Logical('sim_caliberror',.false.)
    rp%sim_suppress_bmode = Ini_Read_Logical('sim_suppress_bmode',.false.)
    rp%sim_component      = Ini_Read_Logical('sim_component',.false.)

    if(rp%sim_signal) then
       rp%cls_file = Ini_Read_String('cls_file')
       rp%sim_lmax = Ini_Read_Int('sim_lmax',2048)
    endif

    if(rp%sim_caliberror) then
       allocate(rp%sim_calib(rp%nchannel))
       allocate(rp%sim_calib_rms(rp%nchannel))
       allocate(rp%sim_offset(rp%nchannel))
       do ff = 1, rp%nchannel
          rp%sim_calib(ff)    = Ini_Read_Real_array('sim_calib',ff,1.)
          rp%sim_calib_rms(ff)    = Ini_Read_Real_array('sim_calib_rms',ff,0.)
          rp%sim_offset(ff)   = Ini_Read_Real_array('sim_offset',ff,0.)
       end do
    end if

    !----------------------------------------------------
    !Parameters relating to the simulations of components
    !----------------------------------------------------
    if(rp%sim_component) then       
       rp%sim_ncomponent = Ini_Read_Int('sim_ncomponent',1)
       allocate(rp%sim_component_list(rp%sim_ncomponent))
       allocate(rp%sim_componentfile(rp%sim_ncomponent))
       allocate(rp%sim_component_frequency(rp%sim_ncomponent))
       allocate(rp%sim_componentfile_factor(rp%sim_ncomponent))
       max_params = 1
       do cc = 1, rp%sim_ncomponent
          rp%sim_component_list(cc)      = Ini_Read_String_array('sim_component',cc)
          rp%sim_componentfile(cc)       = Ini_Read_String_array('sim_componentfile',cc)
          rp%sim_component_frequency(cc) = Ini_Read_Real_array('sim_component_frequency',&
               cc,real(rp%Inst%channel_nu0_ghz(1)))
          rp%sim_componentfile_factor(cc) = Ini_Read_Real_array('sim_componentfile_factor',cc,1.0)
          component = rp%sim_component_list(cc)
!          pp = get_totalforegroundparameters(rp%sim_component_list(cc),1)
          pp = get_totalforegroundparameters(component,1)
          if (pp .gt. max_params) max_params = pp
       end do
       
       !----------------------------------------
       !Allocate memory for component parameters
       !----------------------------------------
       allocate(rp%sim_component_param(rp%sim_ncomponent,max_params))

       !----------------------------------
       !Get simulated component parameters
       !----------------------------------
       do cc = 1, rp%sim_ncomponent
         component = rp%sim_component_list(cc)
!          do pp=1,get_totalforegroundparameters(rp%sim_component_list(cc),1)
          do pp=1,get_totalforegroundparameters(component,1)
             rp%sim_component_param(cc,pp)=& 
                  Ini_Read_Real_array('sim_component'//trim(inttostr(cc))//'_param',pp,0.)
          end do
       end do

       rp%sim_component_perturbations = Ini_Read_Logical('sim_component_perturbations',.false.)

       !----------------------------------------------
       !If yes allocate memory for power spectra names
       !----------------------------------------------
       if(rp%sim_component_perturbations) then 
          rp%component_rand_seed = Ini_Read_Int('component_rand_seed',100000)
          do cc = 1, rp%sim_ncomponent
!            do pp = 1, get_totalforegroundparameters(rp%sim_component_list(cc),1)
            component = rp%sim_component_list(cc)
            do pp = 1, get_totalforegroundparameters(component,1)
             ! Read in scalar power spectrum names for each parameter
               rp%sim_component_cls_file(cc,pp) = &
                    Ini_Read_String_array('sim_component'//trim(inttostr(cc))//'_cls_file_param',pp)
               rp%sim_component_paramvariationsmap_file(cc,pp) = &
                    Ini_Read_String_array('sim_component'//trim(inttostr(cc))//&
                    '_paramvariationsmap_file_param',pp)
            end do
          end do
       endif
    endif

    rp%do_compsep_on_simulated_component = Ini_Read_int('do_compsep_on_simulated_component',-1)
    rp%output_whitened_residuals         = Ini_Read_logical('output_whitened_residuals',.true.)
    rp%output_relative_residuals         = Ini_Read_logical('output_relative_residuals',.true.)
    rp%pixelsofinterestfile              = Ini_Read_String('pixelsofinterestfile')
    rp%convert_units_cmb_to_rj           = Ini_Read_Logical('convert_units_cmb_to_rj',.true.)

    call Ini_Close
    if (rp%feedback .gt. 0) then
       print *, '---------------------------'
       print *, 'Finished reading parameters'
       print *, '---------------------------'
    endif
 end subroutine get_runparams


function stokes_healpix_index(stokes,ss) result(healpix_index)
  character(3), intent(in) :: stokes
  integer, intent(in) :: ss
  integer healpix_index

  select case(trim(stokes))
  case('IQU')
     healpix_index = ss
  case('QU')
     healpix_index = ss + 1
  case('I')
     healpix_index = 1
  case('Q')
     healpix_index = 2
  case('U')
     healpix_index = 3
  end select
 
end function

function stokes2nstokes(stokes) result(nstokes)
  character(3), intent(in) :: stokes
  integer nstokes

  select case(trim(stokes))
  case('IQU')
     nstokes = 3
  case('QU')
     nstokes = 2
  case('I')
     nstokes = 1
  case('Q')
     nstokes = 1
  case('U')
     nstokes = 1
  end select

end function stokes2nstokes


function directory_exists(directory)

  character(*), intent(in) :: directory
  logical directory_exists
  integer unit
  CHARACTER :: delimiter

!  CALL get_environment_variable('DELIMITER',delimiter)
!  print *,'directory exist check for : ',trim(directory)//'/.'
  INQUIRE(DIRECTORY=trim(directory),EXIST=directory_exists)

end function directory_exists


end module RunParameters
