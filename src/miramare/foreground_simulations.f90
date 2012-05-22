module foreground_simulations

  use foregroundscalings
  use runparameters
  use HealpixObj

contains

subroutine get_fg_simulation(component_number_list,ncomp,freq_list, &
                            MapSet,nfreq,runpar,weight,writetocache,&
                            pixlist)
    integer,         intent(in) :: ncomp,nfreq
    integer,         intent(in) :: component_number_list(ncomp)
    real*8,          intent(in) :: freq_list(nfreq)
    Type(HealpixMap)            :: MapSet(nfreq)
    type(RunParams)             :: runpar
    real*8,    optional, intent(in) :: weight(nfreq)
    logical,   optional, intent(in) :: writetocache
    integer*8, intent(in), optional :: pixlist(:)

    integer cc,ccc,ff,vv,perturb_seed,nparam,comp
    integer*8 pix,pp,mm
    Type(HealpixInfo):: HH
    Type(HealpixMap) :: Htemplate,smoothmap,Htemplate_out
    Type(HealpixAlm) :: A
    real(dp) tmpmat(nfreq,1)   
    type(HealpixMap), allocatable :: comp_variations(:)
    real(dp), allocatable         :: comp_params(:)
    type(HealpixPower) comp_variations_power
    character(LEN=120) paramstr,paramstr2
    character(LEN=512) filename
    real*8 calib(nfreq)
    real*8 dummy
    character(len=1)   :: intstring1,intstring2
!    logical dummy_logical

    dummy         = freq_list(1) !Avoid compiler complaints
!    dummy_logical = writetocache !Avoid compiler complaints

    if(.not. present(weight)) then
       do ff = 1, nfreq  
          calib(ff) = 1.0
       end do      
    else
       do ff = 1, nfreq  
          calib(ff) = weight(ff)
       end do      
    endif

    call HealpixInit(HH,runpar%nside_in, 2*runpar%sim_lmax,runpar%want_pol, w8dir='')

    do comp = 1, ncomp
       cc = component_number_list(comp)
       if (runpar%Feedback>0) print *,'Reading component template from ',&
            trim(runpar%sim_componentfile(cc))
       call HealpixMap_read(Htemplate,trim(runpar%sim_componentfile(cc)))
!       nparam = get_totalforegroundparameters(runpar%sim_component_list(cc),1)

       if(.not. runpar%sim_component_perturbations) then
          !---------------------------------
          !Simulate fixed scaling components
          !---------------------------------
          tmpmat = get_mixingmatrix(runpar%sim_component_list(cc),& 
               runpar%sim_component_frequency(cc),1,&
               runpar%sim_component_param(cc,1:nparam),&
               runpar%Inst%bandpass_ghz,runpar%Inst%bandpass_weight,nfreq,&
               runpar%Inst%nabscissa)
          
          do ff = 1, nfreq
             if (runpar%Feedback>0) then
               print *,'Scaling component ',cc,' to channel ',&
                ff,' by ',tmpmat(ff,1)*calib(ff)
             endif
             if (.not. present(pixlist)) then 
                do pp = 0, MapSet(ff)%npix-1
                   do mm = 1, MapSet(ff)%nmaps
                      MapSet(ff)%TQU(pp,mm) = MapSet(ff)%TQU(pp,mm) + &
                           Htemplate%TQU(pp,mm)*tmpmat(ff,1)*calib(ff)* &
                           runpar%sim_componentfile_factor(cc)
                   end do
                end do
             else
                do pp=0,MapSet(ff)%npix-1
                   MapSet(ff)%TQU(pp,:) = MapSet(ff)%TQU(pp,:) + &
                        Htemplate%TQU(pixlist(pp+1),:)*tmpmat(ff,1)*calib(ff)*&
                        runpar%sim_componentfile_factor(cc)
                           
                enddo
             endif

             !--------------------------------------------------
             !Write simulated components to disk at each channel
             !--------------------------------------------------
             call int2string(cc,intstring1)
             call int2string(ff,intstring2)
             filename = trim(runpar%cachedir)//'/'//trim(runpar%rootname)//'_simcomponent'//&
                  trim(intstring1)//'_'//trim(intstring2)//'.fits'
             if (runpar%feedback > 0) print *, 'Writing ',trim(filename)
             call healpixmap_assign(htemplate_out,htemplate)
             Htemplate_out%tqu(:,:) = Htemplate%TQU(:,:)*tmpmat(ff,1)*calib(ff)*runpar%sim_componentfile_factor(cc)
             call HealpixMap_write(Htemplate_out,filename)


          end do
       else
          !-----------------------------------------------------
          !Simulate components with spectral index perturbations
          !-----------------------------------------------------

          !--------------------
          !Work out random seed
          !--------------------
          perturb_seed = runpar%component_rand_seed
          do ccc = 1, cc-1
!            perturb_seed = perturb_seed + &
!              get_totalforegroundparameters(runpar%sim_component_list(ccc),1)             
          end do

          allocate(comp_variations(nparam))
          allocate(comp_params(nparam))
          do vv = 1, nparam
             !Initialise maps
             call HealpixMap_Nullify(comp_variations(vv))
             call HealpixMap_Init(comp_variations(vv), nside2npix(runpar%nside_in),&
                  pol = .false.,nested=orderingisnest(runpar%ordering_in))

             ! Get beta0 value
             comp_variations(vv)%tqu = runpar%sim_component_param(cc,vv)
             
             !Deal with Gaussian isotropic spectral index perturbations
             !C_l of variations are for the dimensionless delta(beta)/beta0
             if(runpar%sim_component_cls_file(cc,vv) .ne. '' ) then
                if (runpar%feedback > 0) print *,'Reading spectral index variation deltaBeta/Beta0 Cl spectrum from ',&
                     trim(runpar%sim_component_cls_file(cc,vv))
                call HealpixPower_ReadFromTextFile(comp_variations_power,&
                     runpar%sim_component_cls_file(cc,vv),&
                     runpar%sim_lmax,pol=.false.,dolens = .false.)             
                if(sum(comp_variations_power%cl(:,1)) .gt. 0.0) then
                   print *,'Simulating component variations'
                   perturb_seed = perturb_seed + 1
                   call HealpixAlm_Sim(A, comp_variations_power,perturb_seed,&
                        HasPhi=.false., dopol = .false.)
                   call HealpixAlm2Map(HH,A,smoothmap,nside2npix(runpar%nside_in))
                   do pp = 0,comp_variations(vv)%npix 
                      comp_variations(vv)%tqu(pp,1) = comp_variations(vv)%tqu(pp,1) + &
                           smoothmap%tqu(pp,1)*runpar%sim_component_param(cc,vv)
                   end do
                endif
             endif

             !Deal with map level spectral index perturbation
             !Variations are spectral index variations maps (beta-beta0) 
             if(runpar%sim_component_paramvariationsmap_file(cc,vv) .ne. '' ) then
                if (runpar%feedback > 0) &
                     print *,'Reading spectral spectral index variation map from from ',&
                     trim(runpar%sim_component_paramvariationsmap_file(cc,vv))
                call HealpixMap_read(smoothmap,runpar%sim_component_paramvariationsmap_file(cc,vv))
                comp_variations(vv)%tqu(:,1) = comp_variations(vv)%tqu(:,1) + smoothmap%tqu(:,1)
             endif

             !Write component parameter maps to disk (even constant maps)
             write (paramstr,*) cc
             write (paramstr2,*) vv
             filename=trim(runpar%cachedir)//'/'//trim(runpar%rootname)//'_component'//&
                  trim(adjustl(paramstr))//'_param'//trim(adjustl(paramstr2))//'.fits'
             if (runpar%feedback > 0) print *, 'Writing ',trim(filename)
             call HealpixMap_write(comp_variations(vv),filename)                  

          end do
          
          !---------------------------------------------------------------
          !Perform extrapolation of template using perturbed fg parameters
          !---------------------------------------------------------------
!          do pp = 0, MapSet(ff)%npix-1,30 !Buggy :  MapSet(ff)%npix is too big.
!          do pp = 0, nside2npix(MapSet(ff)%nside)-1 !Buggy - nside=0
          do pp = 0, Htemplate%npix-1
             do vv=1,nparam
                comp_params(vv) = comp_variations(vv)%tqu(pp,1)
             end do
             tmpmat = get_mixingmatrix(runpar%sim_component_list(cc),& 
                  runpar%sim_component_frequency(cc),1,&
                  comp_params(1:nparam),&
                  runpar%Inst%bandpass_ghz,runpar%Inst%bandpass_weight,nfreq,&
                  runpar%Inst%nabscissa)

             do ff = 1, runpar%nchannel
                if (.not. present(pixlist)) then
                   pix = pp
                else
                   pix = pixlist(pp+1)
                endif
                if (pix .lt. MapSet(ff)%npix) then
                   MapSet(ff)%TQU(pix,:) = MapSet(ff)%TQU(pix,:) + &
                        Htemplate%TQU(pix,:)*tmpmat(ff,1)*calib(ff)*&
                        runpar%sim_componentfile_factor(cc)

                endif
             end do
          end do
          deallocate(comp_variations)
          deallocate(comp_params)
       endif
    end do
    call HealpixMap_Free(Htemplate)
    call HealpixFree(HH)

  end subroutine get_fg_simulation
  

  subroutine get_cmb_simulation(freq_list,MapSet,nfreq,runpar,weight,writetocache,pixlist)
    integer, intent(in) :: nfreq
    real*8, intent(in):: freq_list(nfreq)
    Type(HealpixMap) :: MapSet(nfreq)
    type(RunParams) :: runpar
    real*8, optional, intent(in):: weight(nfreq)
    logical, optional, intent(in) :: writetocache    
    integer*8, intent(in), optional :: pixlist(:)

    integer ff   
    Type(HealpixInfo)  :: HH
    type(HealpixPower) Psim 
    Type(HealpixAlm)   :: A,smoothA
    real(dp) fwhm_deg
    Type(HealpixMap) :: smoothmap
    character(LEN=512) filename
    real*8 calib(nfreq)
    logical do_sim
    integer*8 pp,mm
    real(dp) :: tmpmat(nfreq,1)
    character(32), dimension(:), allocatable :: component
    real*8 :: ref_freq(1)
    real*8 :: fg_param(1)
    real*8 dummy

    dummy=freq_list(1) !Avoid compiler warnings
    
    allocate(component(1))
    component(1) = 'cmb'
    ref_freq     = 1.
    fg_param     = 1.

    do_sim   = .true.
    filename = trim(runpar%cachedir)//'/'//trim(runpar%rootname)//'_cmb_1.fits'

    !----------------------------------------------------
    !Read in cached CMB simulation or perform simulations
    !----------------------------------------------------
    if(present(writetocache)) then
       if(writetocache) then
          if(FileExists(filename)) then 
             if (runpar%feedback > 0) print *, 'Reading cached file ',trim(filename)
             call HealpixMap_read(smoothMap,filename)
             do_sim = .false.
          endif
       endif
    endif

    if(do_sim) then
       call HealpixInit(HH,runpar%nside_in, 2*runpar%sim_lmax,runpar%want_pol, w8dir='')
       if (runpar%Feedback>0) print *,' Simulating the CMB with seed = ',runpar%signal_rand_seed
       call HealpixPower_ReadFromTextFile(Psim,runpar%cls_file,runpar%sim_lmax,pol=runpar%want_pol,dolens = .false.)

       call HealpixAlm_Sim(A, Psim, runpar%signal_rand_seed,HasPhi=.false., dopol = runpar%want_pol)    

       smoothA  = A
       fwhm_deg = runpar%fwhm_arcmin(1)/60.
       if (runpar%Feedback>0) print *,' Smoothing CMB alm to ',real(fwhm_deg*60.),trim(' arcmin')
       call HealpixAlm_Smooth(smoothA,fwhm_deg)

       !-------------------------
       !Option to supress B modes
       !-------------------------
       if(runpar%sim_suppress_bmode) then
          if (runpar%Feedback>0) print *,'Setting B mode to zero'
          smoothA%TEB(3,:,:) = 0.
       endif
!       smoothA%TEB(2,:,:) = 0.

       if (runpar%Feedback>0) print *,' Synthesising CMB map'
       call HealpixAlm2Map(HH,smoothA,smoothMap,nside2npix(runpar%nside_in))
       
       if (orderingisnest(runpar%ordering_in)) then
          print *,' Reordering CMB map to nested'
          call HealpixMap_ForceNest(smoothmap)
       endif

       call HealpixFree(HH)

       
       !----------------------------------
       !Cache the smoothed CMB map to disk
       !----------------------------------
       if (present(writetocache)) then
          if(writetocache) then
             if (runpar%feedback > 0) print *, 'Writing ',trim(filename)
             call HealpixMap_write(smoothMap,trim(filename))
          endif
       endif
    endif

    if(.not. present(weight)) then
       do ff = 1, nfreq  
          calib(ff) = 1.0
       end do      
    else
       do ff = 1, nfreq  
          calib(ff) = weight(ff)
       end do      
    endif
    
    !-------------------------------
    !Convert CMB to RJ units and add
    !------------------------------- 
    tmpmat =  get_mixingmatrix(component,ref_freq,1,fg_param,&
         runpar%Inst%bandpass_ghz,runpar%Inst%bandpass_weight,nfreq,&
         runpar%Inst%nabscissa)

    
    do ff = 1, nfreq
       if (.not. present(pixlist)) then
          do pp = 0, MapSet(ff)%npix-1
             do mm = 1, MapSet(ff)%nmaps
                MapSet(ff)%TQU(pp,mm) = MapSet(ff)%TQU(pp,mm) + &
                     smoothMap%TQU(pp,mm)*calib(ff)*tmpmat(ff,1)
             enddo
          enddo
       else
          do pp = 0, MapSet(ff)%npix-1
             do mm = 1, MapSet(ff)%nmaps
                MapSet(ff)%TQU(pp,mm) = MapSet(ff)%TQU(pp,mm) + &
                     smoothMap%TQU(pixlist(pp+1),mm)*calib(ff)*tmpmat(ff,1)
             end do
          enddo
       endif
         
    end do
    call HealpixMap_Free(smoothMap)

    deallocate(component)
    
  end subroutine get_cmb_simulation
  

end module foreground_simulations
