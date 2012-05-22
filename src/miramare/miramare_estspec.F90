program EstSpec
        use IniFile
        use MonteCarlo
        use ParamDef
        use settings
	use Likelihoodcmbforeground
        use CalcLike
        use EstCovmatModule
        use ConjGradModule
        use CMBData
        use InstrumentData
        use RunParameters
        implicit none
        
        character(LEN=Ini_max_string_len) InputFile, LogFile
        logical bad
        integer i, numtoget, action,rr,cc,index,status
        character(LEN=Ini_max_string_len) baseroot,fname
        character(LEN=512) numstr
        Type(ParamSet) Params, EstParams
        real delta_loglike

        character(len=*), parameter :: CODE = "MIRAMARE_ESTSPEC"

#ifdef MPI
        integer ierror
        call mpi_init(ierror)
        if (ierror/=MPI_SUCCESS) stop 'MPI fail: rank'
#endif


#ifndef MPI 
        InputFile = GetParam(1)
        if (InputFile == '') call DoStop('No parameter input file')
#endif
        numstr = GetParam(2)

#ifdef MPI 
        if (instance /= 0) call DoStop('With MPI should not have second parameter')
        call mpi_comm_rank(mpi_comm_world,MPIrank,ierror)
        instance = MPIrank +1 !start at 1 for chains                                             
        write (numstr,*) instance
        if (ierror/=MPI_SUCCESS) call DoStop('MPI fail')
        call mpi_comm_size(mpi_comm_world,MPIchains,ierror)
        if (instance == 1) then
           print *, 'Number of MPI processes:',mpichains
           InputFile = GetParam(1)
        end if
        CALL MPI_Bcast(InputFile, LEN(InputFile), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierror) 
#endif


        !-----------------------------------------------------------
        !Read in component separation parameters from parameter file     
        !-----------------------------------------------------------
        call get_runparams(Inputfile,rp)

        !---------------------------------------------
        !Read in other parameters (eg MCMC parameters)
        !---------------------------------------------
        call Ini_Open(InputFile, 1, bad, .false.)
        if (bad) call DoStop('Error opening parameter file')
        Ini_fail_on_not_found = .false.

        propose_scale = Ini_Read_Real('propose_scale',2.4)


#ifdef MPI 
        MPI_StartTime    = MPI_WTime()
        MPI_R_Stop       = Ini_Read_Real('MPI_Converge_Stop',MPI_R_Stop)
        MPI_LearnPropose = Ini_Read_Logical('MPI_LearnPropose',.true.)  
        if (MPI_LearnPropose) then
           MPI_R_StopProposeUpdate    = Ini_Read_Real('MPI_R_StopProposeUpdate',0.)
           MPI_Max_R_ProposeUpdate    = Ini_Read_Real('MPI_Max_R_ProposeUpdate',2.)
           MPI_Max_R_ProposeUpdateNew = Ini_Read_Real('MPI_Max_R_ProposeUpdateNew',30.)
        end if
        MPI_Check_Limit_Converge = Ini_Read_Logical('MPI_Check_Limit_Converge',.false.)
        MPI_StartSliceSampling   = Ini_Read_Logical('MPI_StartSliceSampling',.false.)
        if (MPI_Check_Limit_Converge) then
           MPI_Limit_Converge     = Ini_Read_Real('MPI_Limit_Converge',0.025)
           MPI_Limit_Converge_Err = Ini_Read_Real('MPI_Limit_Converge_Err',0.3)
           MPI_Limit_Param        = Ini_Read_Int('MPI_Limit_Param',0)
        end if      
#endif

        Ini_fail_on_not_found = .true.        
        baseroot              = Ini_Read_String('file_root')        
        rootname              = trim(rp%out_dir)//'/'//trim(baseroot)
        FileChangeIniAll      = trim(rootname)//'.read'

        if (instance /= 0 .and. mpichains .gt. 1) then
           rootname = trim(rootname)//'_'//trim(adjustl(numstr))
        end if
        

        !------------------------------------------------------
        !Random seed for MCMC starting point in parameter space
        !------------------------------------------------------
#ifdef MPI
        call InitRandom(rp%mcmc_rand_seed+mpirank)
#else
        call InitRandom(rp%mcmc_rand_seed)
#endif

        new_chains   = .true. 
        action       = Ini_Read_Int('action',0)
        chain_format = Ini_Read_Int('chain_format',1)

        FeedBack      = Ini_Read_Int('feedback',0)
        FileChangeIni = trim(rootname)//'.read'
        if (action /= 1) then

           LogFile = trim(rootname)//'.log'
        
           if (LogFile /= '') then
              logfile_unit = 49
              call CreateOpenTxtFile(LogFile,logfile_unit,.not. new_chains)
           else
              logfile_unit = 0
           end if
        
           outfile_unit = 48
           fname        = trim(rootname)//'.txt'
           if (new_chains) call CreateTxtFile(fname,outfile_unit)

           
           indep_sample = Ini_Read_Int('indep_sample')
           if (indep_sample /=0) then
              indepfile_unit = 47
              fname          = trim(rootname)//'.data' 
              call CreateOpenFile(fname,indepfile_unit,'unformatted',.not. new_chains)
           end if

           Ini_fail_on_not_found = .false.
           burn_in         = Ini_Read_Int('burn_in',0)     
           sampling_method = Ini_Read_Int('sampling_method',sampling_metropolis)
           if (sampling_method > 6 .or. sampling_method<1) call DoStop('Unknown sampling method')
           if (sampling_method==4) directional_grid_steps = Ini_Read_Int('directional_grid_steps',20)
        end if

        Temperature           = Ini_Read_Real('temperature',1.)
        num_threads           = Ini_Read_Int('num_threads',0)
        !$ if (num_threads /=0) call OMP_SET_NUM_THREADS(num_threads)
        delta_loglike         = Ini_Read_Real('delta_loglike',2.)
        Ini_fail_on_not_found = .true.
        numtoget              = Ini_Read_Int('samples')


        !--------------------------------------------------------
        !Reads in the mask, or makes a dummy mask containing ones
        !--------------------------------------------------------
        if (rp%maskfile .eq. '') then
           call MakeDummyMask(Hmask,rp%nside_in,getbadpix=.true.)           
           rp%ordering_in = ord_ring !Assume ring ordering
        else
           call ReadMask(Hmask,trim(rp%maskfile),getbadpix=.true.)
           rp%nside_in    = Hmask%nside
           rp%ordering_in = Hmask%ordering 
        endif


        !-------------------------
        !Reads in the regional map
        !-------------------------
        if(rp%regionfile .eq. '') then  !Assume single region case
           call ReadRegionalMask(Region,trim(rp%maskfile),rp%nregion,rp%npixtot,Hmask)
        else
           call ReadRegionalMask(Region,trim(rp%regionfile),rp%nregion,rp%npixtot,Hmask)
        endif
        print *, ' Number of regions    = '//trim(inttostr(rp%nregion))
        print *, ' Number of pixels     = '//trim(inttostr(rp%npixtot))

        !-----------------------------------
        !Copy component model to each region
        !-----------------------------------
        do rr = 1, rp%nregion
           region(rr)%ncomponent       = rp%ncomponent
           allocate(region(rr)%component_list(rp%ncomponent))
           region(rr)%component_list   = rp%component_list
           region(rr)%nparam_component = rp%nparam_component
        enddo

        num_params_used = sum(region(:)%nparam_component)
	if (rp%calibration_marginalisation) num_params_used = num_params_used + rp%nchannel
	print *, ' Number of parameters = '//trim(inttostr(num_params_used))
	num_params_model = num_params_used
	if (chain_format .eq. 2) write(outfile_unit,'(I5)') num_params_model + 2

	call Allocate_ParamSet(Params,num_params_model)
	call Allocate_ParamSet(LastParams,num_params_model)
	call Allocate_ParamSet(EstParams,num_params_model)
	call Allocate_ParamScale(Scales,num_params_model)

        !--------------------------------------
        !Copies parameter ranges to all regions
        !--------------------------------------
	if(rp%calibration_marginalisation) then
           call get_parameter_ranges_and_priors(EstParams,Scales,rp%nparam_component+rp%nchannel)
        else
           call get_parameter_ranges_and_priors(EstParams,Scales,rp%nparam_component)
        endif
        index = rp%nparam_component
        do rr = 2, rp%nregion
           do cc = 1, rp%nparam_component
              index               = index + 1
              Scales%PMin(index)  = Scales%PMin(cc)
              Scales%PMax(index)  = Scales%PMax(cc)
              Scales%PWidth(index)= Scales%PWidth(cc)
              Scales%Pmean(index) = Scales%Pmean(cc)
              Scales%Pvar(index)  = Scales%Pvar(cc)
              EstParams%P(index)  = EstParams%P(cc)
           end do
        end do

        !--------------------------------------------
        !Set up MCMC parameter indices for likelihood
        !--------------------------------------------
        do rr = 1,rp%nregion
           region(rr)%lastparam_index  = sum(region(1:rr)%nparam_component)
           region(rr)%firstparam_index = region(rr)%lastparam_index - region(rr)%nparam_component + 1
        end do

        call Ini_Close

        if(any(Scales%pvar > 0.)) then
          use_prior = .true.
          if (Feedback>0) print *,'Using prior constraints on parameters'
        endif

        !------------------------------------------------------
        !Reads in the data and weights, or performs simulations
        !------------------------------------------------------
        call Likelihood_init


        !-----------------------------------------------------------------------
        !'Preweight' the data for computational speed up of likelihood functions
        !-----------------------------------------------------------------------
        if(rp%likelihood_label .eq. 4 .or. rp%likelihood_label .eq. 5 &
             .or. rp%likelihood_label .eq. 7 .or. rp%likelihood_label .eq. 8 &
             .or. rp%likelihood_label .eq. 11) then
           call Multiply_Ninv_data
        end if
        
        !--------------------------------------------------------------
        !Estimates ln(like) offset as given by first call to likelihood
        !--------------------------------------------------------------
        rp%lnlike_offset = GuessLnLikeOffset(Estparams)
#ifdef MPI        
        CALL MPI_Bcast(rp%lnlike_offset, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierror)
#endif          
        if (MPIRank == 0) then
           if (rp%Feedback > 0) write(*,*) 'Setting LnLike_offset = ',rp%lnlike_offset                
        endif
 
        !---------------------------------------------------------
        !Works out which foreground MCMC parameters are begin used
        !and which parameters are fixed.
        !---------------------------------------------------------
        call Initialize_conjgrad_wrapper(EstParams,Scales)
	print *,'Total number of parameters in MCMC = ',num_params_used


        allocate(propose_matrix(num_params_used, num_params_used))
        propose_matrix(:,:) = 0.0
        if (MPIRank == 0) then

           !-----------------------------------------------
           ! Find likelihood maximum using a descent method
           !-----------------------------------------------
           if (rp%Feedback>0) write (*,*) '------------------------'
           if (rp%Feedback>0) write (*,*) 'Starting best-fit finder'
           if (rp%Feedback>0) write (*,*) '------------------------'
           call conjgrad_wrapper(EstParams,delta_loglike,status) 
           Params = EstParams     !Use best fit params
           if (rp%Feedback>0) write (*,*) '---------------------------'
           if (rp%Feedback>0) write (*,*) 'Best fit parameters values:'
           if (rp%Feedback>0) write (*,*) '---------------------------'

           call CreateTxtFile(trim(rootname)//'.minimum',tmp_file_unit)
           write (tmp_file_unit,*) 'Best fit -log(Like) found: ', Bestfit_loglike
           write (tmp_file_unit,*) '' 
           do i = 1, num_params_used
              if (rp%Feedback>0) write (*,*) params_used(i), ' : ', Params%P(params_used(i))
              write (tmp_file_unit,*)        params_used(i), Params%P(params_used(i))
           end do
           close(tmp_file_unit)

           if (action == 2) then
              if (rp%Feedback>0) then
                 write(*,*) 'Have estimated the minimum, now exiting since action=2'
                 write(*,*) 'Wrote the minimum to file ',trim(rootname)//'.minimum'
              end if
              call DoStop
           end if

           !---------------------------------------------------------------
           ! Estimate covariance of likelihood for the MCMC proposal matrix
           !---------------------------------------------------------------
           if (rp%Feedback>0) write (*,*) '------------------------------------------'
           if (rp%Feedback>0) write (*,*) 'Now estimating propose matrix from Hessian'
           if (rp%Feedback>0) write (*,*) '------------------------------------------'
           ! By default the grid used to estimate the covariance matrix has spacings
           ! such that deltaloglike ~ 4 for each parameter.               
           propose_matrix     = EstCovmat(Params,4.,status)
           has_propose_matrix = status > 0
           if (rp%Feedback>0) write (*,*) ' Estimated covariance matrix:'
           call WriteSqMatrix(trim(rootname) //'.local_invhessian',propose_matrix, num_params_used)
           if (rp%Feedback>0) write(*,*) ' Wrote the local inv Hessian to file ',&
                trim(rootname)//'.local_invhessian'
        endif


#ifdef MPI        
        CALL MPI_Bcast(params%p,num_hard, MPI_REAL, 0, MPI_COMM_WORLD, ierror)
        CALL MPI_Bcast(propose_matrix, size(propose_matrix), MPI_REAL, 0, MPI_COMM_WORLD, ierror)
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
#endif

        call SetProposeMatrix

#ifdef RUNIDLE
        call SetIdle
#endif 
         
        if (action == 0) then
           if (rp%Feedback > 0) write (*,*) '--------------------'          
           if (rp%Feedback > 0) write (*,*) 'Starting Monte-Carlo'          
           if (rp%Feedback > 0) write (*,*) '--------------------'          
           if (rp%Feedback > 0) write (*,*) 'Writing chain to ',trim(fname)
           rp%modulo_npix=1

           !--------------------------------------------------------------
           !Estimates ln(like) offset as given by first call to likelihood
           !--------------------------------------------------------------
           rp%lnlike_offset = GuessLnLikeOffset(Params)
#ifdef MPI
           CALL MPI_Bcast(rp%lnlike_offset, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierror)
#endif
           if (MPIRank == 0) then
              if (rp%Feedback > 0) write(*,*) ' Setting LnLike_offset = ',rp%lnlike_offset
           endif

           !---------
           ! Run MCMC
           !---------     
           call MCMCSample(Params, numtoget)           
           if (rp%Feedback > 0) write (*,*) '--------------------'                        
           if (rp%Feedback > 0) write (*,*) 'Finished Monte-Carlo'           
           if (rp%Feedback > 0) write (*,*) '--------------------'                        
           if (logfile_unit /=0) close(logfile_unit)
           if (indepfile_unit /=0) close(indepfile_unit)           
           close(outfile_unit)
        else
           call DoStop('undefined action')
        end if
        
        call DoStop
         
end program EstSpec
