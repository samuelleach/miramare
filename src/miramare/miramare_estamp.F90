program EstAmp

      use mpi_stop
      use ParamDef
      use AMLutils
      use IniFile
      use settings
      use HealpixObj
      use HealpixVis
      use foregroundscalings
      use Likelihoodcmbforeground 
      use CMBData
      use InstrumentData
      use runparameters
      use MatrixUtils
      
      character(LEN=Ini_max_string_len) InputFile
      character(LEN=256)  rootdirname,infile,outfile, paramstr
      integer idummy,cc,ccc,dd,ff,gg,i,ii,ss,map_index,bb,rr,ncomp
      integer npix,stokes_index,noisemat_index,nmap,nsing
      integer firstparam_index,lastparam_index
      logical usemeanfile

      integer*8 pp,thispix,pix_index1,pix_index2
      integer   param_index

      real mean,sddev,minfreq,maxfreq
      real, dimension(:), allocatable :: omegavec,omegavec_sddev
      real*8 mean_fgparam(num_hard),var_fgparam(num_hard)
      Type(ParamSet) Params

      type(HealpixMap)   :: map_M,map_comp,map_fgparams,map_cmbin
      type(HealpixMap)   :: map_write, map_residual, map_resinvN, map_resinvd
      Type(HealpixInfo)  :: HH
      integer :: mpi_division_method = division_equalrows
      logical bad,somecomponentsareunpol
      real(dp) gammasquared
      real undef      

      real(dp), dimension(:,:), allocatable :: pomat,tmpmat,tmpmat_reduced,A,At,AM,AMAt
      real(dp), dimension(:,:), allocatable :: aTInvN,N,invdiagN,pupuT,AMAtinvN,Iden,chi2mat
      real(dp), dimension(:,:), allocatable :: Nhat,invNhat,invuNu,M,VT
      real(dp), dimension(:),   allocatable :: smap,dmap,resmap,resmapinvn,lnlike,lnlike_mod,chi2est,D
      real(dp), dimension(:,:), allocatable :: Omega
      real(dp), dimension(:,:), allocatable :: singVecttmp
      real(dp), dimension(:),   allocatable :: compMap,invdiagNU,pu,U,myvect

      !------------------------------------------
      !Spectral index uncertainty correction data
      !------------------------------------------
      real(dp), dimension(:,:,:), allocatable :: tmpmat_partialbeta
      real(dp), dimension(:,:),   allocatable :: nbetabeta, pmatpartialTinvN
      real(dp), dimension(:),     allocatable :: pmatpartialsmap
      real(dp), dimension(:,:),   allocatable :: delta,Mdelta

      !---------------------------------------------
      !Parameters for frequency scaling output files
      !---------------------------------------------
      real(dp), dimension(:,:), allocatable     :: As,freq,weight
      real(dp), dimension(:,:,:), allocatable     :: pixel_data,pixel_rms,dAdbeta
      real(dp), dimension(:,:,:,:), allocatable   :: pixel_As
      real(dp), dimension(:,:,:,:,:), allocatable :: pixel_dAdbetas
      logical output_frequency_scalings,frequency_scalings_calculated
      integer nfreq,nabscissa,pix
      character(LEN =30) format
      character(len=2)   :: intstring

      integer*8         npixel_of_interest!,pixel_of_interest
      type(HealpixMask) pixofinterest
      integer, dimension(:), allocatable     :: ordered_pixofinterest

!      integer mpiid,mpisize
      character(len=*), parameter :: CODE = "MIRAMARE_ESTAMP"

#ifdef MPIPIX
      call mpi_init(i)
      mpi_division_method = 3
#endif

!      call GetMpiStat(MpiID, MpiSize)

      undef = -1.6375e30
      !----------------------
      !Read in parameter file
      !----------------------
      InputFile = GetParam(1)
      call Ini_Open(InputFile, 1, bad, .false.)
      if (bad) call DoStop('Error opening parameter file')
      Ini_fail_on_not_found = .false.
      if (InputFile == '') call DoStop('No parameter input file')      
      call get_runparams(Inputfile,rp)
 
      !--------------------------------------------------------     
      !Reads in the mask, or makes a dummy mask containing ones
      !--------------------------------------------------------     
      if (rp%maskfile .eq. '') then
         call MakeDummyMask(Hmask,rp%nside_in,getbadpix=.true.)           
      else
         call ReadMask(Hmask,trim(rp%maskfile),getbadpix=.true.)
         rp%nside_in    = Hmask%nside
         rp%ordering_in = Hmask%ordering
      endif

      npix = nside2npix(rp%nside_in)

      !--------------------------------
      ! Read in pixels of interest mask
      !--------------------------------
      npixel_of_interest = 0
      if(rp%pixelsofinterestfile .ne. '') then
         if (rp%Feedback>0) print *,' Reading pixels of interest: ',trim(rp%pixelsofinterestfile)
         call ReadMask(PixofInterest,trim(rp%pixelsofinterestfile))
         npixel_of_interest = PixofInterest%npix
         allocate(ordered_pixofinterest(npixel_of_interest))
      end if


      !-----------------------------------------------------
      !NEED A BETTER WAY OF HANDLING REGIONAL MODEL and mask
      !-----------------------------------------------------
      !-------------------------
      !Reads in the regional map
      !-------------------------
      if(rp%regionfile .eq. '') then !Assume single region case
         call ReadRegionalMask(Region,trim(rp%maskfile),rp%nregion,rp%npixtot,Hmask)
      else
         call ReadRegionalMask(Region,trim(rp%regionfile),rp%nregion,rp%npixtot,Hmask)
      endif
      if (rp%Feedback>0) print *, 'Number of regions = ',rp%nregion
      if (rp%Feedback>0) print *, 'Number of pixels  = ',rp%npixtot

      do rr = 1,rp%nregion
         region(rr)%ncomponent       = rp%ncomponent
         allocate(region(rr)%component_list(rp%ncomponent))
         region(rr)%component_list   = rp%component_list
         region(rr)%nparam_component = rp%nparam_component
      enddo

      !------------------
      !Initialise Healpix
      !------------------
      if (rp%Feedback>0) print *,'Initialising MPI Healpix'
      call HealpixInit(HH,rp%nside_in, 2*rp%nside_in,rp%want_pol,&
           w8dir='', method= mpi_division_method)
      if(HH%MpiID .eq. 0) then  !if we are main thread

         rootdirname = concat(rp%out_dir,rp%rootname)

         !-------------------------------------------------------------
         !Read in .mean file (output from miramare_getdist):
         !Contains estimates of mixing matrix parameters and error bars
         !-------------------------------------------------------------
         if (rp%nparam_component .gt. 0) then
            infile      = trim(rootdirname)//'.mean'
            usemeanfile = fileexists(infile)
            select case(usemeanfile)
            case(.true.)
               if (rp%Feedback>0) print *,'Reading: ',trim(infile)
               open(unit=50,file=infile,form='formatted')
               param_index = 0
               do rr = 1, rp%nregion
                  do ii = 1, region(rr)%nparam_component
                     param_index               = param_index + 1
                     read(50,*)                  idummy, mean,sddev
                     mean_fgparam(param_index) = mean
                     var_fgparam(param_index)  = sddev**2
!                     if (rp%Feedback>0) write(*, '(A,I,A,f9.4,A,f9.4)') 'Parameter ',param_index,&
!                          trim(' = '),mean_fgparam(param_index),trim(' +/- '),sqrt(var_fgparam(param_index))
                     if (rp%Feedback>0) write(*, '(A,I3,A,f9.4,A,f9.4)') 'Parameter ',param_index,&
                          trim(' = '),mean_fgparam(param_index),trim(' +/- '),sqrt(var_fgparam(param_index))
                  end do
               end do
            case(.false.)
               if (rp%Feedback>0) print *,'No .mean file found: ',trim(infile)
               if (rp%Feedback>0) print *,'Using spectral parameter values from parameter file'
               call Allocate_ParamSet(Params,rp%nparam_component)
               call Allocate_ParamScale(Scales,rp%nparam_component)
               call Ini_Open(InputFile, 1, bad, .false.)
               call get_parameter_ranges_and_priors(Params,Scales,rp%nparam_component)
               call Ini_Close
               param_index = 0
               do ii = 1, region(1)%nparam_component
                  param_index               = param_index + 1
                  mean_fgparam(param_index) = Scales%center(ii)
                  var_fgparam(param_index)  = Scales%Pwidth(ii)**2
                  if (rp%Feedback>0) write(*, '(A,f9.4,A,f9.4)') 'Parameter '//trim(inttostr(param_index))//&
                       trim(' = '),mean_fgparam(param_index),trim(' +/- '),sqrt(var_fgparam(param_index))
               end do               
            end select
         end if

         !-------------------------------------------
         !Read in estimates of calibration parameters
         !-------------------------------------------
         if(rp%calibration_marginalisation) then
            allocate(omegavec(rp%nchannel),omegavec_sddev(rp%nchannel))
            allocate(Omega(rp%nchannel,rp%nchannel))
            do ff = 1, rp%nchannel
               read(50,*) idummy, omegavec(ff),omegavec_sddev(ff)
            end do
            Omega = 0.0
            do ff = 1, rp%nchannel
               Omega(ff,ff) = omegavec(ff)
            end do
         endif
         close(50)
     
         do ff = 1, rp%nchannel
            if(rp%calibration_marginalisation) then
               if (rp%Feedback>0) print *,'omega(',ff,') = ',omegavec(ff),' +/- ',omegavec_sddev(ff)
            endif
            if (rp%Feedback>0) then
               if (rp%calib(ff) .ne. 1.0)  print *,'calib('//trim(inttostr(ff))//') = ',rp%calib(ff)
               if (rp%offset(ff) .ne. 0.0) print *,'offset('//trim(inttostr(ff))//') = ',rp%offset(ff)
            endif
         end do


         !---------------------------------------------------------
         !Write foreground parameters to disk in Healpix Map format
         !---------------------------------------------------------
         if (rp%nparam_component .gt. 0) then
            call HealpixMap_Init(map_fgparams,npix,nmaps=rp%nparam_component*2,&
                 nested=orderingisnest(rp%ordering_in))
            map_fgparams%tqu(:,:) = undef
            param_index = 0
            do rr = 1,rp%nregion
               do ii=1,region(rr)%nparam_component
                  param_index = param_index + 1
                  map_index   = 1 + 2*(ii-1) 
                  do pp = 1,region(rr)%npix
                     thispix                               = region(rr)%pix(pp)
                     map_fgparams%tqu(thispix,map_index)   = mean_fgparam(param_index)
                     map_fgparams%tqu(thispix,map_index+1) = sqrt(var_fgparam(param_index))                   
                  end do
               end do
            enddo
            outfile = trim(rootdirname)//'_fgparams.fits'
            if (rp%Feedback>0) print *,'Writing foreground parameters to ',trim(outfile)
            call HealpixMap_Write(map_fgparams, '!'//trim(outfile),units='')
            call HealpixMap_free(map_fgparams)
         end if
     
         !------------------------
         !Read in data and weights
         !------------------------
         if (rp%Feedback>0) print *,'Reading data and weights'
         call Likelihood_init

         !-------------------
         !Set up Healpix maps
         !-------------------
         call HealpixMap_Init(map_comp,npix,rp%ncomponent*rp%nstokes,&
                 nested=orderingisnest(rp%ordering_in))
         map_comp%tqu(:,:) = undef
         call HealpixMap_Init(map_residual,npix,rp%nchannel*rp%nstokes,&
              nested=orderingisnest(rp%ordering_in))
         map_residual%tqu  = undef
         if(rp%output_whitened_residuals) then
            call HealpixMap_Init(map_resinvn,npix,rp%nchannel*rp%nstokes,&
                 nested=orderingisnest(rp%ordering_in))
            map_resinvn%tqu = undef
         endif
         if(rp%output_relative_residuals) then
            call HealpixMap_Init(map_resinvd,npix,rp%nchannel*rp%nstokes,&
                 nested=orderingisnest(rp%ordering_in))
            map_resinvd%tqu = undef
         endif
         nmap = (rp%ncomponent*(rp%ncomponent+1)/2*rp%nstokes)
         call HealpixMap_Init(map_M,npix,nmaps=nmap,&
              nested=orderingisnest(rp%ordering_in))

         !--------------------------------------------------------
         !Modify the output fileroot if doing component separation
	 !on a simulated component + noise
         !--------------------------------------------------------
         if(rp%do_compsep_on_simulated_component .gt. -1) then
            write (paramstr,*) rp%do_compsep_on_simulated_component
            rootdirname = concat(rootdirname,'_'//trim(adjustl(paramstr)))
         endif

         !----------------------------------------
         !Set up data for frequency scaling output
         !----------------------------------------
         if (npixel_of_interest .ge. 1) then
            nfreq     = 1000
            nabscissa = 1
            allocate(dAdbeta(rp%nparam_component,nfreq,rp%ncomponent))
            allocate(As(nfreq,rp%ncomponent))
            allocate(freq(nabscissa,nfreq))
            allocate(weight(nabscissa,nfreq))
            allocate(pixel_data(npixel_of_interest,rp%nstokes,rp%nchannel))
            allocate(pixel_rms(npixel_of_interest,rp%nstokes,rp%nchannel))
            allocate(pixel_As(npixel_of_interest,rp%nstokes,nfreq,rp%ncomponent))
            allocate(pixel_dAdbetas(npixel_of_interest,rp%nstokes,rp%nparam_component,nfreq,rp%ncomponent))
            weight(:,:) = 1.0
            minfreq     = minval(rp%Inst%channel_nu0_ghz(:))
            maxfreq     = maxval(rp%Inst%channel_nu0_ghz(:))
            do ff = 1,nfreq
               freq(1,ff) = minfreq*0.8 + float(ff-1)/float(nfreq-1)*(maxfreq-minfreq)*1.2
            enddo
            frequency_scalings_calculated = .false.
            output_frequency_scalings     = .true.            
         else
            output_frequency_scalings     = .false.
         endif
            
         !----------------------------------------------
         !Generalised least squares component separation
         !Stompor et al'08 
         !----------------------------------------------         
         allocate(pomat(rp%nchannel,rp%ncomponent))
         allocate(tmpmat(rp%nchannel,rp%ncomponent))
         allocate(tmpmat_reduced(rp%nchannel,rp%npolcomponent))
         allocate(aTInvN(rp%ncomponent,rp%nchannel))
         allocate(dmap(rp%nchannel))
         allocate(resmap(rp%nchannel))
         allocate(resmapinvN(rp%nchannel))
         allocate(lnlike(rp%nchannel))
         allocate(lnlike_mod(rp%nchannel))
         allocate(invdiagN(rp%nchannel,rp%nchannel))
         allocate(N(rp%nchannel,rp%nchannel))
         allocate(invNhat(rp%nchannel,rp%nchannel))
         allocate(Nhat(rp%nchannel,rp%nchannel))
         allocate(VT(rp%nchannel,rp%nchannel))
         allocate(D(rp%nchannel))
         allocate(AMAt(rp%nchannel,rp%nchannel))
         allocate(AMAtinvN(rp%nchannel,rp%nchannel))
         allocate(iden(rp%nchannel,rp%nchannel))
         allocate(chi2mat(rp%nchannel,rp%nchannel))
         allocate(chi2est(rp%nchannel))
         if (rp%offset_marginalisation) then
            allocate(invdiagNU(rp%nchannel))
            allocate(pu(rp%ncomponent))
            allocate(pupuT(rp%ncomponent,rp%ncomponent))
            allocate(U(rp%nchannel))
            allocate(compMap(rp%nchannel))
            allocate(invuNu(rp%nchannel,rp%nchannel))
            allocate(myvect(rp%nchannel))        
            allocate(singvecttmp(rp%ncomponent,rp%nchannel))
         end if
!         if(rp%noise_matrix_correction) then
         allocate(Nbetabeta(rp%nparam_component,rp%nparam_component))
         allocate(tmpmat_partialbeta(rp%nparam_component,rp%nchannel,rp%ncomponent))
         allocate(pmatpartialTinvN(rp%ncomponent,rp%nchannel))
         allocate(pmatpartialsmap(rp%nchannel))
         allocate(delta(rp%nparam_component,rp%ncomponent))
         allocate(Mdelta(rp%nparam_component,rp%ncomponent))
!         endif

         iden(:,:) = 0.d0
         do ff = 1, rp%nchannel
            iden(ff,ff) = 1d0
         end do         

         !-----------------
         !LOOP OVER REGIONS !Fine as long as same component model on all regions
         !-----------------
         lnlike(:)     = 0.d0
         lnlike_mod(:) = 0.d0
         chi2est(:)    = 0.d0
         do rr = 1, rp%nregion
         lastparam_index  = sum(region(1:rr)%nparam_component)
         firstparam_index = lastparam_index - region(rr)%nparam_component + 1

         !--------------------
         !Set up mixing matrix
         !--------------------
         pomat(:,:) = 0.0
         pomat = get_mixingmatrix(region(rr)%component_list,rp%component_reference_frequency_list,&
              region(rr)%ncomponent,mean_fgparam(firstparam_index:lastparam_index),&
              rp%Inst%bandpass_ghz,rp%Inst%bandpass_weight,rp%nchannel,rp%Inst%nabscissa)
         
         !-------------------------------------------------------------------------
         !Set up rescaled mixing matrix using maximum likelihood calibration values
         !-------------------------------------------------------------------------
         if(rp%calibration_marginalisation) then
            tmpmat = matmul(Omega,pomat)
         else
            tmpmat = pomat
         endif

         !----------------------------------------------------------------------------------------
         !Set up mixing matrix and its partial derivative for optional output of frequency scaling
         !----------------------------------------------------------------------------------------
         if (output_frequency_scalings) then
            As = get_mixingmatrix(region(rr)%component_list,rp%component_reference_frequency_list,&
                 region(rr)%ncomponent,mean_fgparam(firstparam_index:lastparam_index),&
                 freq,weight,nfreq,nabscissa)
            if(rp%calibration_marginalisation) As = matmul(Omega,As)
            
            do bb = 1, region(rr)%nparam_component
               dAdbeta(bb,1:nfreq,1:region(rr)%ncomponent) = &
                    get_mixingmatrix_partialderivative(region(rr)%component_list,& 
                    rp%component_reference_frequency_list,& 
                    region(rr)%ncomponent,mean_fgparam(firstparam_index:lastparam_index),& 
                    freq,nfreq,bb) 
            end do
         endif

         !------------------------------------------------------------------------
         !Get first derivative of mixing matrix needed for noise matrix correction
         !------------------------------------------------------------------------         
!         if(rp%noise_matrix_correction) then
         nbetabeta(:,:) = 0.0
         do bb = 1, region(rr)%nparam_component
            nbetabeta(bb,bb) = var_fgparam(bb)
            tmpmat_partialbeta(bb,1:rp%nchannel,1:region(rr)%ncomponent) = &
                 get_mixingmatrix_partialderivative(region(rr)%component_list,& 
                 rp%component_reference_frequency_list,& 
                 region(rr)%ncomponent,mean_fgparam(firstparam_index:lastparam_index),& 
                 rp%Inst%channel_nu0_ghz,rp%nchannel,bb)
         end do
!         end if

         if (rp%Feedback>0) print *,'Performing GLS component separation on region '//trim(inttostr(rr))
         !---------------------------
         !Loop over Stokes parameters
         !---------------------------
         do ss = 1, rp%nstokes
            !----------------------
            !Write A matrix to disk
            !----------------------
            if (rp%nregion .gt. 1) then
               outfile = trim(rp%out_dir)//'/aMat.map'//trim(inttostr(ss))//&
                    'region'//trim(inttostr(rr))//'.'//trim(rp%rootname)//'.txt'
            else
               outfile = trim(rp%out_dir)//'/aMat.map'//trim(inttostr(ss))//&
                    '.'//trim(rp%rootname)//'.txt'
            endif
            if (rp%Feedback>0) print *,' Writing mixing matrix to ',trim(outfile)
            call Matrix_write(outfile,tmpmat)
!            open(unit=50,file=trim(outfile),form='formatted',status='replace')
!            write(50,*)  tmpmat
!            close(50)

            !---------------------------------------------
            !Write partial derivatives of A matrix to disk
            !---------------------------------------------
            do bb = 1, region(rr)%nparam_component
               if (rp%nregion .gt. 1) then
                  outfile = trim(rp%out_dir)//'/daMat.param'//trim(inttostr(bb))//&
                       '.map'//trim(inttostr(ss))//'.region'//trim(inttostr(rr))//'.'//&
                       trim(rp%rootname)//'.txt'
               else
                  outfile = trim(rp%out_dir)//'/daMat.param'//trim(inttostr(bb))//&
                       '.map'//trim(inttostr(ss))//'.'//trim(rp%rootname)//'.txt'
               endif
               if (rp%Feedback>0) print *,' Writing partial derivatives of mixing matrix to ',trim(outfile)
               call Matrix_write(outfile,tmpmat_partialbeta(bb,:,:))
!               open(unit=50,file=trim(outfile),form='formatted',status='replace')
!               write(50,*)  tmpmat_partialbeta(bb,:,:)
!               close(50)
            end do

            !--------------------------------------------------------------------------
            ! Set to zero the Q and U mixing matrix elements for unpolarized components
            !--------------------------------------------------------------------------
            if (rp%stokes_list(ss) .gt. 1 .and. rp%npolcomponent .lt. rp%ncomponent) then
               do cc = 1, rp%npolcomponent
                  tmpmat_reduced(:,cc)= tmpmat(:,rp%polcomponent_list(cc))
               end do
               allocate(A(rp%nchannel,rp%npolcomponent))
               allocate(AM(rp%nchannel,rp%npolcomponent))
               allocate(smap(rp%npolcomponent))
               allocate(M(rp%npolcomponent,rp%npolcomponent))
               A(:,:)                 = tmpmat_reduced(:,:)
               ncomp                  = rp%npolcomponent
               somecomponentsareunpol = .true.
            else
               allocate(A(rp%nchannel,rp%ncomponent))
               allocate(At(rp%ncomponent,rp%nchannel))
               allocate(AM(rp%nchannel,rp%ncomponent))
               allocate(smap(rp%ncomponent))
               allocate(M(rp%ncomponent,rp%ncomponent))
               A(:,:)                 = tmpmat(:,:)
               ncomp                  = rp%ncomponent
               somecomponentsareunpol = .false.
           endif
           
           stokes_index   = ss
           noisemat_index = (ss-1)*region(rr)%ncomponent*(region(rr)%ncomponent+1)/2           
           if(rp%offset_marginalisation) then
              if(rp%nregion .gt. 1) then 
                 print *, ' miramare_estamp2.F90 ERROR: Offset marginalisation only works for single region case.'
                 stop
              endif
              call offsetmarg_computations(tmpmat,rp%nchannel,Hmask%npix, &
                   rp%ncomponent,stokes_index,SingVecttmp,U,nsing,compMap,invuNu)
           end if
        
           !-------------------------------
           !Loop over pixels of each region
           !-------------------------------
            pix_index1 = 1
            pix_index2 = region(rr)%npix
            !Getting ready for reduced memory mode.
!            pix_index1 = region(rr)%pix_index1   
!            pix_index2 = region(rr)%pix_index1 + region(rr)%npix - 1
!           do pp = 1, region(rr)%npix
            pix   = 0
            do pp = pix_index1, pix_index2
               thispix = region(rr)%pix(pp)
              !-------------------
              !Set up noise matrix
              !-------------------
              invdiagN(:,:) = 0.0
              do ff = 1, rp%nchannel     
                 invdiagN(ff,ff) = MapWeightArray(ff)%TQU(thispix,rp%stokes_list(ss))
                 dmap(ff)        = MapDataArray(ff)%TQU(thispix,rp%stokes_list(ss))           
              end do
           
              !--------------------------
              !Set up (A^t N^-1 A) matrix
              !--------------------------
              if(rp%offset_marginalisation) then
                 !---------------------------
                 !With offset marginalisation
                 !---------------------------
                 invdiagNU = matmul(invdiagN,U) 
                 myvect    = matmul(invuNu,invdiagNU) 
                 do ff=1, rp%nchannel     
                    do gg=1, rp%nchannel     
                       invdiagN(ff,gg) = invdiagN(ff,gg) - myvect(ff)*myvect(gg)
                    end do
                 end do
                 aTInvN = matmul(transpose(tmpmat),invdiagN)
                 M      = matmul(aTInvN,tmpmat)
                 pu     = matmul(transpose(tmpmat),U)
                 do ff = 1, rp%ncomponent    
                    do gg = 1, rp%ncomponent     
                       pupuT(ff,gg) = pu(ff)*pu(gg)
                    end do
                 end do
                 M = M + gammasquared*pupuT
              else
                 !--------------------------------------------
                 !Standard case without offset marginalisation
                 !--------------------------------------------
                 At     = transpose(A)
                 aTInvN = matmul(At,invdiagN)
                 M      = matmul(aTInvN,A)
              endif

              !-----------------------------------------
              !Component error covariance matrix Eq.(12)
              !-----------------------------------------           
              call Matrix_Inverse(M)

              !-------------------------
              !GLS component maps Eq.(8)
              !-------------------------           
              smap = matmul(M,matmul(aTInvN,dmap))

              !------------------------------------------------------------------------------------
              !Get frequency scalings of components and their first derivatives at a selected pixel
              !------------------------------------------------------------------------------------
!              if (any(rp%pixel_of_interest(:) .eq. thispix)) then
              if (any(pixofinterest%good(:) .eq. thispix)) then
                 pix                        = pix + 1
                 ordered_pixofinterest(pix) = thispix                 
                 do ff = 1, nfreq
                    do cc = 1,rp%ncomponent
                       pixel_As(pix,ss,ff,cc) = As(ff,cc)*smap(cc)
                       if (abs(pixel_As(pix,ss,ff,cc)) .lt. 1d-98) pixel_As(pix,ss,ff,cc) = 1d-98 
                    end do
                 end do
                 do ff = 1,rp%nchannel                   
                    pixel_data(pix,ss,ff) = dmap(ff) 
                    pixel_rms(pix,ss,ff)  = 1./sqrt(invdiagN(ff,ff))
                 end do
                 do bb = 1,region(rr)%nparam_component
                    do ff = 1,nfreq
                       do cc = 1,rp%ncomponent
                          pixel_dadbetas(pix,ss,bb,ff,cc) = dAdbeta(bb,ff,cc)*smap(cc)
                       end do
                    end do
                 end do
                 frequency_scalings_calculated = .true.
              endif

              !----------------------------------------------
              !Calculate residual (d-As) maps and naive chi^2
              !----------------------------------------------
              resmap(:)  = dmap(:) - matmul(A,smap)
              resmapinvN = matmul(resmap,invdiagN)
              lnlike(:)  = lnlike(:) + resmapinvN(:)*resmap(:)

              !------------------------
              !Calculate modified chi^2
              !------------------------
!              AM            = matmul(A,M)
!              AMAt          = matmul(AM,At)
!              N             = invDiagN
!              call Matrix_inverse(N)
!              Nhat(:,:)     = N(:,:) - AMAt(:,:)
!              print *,nhat
!              call Matrix_SVD(Nhat,rp%nchannel,rp%nchannel,D,VT)
!              print *,nhat
!              print *,D
!              stop
!              invNhat       = Nhat
!              call Matrix_inverse(invNhat)

!              do ff = 1, rp%nchannel                        
!                 resmapinvN(ff)   = resmap(ff)/Nhat(ff,ff)
!!                 lnlike_mod(ff)   = lnlike_mod(ff) + resmapinvN(ff)*resmap(ff)
!              end do
!              lnlike_mod(:)   = lnlike_mod(:) + resmapinvN(:)*resmap(:)

!              AMAtinvN = matmul(AM,aTInvN)              
!              chi2mat  = matmul(matmul(invdiagN,(iden-AMAtinvN)),&
!                   matmul(invdiagN,transpose(iden-AMAtinvN)))
!              do ff = 1, rp%nchannel                        
!                 chi2est(ff) = chi2est(ff) + chi2mat(ff,ff)
!              enddo
              

              !------------------------------------------------------------
              !Calculate corrections to components noise map given by Eq A4
              !------------------------------------------------------------
              if(rp%noise_matrix_correction) then
                 do bb = 1, region(rr)%nparam_component
                    pmatpartialsmap  = matmul(tmpmat_partialbeta(bb,:,:),smap)
                    pmatpartialTinvN = matmul(transpose(tmpmat_partialbeta(bb,:,:)),invdiagN)
                    delta(bb,:)      = matmul(aTinvN,pmatpartialsmap) - matmul(pmatpartialTinvN,resmap)
                 enddo
                 Mdelta = matmul(M,delta)
                 !Problem here for the two component dust model, with syn and cmb
                 M(:,:) = M(:,:) + matmul(matmul(Mdelta,nbetabeta),transpose(Mdelta))
              endif
              !----------------------------------------------------------
              !Transfer components, residuals and weights to Healpix Maps
              !----------------------------------------------------------
              do cc = 1, ncomp
                 if (somecomponentsareunpol) then
                    ccc = rp%polcomponent_list(cc)
                 else
                    ccc = cc
                 endif
                 map_index                       = ss + rp%nstokes*(ccc-1)
                 map_comp%TQU(thispix,map_index) = smap(cc)
              end do
              
              do ff = 1, rp%nchannel
                 map_index                              = ss + rp%nstokes*(ff-1)
                 map_residual%TQU(thispix,map_index)    = resmap(ff)
                 if(rp%output_whitened_residuals) then
                    map_resinvn%TQU(thispix,map_index)  = resmap(ff)*sqrt(invdiagN(ff,ff))
                 endif
                 if(rp%output_relative_residuals) then
                    map_resinvd%TQU(thispix,map_index)  = resmap(ff)/dmap(ff)
                 endif
              end do
           
              map_index=noisemat_index
              do cc = 1, ncomp
                 do dd = cc, ncomp
                    map_index                    = map_index + 1
                    map_M%TQU(thispix,map_index) = M(cc,dd)
                 end do
              end do

           end do       ! Loop over pixels
           deallocate(A,At,smap,M,AM)
        end do          !Loop over Stokes parameters
        end do          !Loop over regions

        !--------------------
        !DEALLOCATE DATA MAPS
        !--------------------
        do cc = 1, rp%nchannel
           call HealpixMap_free(MapDataArray(cc))
           call HealpixMap_free(MapWeightArray(cc))
        enddo
        deallocate(MapDataArray)
        deallocate(MapWeightArray)

        !-----------------
        !Chi^2 information
        !-----------------
        if (rp%Feedback > 0) then
           print *,'npix                       = ',Hmask%npix
           print *,'npix*nstokes*(nchan-ncomp) = ',rp%nstokes*Hmask%npix*(rp%nchannel-rp%ncomponent)
           print *,'chi^2                      = ',sum(lnlike)
!           print *,'Modified chi^2             = ',sum(lnlike_mod)
           print *,'Breakdown of chi^2 '!("Naive" / "Modified"):'
           do ff = 1,rp%nchannel
              print *,' Channel '//trim(inttostr(ff))//', chi^2 =',real(lnlike(ff))!,lnlike_mod(ff)
           enddo
        endif
       
        !----------------------------------------
        !Transfer components to a new Healpix Map
        !----------------------------------------
        call HealpixMap_Init(map_write,npix,rp%npol,nested=orderingisnest(rp%ordering_in))  
        map_write%tqu(:,:) = undef
               
         !-------------------------------------
         !Output and analysis of component maps
         !-------------------------------------
         do cc = 1 , rp%ncomponent
            if (rp%nstokes .gt. 1) then
               do ss = 1, rp%nstokes
                  map_index = ss + rp%nstokes*(cc-1)
                  do pp = 1, Hmask%npix 
                     map_write%TQU(Hmask%good(pp),stokes_healpix_index(rp%stokes,ss)) = &
                          map_comp%TQU(Hmask%good(pp),map_index)
                  enddo
               end do
            else
               do pp = 1, Hmask%npix 
                  map_write%TQU(Hmask%good(pp),1) = map_comp%TQU(Hmask%good(pp),cc)
               end do
            end if
            
            !------------------------
            !Write components to disk
            !------------------------
            call int2string(cc,intstring)
            outfile = trim(rootdirname)//'_comp_'//trim(intstring)//'.fits'            
	    if (Hmask%nbad .gt. 0) then
      	       map_write%tqu(Hmask%bad(1:Hmask%nbad),1:rp%npol) =  undef !0.0
            endif
            if (rp%Feedback>0) print *,' Writing component to ',trim(outfile)        
            call HealpixMap_Write(map_write, '!'//trim(outfile))
        
         end do              !Loop over components


         !---------------------------------------------------
         !Output of CMB residual map (output CMB - input CMB)
         !---------------------------------------------------
         if(rp%sim_signal .eqv. .true.) then 
            do cc = 1 , rp%ncomponent
               if (rp%component_list(cc) .eq. 'cmb') then
                  write (paramstr,*) 1
                  infile = trim(rp%cachedir)//'/'//trim(rp%rootname)// & 
                       '_cmb_'//trim(adjustl(paramstr))//'.fits'
                  call HealpixMap_read(map_cmbin,infile)
                  map_write%TQU(:,:)=0.0 
                  do ss = 1, rp%nstokes
                     map_index = ss + rp%nstokes*(cc - 1)
                     do pp = 1, Hmask%npix                         
                        map_write%TQU(Hmask%good(pp),rp%stokes_list(ss)) = &
                             map_comp%TQU(Hmask%good(pp),map_index) - &
                             map_cmbin%TQU(Hmask%good(pp),rp%stokes_list(ss))
                     enddo
                  end do
                  outfile = trim(rootdirname)//'_'//trim(rp%component_list(cc))//'_residual.fits'
                  if (rp%Feedback>0) print *,' Writing component residual (output - input) to ',trim(outfile)        
                  call HealpixMap_write(map_write,'!'//trim(outfile))
               endif
            enddo
         endif

         !------------------------------
         !Output of residual maps (d-As)
         !------------------------------
         do ff = 1,rp%nchannel    !loop over channels
            call int2string(ff,intstring)
            outfile = trim(rootdirname)//'_residual_'//intstring//'.fits'
            if (rp%nstokes .gt. 1) then
               do ss = 1, rp%nstokes
                  map_index = ss + rp%nstokes*(ff-1)
                  do pp = 1, Hmask%npix                                             
                     map_write%TQU(Hmask%good(pp),stokes_healpix_index(rp%stokes,ss)) = &
                          map_residual%TQU(Hmask%good(pp),map_index)
                  end do
               end do
            else
               do pp = 1, Hmask%npix                         
                  map_write%TQU(Hmask%good(pp),1) = map_residual%TQU(Hmask%good(pp),ff)
               end do
            end if
            if (rp%Feedback>0) print *,' Writing channel residual (d - As) to ',trim(outfile)        
            call HealpixMap_Write(map_write, '!'//trim(outfile))
               
            !-----------------------------------------------------------
            !Output of residual maps in units of the noise
            !-----------------------------------------------------------
            if(rp%output_whitened_residuals) then
               outfile = trim(rootdirname)//'_whitenedresidual_'//intstring//'.fits'
               if (rp%nstokes .gt. 1) then
                  do ss = 1, rp%nstokes
                     map_index = ss + rp%nstokes*(ff-1)
                     do pp = 1, Hmask%npix                                             
                        map_write%TQU(Hmask%good(pp),stokes_healpix_index(rp%stokes,ss)) = &
                             map_resinvn%TQU(Hmask%good(pp),map_index)
                     end do
                  end do
               else
                  do pp = 1, Hmask%npix                                             
                     map_write%TQU(Hmask%good(pp),1) = map_resinvn%TQU(Hmask%good(pp),ff)
                  end do
               end if
               if (rp%Feedback>0) print *,' Writing whitened residuals to ',trim(outfile)        
               call HealpixMap_Write(map_write, '!'//trim(outfile))
            endif

            !--------------------------------------------
            !Output of residual maps in units of the data
            !--------------------------------------------
            if(rp%output_relative_residuals) then
               outfile = trim(rootdirname)//'_relativeresidual_'//intstring//'.fits'
               if (rp%nstokes .gt. 1) then
                  do ss = 1, rp%nstokes
                     map_index = ss + rp%nstokes*(ff-1)
                     do pp = 1, Hmask%npix                                             
                        map_write%TQU(Hmask%good(pp),stokes_healpix_index(rp%stokes,ss)) = &
                             map_resinvd%TQU(Hmask%good(pp),map_index)
                     end do
                  end do
               else
                  do pp = 1, Hmask%npix                                             
                     map_write%TQU(Hmask%good(pp),1) = map_resinvd%TQU(Hmask%good(pp),ff)
                  end do
               end if
               if (rp%Feedback>0) print *,' Writing relative residuals to ',trim(outfile)        
               call HealpixMap_Write(map_write, '!'//trim(outfile))
            endif


         enddo !Loop over channels


         !------------------------------------
         !Output of component noise covariance
         !------------------------------------
         outfile = trim(rootdirname)//'_Nss.fits'
         nmap    = (rp%ncomponent*(rp%ncomponent+1)/2*rp%nstokes)     
  	 if (rp%Feedback>0) print *,'Writing component noise covariance to '//trim(outfile) 
         if (Hmask%nbad .gt. 0) map_M%tqu(Hmask%bad(:),1:nmap)=undef
         call HealpixMap_Write(map_M, '!'//trim(outfile))
      endif                     !MpiID=0


      !-----------------------------------------------------------
      !Output of component frequency scalings for a selected pixel
      !-----------------------------------------------------------
      if(output_frequency_scalings .and. frequency_scalings_calculated) then         
         do pp = 1, pix
            outfile = trim(rp%out_dir)//'/'//trim(rp%rootname)// & 
                 '_scalings_'//trim(inttostr(ordered_pixofinterest(pp)))//'.dat'         
!            if (rp%Feedback>0) print *,'Writing component frequency scalings to ',trim(outfile)
            open(unit=50,file=trim(outfile),form='formatted',status='replace')
            format = trim(numcat('(1E16.7 , ',rp%nstokes*rp%ncomponent))//'E16.7)'
            do ff = 1,nfreq
               select case(rp%nstokes)
               case(1)
                  write(50,fmt=format)  freq(1,ff), pixel_As(pp,1,ff,1:rp%ncomponent)
               case(2)
                  write(50,fmt=format)  freq(1,ff), pixel_As(pp,1,ff,1:rp%ncomponent), pixel_As(pp,2,ff,1:rp%ncomponent)
               case(3)
                  write(50,fmt=format)  freq(1,ff), pixel_As(pp,1,ff,1:rp%ncomponent), pixel_As(pp,2,ff,1:rp%ncomponent),&
                       pixel_As(pp,3,ff,1:rp%ncomponent)
               end select
            enddo
            close(50)
            
            !------------------------------------------------
            !Output of first derivatives of the mixing matrix
            !------------------------------------------------
            do bb = 1,rp%nparam_component
               outfile = trim(rp%out_dir)//'/'//trim(rp%rootname)// & 
                    '_derivatives_'//trim(inttostr(bb))//'_'//trim(inttostr(ordered_pixofinterest(pp)))//'.dat'         
!               if (rp%Feedback>0) print *,'Writing component frequency scalings to ',trim(outfile)
               open(unit=50,file=trim(outfile),form='formatted',status='replace')
               format = trim(numcat('(1E16.7 , ',rp%nstokes*rp%ncomponent))//'E16.7)'
               do ff = 1,nfreq
                  select case(rp%nstokes)
                  case(1)
                     write(50,fmt=format) freq(1,ff),pixel_dAdbetas(pp,1,bb,ff,:)
                  case(2)
                     write(50,fmt=format) freq(1,ff),pixel_dAdbetas(pp,1,bb,ff,1:rp%ncomponent),&
                          pixel_dAdbetas(pp,2,bb,ff,1:rp%ncomponent)
                  case(3)
                     write(50,fmt=format) freq(1,ff),pixel_dAdbetas(pp,1,bb,ff,1:rp%ncomponent),&
                          pixel_dAdbetas(pp,2,bb,ff,1:rp%ncomponent),pixel_dAdbetas(pp,3,ff,bb,1:rp%ncomponent)
                  end select
               enddo
            end do
            close(50)
            
            !------------------------------------
            !Output of data at the selected pixel
            !------------------------------------
            outfile = trim(rp%out_dir)//'/'//trim(rp%rootname)// & 
                 '_data_'//trim(inttostr(ordered_pixofinterest(pp)))//'.dat'         
!            if (rp%Feedback>0) print *,'Writing pixel data to ',trim(outfile)
            open(unit=50,file=trim(outfile),form='formatted',status='replace')
            format = trim(numcat('(1E16.7 , ',2*rp%nstokes))//'E16.7)'
            do ff = 1,rp%nchannel
               select case(rp%nstokes)
               case(1)
                  write(50,fmt=format)  rp%Inst%channel_nu0_ghz(ff), pixel_data(pp,1,ff), pixel_rms(pp,1,ff)
               case(2)
                  write(50,fmt=format)  rp%Inst%channel_nu0_ghz(ff), pixel_data(pp,1,ff), pixel_rms(pp,1,ff),&
                       pixel_data(pp,2,ff), pixel_rms(pp,2,ff)
               case(3)
                  write(50,fmt=format)  rp%Inst%channel_nu0_ghz(ff), pixel_data(pp,1,ff), pixel_rms(pp,1,ff),&
                       pixel_data(pp,2,ff), pixel_rms(pp,2,ff), pixel_data(pp,3,ff), pixel_rms(pp,3,ff)
               end select
            enddo
            close(50)
         end do
      endif


#ifdef MPIPIX
      call HealpixFree(HH)
      call mpi_finalize(i)
#endif
  
end program EstAmp
