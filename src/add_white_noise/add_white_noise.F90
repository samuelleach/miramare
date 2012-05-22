program add_white_noise

      use IniFile
      use mpi_stop
      use HealpixObj
      use HealpixVis
      use Random
      
      character(LEN=Ini_max_string_len) InputFile
      logical bad

      type(HealpixMap)   :: map

      real undef      
      integer pp,mm

      !Parameters
      integer noise_seed
      real noise_rms,qu_rms
      logical verbose
      character(LEN=256) :: infile
      character(LEN=256) :: outfile
      character(LEN=256) :: inweightfile
      character(LEN=256) :: outweightfile
      character(LEN=256) :: outrmsfile

      character(len=*), parameter :: CODE = "ADD_WHITE_NOISE"
      real dummy
      
      undef= -1.6375e30
      !----------------------
      !Read in parameter file
      !----------------------
      InputFile = GetParam(1)
      call Ini_Open(InputFile, 1, bad, .false.)
      if (bad) then
         print *,CODE,': Add white noise to T of a T or TQU map. '
         print *,'---------------------'
         print *,'Parameter file format'
         print *,'---------------------'
         print *,'verbose    = T or F (default = F).'
         print *,'infile     = Healpix map filename.'
         print *,'outfile    = Output Healpix map filename.'
         print *,'noise_seed = Noise random seed number.'
         print *,'noise_rms  = Noise RMS value (default = 0.).'
         print *,'qu_rms     = Noise RMS value for Q and U (default = 0.).'
         print *,'inweightfile  = (Optional) Input weight (inverse variance) filename.'
         print *,'outweightfile = (Optional) Output corrected weight filename.'
         print *,'outrmsfile    = (Optional) Output corrected 1 or 3 column RMS (weight^-1) filename.'
         print *,'---------------------'
         call DoStop(trim(concat(CODE,': Error opening parameter file')))
      endif
      Ini_fail_on_not_found = .false.
      if (InputFile == '') call DoStop('No parameter input file')      

      verbose       = Ini_Read_logical('verbose',.false.)
      noise_seed    = Ini_Read_Int('noise_seed',-1)
      noise_rms     = Ini_Read_Real('noise_rms',0.)
      qu_rms        = Ini_Read_Real('qu_rms',0.)
      infile        = Ini_Read_String('infile')
      outfile       = Ini_Read_String('outfile')
      inweightfile  = Ini_Read_String('inweightfile')
      outweightfile = Ini_Read_String('outweightfile')
      outrmsfile    = Ini_Read_String('outrmsfile')

      if (noise_seed .eq. -1) then
        print *,trim(concat(trim(CODE),': ERROR: Must set noise_seed in parameter file.'))
      endif
      
      if(verbose) print *,trim(CODE),': running'


      if(verbose) print *,trim(CODE),': reading ',trim(infile)
      call HealpixMap_read(map,trim(infile))

      
      call InitRandom(noise_seed)
      dummy = gaussian1()
      do pp = 0, map%npix-1
         map%tqu(pp,1) = map%tqu(pp,1) + Gaussian1()*noise_rms
      end do

      if (qu_rms .gt. 0.) then
         do pp = 0, map%npix-1
            map%tqu(pp,2) = map%tqu(pp,2) + Gaussian1()*qu_rms
            map%tqu(pp,3) = map%tqu(pp,3) + Gaussian1()*qu_rms         
         enddo
      end if

      if(verbose) print *,trim(CODE),': writing ',trim(outfile)
      call HealpixMap_write(map,trim(outfile),units=' ')
     

      !-----------------------------------
      ! Output modified output weight file
      !-----------------------------------
      if(outweightfile .ne. '') then

         if (inweightfile .ne. '') then 
            if(verbose) print *,trim(CODE),': reading ',trim(inweightfile)
            call HealpixMap_read(map,trim(inweightfile))

            if (map%nmaps .gt. 3) then
               print *,trim(CODE),' : not yet implemented handing of 6 column weight matrices.'
               stop
            end if

            !Convert weight to variance
            do pp = 0, map%npix-1
               do mm = 1, map%nmaps
                  map%tqu(pp,mm) = 1./ map%tqu(pp,mm)
               end do
            end do

         else
            do pp = 0, map%npix-1
               do mm = 1, map%nmaps
                  map%tqu(pp,mm) = 0. !Zero variance map
               end do
            end do
         endif


         !Add variance
         do pp = 0, map%npix-1
!            map%tqu(pp,1) = ( 1./map%tqu(pp,1) + noise_rms**2)**(-1)
            map%tqu(pp,1) = ( map%tqu(pp,1) + noise_rms**2)**(-1)
         end do
         if (qu_rms .gt. 0.) then
            do pp = 0, map%npix-1               
               map%tqu(pp,2) = ( map%tqu(pp,2) + qu_rms**2)**(-1)
               map%tqu(pp,3) = ( map%tqu(pp,3) + qu_rms**2)**(-1)
            enddo
         endif

         if(outweightfile .ne. '') then
            if(verbose) print *,trim(CODE),': writing ',trim(outweightfile)
            call HealpixMap_write(map,trim(outweightfile),units=' ')
         end if

         if(outrmsfile .ne. '') then
            do pp = 0, map%npix-1
               map%tqu(pp,:) = sqrt( 1./map%tqu(pp,:))
            end do            
            if(verbose) print *,trim(CODE),': writing ',trim(outrmsfile)
            call HealpixMap_write(map,trim(outrmsfile),units=' ')
         end if
      end if
      


    end program add_white_noise
