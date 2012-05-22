program ud_grade_mask

      use IniFile
      use mpi_stop
      use HealpixObj
      use HealpixVis
      use pix_tools
      
      character(LEN=Ini_max_string_len) InputFile
      logical bad

      type(HealpixMap)   :: maskin
      type(HealpixMap)   :: maskout

      real undef      
      integer npix,pp     
      real tol


      !Parameters
      logical verbose
      character(LEN=100) :: infile
      character(LEN=100) :: outfile
      character(LEN=256) :: ordering_out

      integer nside

      character(len=*), parameter :: CODE = "UD_GRADE_MASK"
      integer ord_in,ord_out
      
      
      undef= -1.6375e30
      !----------------------
      !Read in parameter file
      !----------------------
      InputFile = GetParam(1)
      call Ini_Open(InputFile, 1, bad, .false.)
      if (bad) then
         print *,CODE,': Up or degrading of healpix masks - rejects boundary pixels.'
         print *,'---------------------'
         print *,'Parameter file format'
         print *,'---------------------'
         print *,'verbose      = T or F (default = F).'
         print *,'infile       = Input Healpix mask filename.'
         print *,'outfile      = Output Healpix mask filename.'
         print *,'nside        = Output Healpix mask nside value (default = 512)'
         print *,'ordering_out = Output Healpix mask ordering (nest or ring, default = input ordering)'
         print *,'---------------------'
         call DoStop(trim(concat(CODE,': Error opening parameter file')))
      endif
      Ini_fail_on_not_found = .false.
      if (InputFile == '') call DoStop('No parameter input file')      

      verbose      = Ini_Read_logical('verbose',.false.)
      infile       = Ini_Read_String('infile')
      outfile      = Ini_Read_String('outfile')
      nside        = Ini_Read_Int('nside',512)
      ordering_out = Ini_Read_String('ordering_out')

      if(verbose) print *,trim(CODE),': running'


      select case (infile .eq. '') 
      case(.true.)
         !Make a dummy mask with all 1s.
         if(verbose) print *,trim(CODE),': making dummy mask at nside ',nside
         npix = nside2npix(nside)
         call HealpixMap_Init(maskout, npix, 1)
         maskout%TQU(:,1) = 1.0
      case(.false.)
         if(verbose) print *,trim(CODE),': reading ',trim(infile)
         call HealpixMap_read(maskin,trim(infile))
         !Check is a mask
         tol = 1.d-3 
         do pp = 0, maskin%npix-1
            if(abs(maskin%tqu(pp,1) - 1) .gt. tol) then
               if(abs(maskin%tqu(pp,1)) .gt. tol) then
                  print *,trim(CODE),': ERROR :',trim(infile),'contains values other than 1 or 0'
                  stop
               end if
            end if
            ord_in              = maskin%ordering
         end do
                  
         !-------------
         !Ud grade mask
         !-------------
         if(verbose) print *,trim(CODE),': ud_grading mask from nside ',&
              maskin%nside,' to ',nside
         call HealpixMap_udgrade(maskin, maskout, nside,pessimistic=.true.)

         !-------------
         !Clean up mask
         !-------------
         where( maskout%TQU(:,1) .lt. 1.0)  maskout%TQU(:,1)=0.
      end select

      
      !-----------
      !Reorder map
      !-----------
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
      
      select case(ord_out)
      case(ord_ring)
         if(verbose) print *,trim(CODE)//': Forcing output mask to ring ordering.'
         call HealpixMap_ForceRing(maskout)
      case(ord_nest)
         if(verbose) print *,trim(CODE)//': Forcing output mask to nest ordering.'
         call HealpixMap_ForceNest(maskout)
      end select
      


      if(verbose) print *,trim(CODE),': writing ',trim(outfile)
      call HealpixMap_write(maskout,trim(outfile),units=' ')

      

    end program ud_grade_mask
