program newmodule

      use IniFile
      use mpi_stop
      use HealpixObj
      use HealpixVis
      
      character(LEN=Ini_max_string_len) InputFile
      logical bad

      Type(HealpixInfo)  :: HH
      type(HealpixMap)   :: map
      integer :: mpi_division_method = division_equalrows

      real undef      
      integer nside

      !Parameters
      integer lmax
      logical verbose
      real target_fwhm_arcmin
      logical want_pol
      character(LEN=100) :: rootname
      character(LEN=100) :: inmap
      character(LEN=100) :: outmap

      character(len=*), parameter :: CODE = "NEWMODULE"
      
!      debugmsgs=1
      
#ifdef MPIPIX
      integer i
      call mpi_init(i)
      mpi_division_method = 3
#endif
      undef= -1.6375e30
      !----------------------
      !Read in parameter file
      !----------------------
      InputFile = GetParam(1)
      call Ini_Open(InputFile, 1, bad, .false.)
      if (bad) then 
         print *,CODE,': description here'
         print *,'---------------------'
         print *,'Parameter file format'
         print *,'---------------------'
         print *,'verbose = T or F (default = F).'
         print *,'---------------------'
         call DoStop(trim(concat(CODE,': Error opening parameter file')))
      endif
      Ini_fail_on_not_found = .false.
      if (InputFile == '') call DoStop('No parameter input file')      

      verbose = Ini_Read_logical('verbose',.false.)
      lmax = Ini_Read_Int('lmax',750)
      target_fwhm_arcmin = Ini_Read_Real('target_fwhm_arcmin',60.)
      want_pol = Ini_Read_Logical('want_pol',.false.)      
      rootname = Ini_Read_String('file_root')
      inmap = Ini_Read_String('inmap')
      outmap = Ini_Read_String('outmap')

      if(verbose) print *,trim(concat(trim(CODE),': running'))

      call HealpixInit(HH,nside, 2*lmax,want_pol,w8dir='',&
           method= mpi_division_method)


      call HealpixMap_read(map,trim(inmap))
      call HealpixMap_write(map,trim(outmap),units='uK_CMB')

      
      if(HH%MpiID .eq. 0) then  !if we are main thread


      endif                     !MpiID=0





#ifdef MPIPIX
      call HealpixFree(HH)
      call mpi_finalize(i)
#endif
  
    end program Newmodule
