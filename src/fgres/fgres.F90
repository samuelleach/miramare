program fgres
  ! Code to implement residual foreground power spectrum expressions
  ! of Tucci et al 04

      use IniFile
      use mpi_stop
      use HealpixObj
      use HealpixVis
      
      character(LEN=Ini_max_string_len) InputFile
      logical bad

      Type(HealpixInfo)  :: HH
      integer :: mpi_division_method = division_equalrows
      integer i,l1,ll,l2min,l2max,offs,m3

      real undef      
      integer nside

      !Parameters
      integer lmax
      logical verbose
      logical want_pol
      character(len=256) C_beta,C_foreground,C_noise,c_residual
      real nu_template,nu_residual 
      real beta_ave

      character(len=*), parameter :: CODE = "FGRES"
      real(dp), dimension(:), allocatable :: threejs   !(lmax*2+1)
      real(dp) tmp,tmp2
      Type(HealpixPower) Presidual,Pbeta,Pforeground,Pnoise

      
#ifdef MPIPIX
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

      verbose       = Ini_Read_logical('verbose',.false.)
      lmax          = Ini_Read_Int('lmax',750)
      c_beta        = Ini_Read_string('c_beta')
      c_foreground  = Ini_Read_string('c_foreground')
      c_noise       = Ini_Read_string('c_noise')
      c_residual    = Ini_Read_string('c_residual')
      nu_template   = Ini_Read_Real('nu_template',410.)
      nu_residual   = Ini_Read_Real('nu_residual',150.)
      beta_ave      = Ini_Read_Real('beta_ave',1.6)

      if(verbose) print *,trim(concat(trim(CODE),': running'))
      
      call HealpixInit(HH,nside, 2*lmax,want_pol,w8dir='',&
           method= mpi_division_method)

      allocate(threejs((lmax)**2+1))
      call HealpixPower_Init(Presidual,lmax,.false.)

!      do ll=0, lmax
!         Pbeta%cl(ll,1)       = c_beta*(2.*float(ll)+1.)/float(ll)**3.
!         Pforeground%cl(ll,1) = c_foreground*(2.*float(ll)+1.)/float(ll)**2.
!      end do

      call HealpixPower_readfromtextfile(Pbeta,c_beta,lmax,rescale=.false.)
      call HealpixPower_readfromtextfile(Pforeground,c_foreground,lmax,rescale=.false.)
      call HealpixPower_readfromtextfile(Pnoise,c_noise,lmax,rescale=.false.)

      do ll=0, lmax
         print *,ll
         tmp = 0.
         do l1=0, ll
         
            if (ll .ge. 2 .and. l1 .ge. 2) then
               call GetThreeJs(threejs,ll,l1,2,-2) 
            endif
            m3     = 0
            l2min  = max(abs(ll-l1),abs(m3))
            offs   = l2min - abs(ll-l1)
            l2max  = min(ll+l1,lmax)
            
            tmp2 = sum(Pbeta%cl(l2min:l2max,1)*threejs(1+offs:1+offs+l2max-l2min)*&
                 threejs(1+offs:1+offs+l2max-l2min))
            tmp  = tmp + Pforeground%cl(l1,1)*tmp2
         enddo
         Presidual%cl(ll,1) = tmp/16./3.14159*log(nu_template/nu_residual)**2
         Presidual%cl(ll,2) = Pnoise%cl(ll,2) * (nu_residual/nu_template)**(-2.*beta_ave) 
      enddo

      call healpixpower_writecl(Presidual,trim(c_residual))
      
      
      if(HH%MpiID .eq. 0) then  !if we are main thread
      endif                     !MpiID=0


#ifdef MPIPIX
      call HealpixFree(HH)
      call mpi_finalize(i)
#endif
  
    end program Fgres
