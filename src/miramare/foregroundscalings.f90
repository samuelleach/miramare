module foregroundscalings
  use healpix_types, ONLY: SP,DP,I4B,SPC
  use miramare_utils
  implicit none

  !-----------
  !Global data
  !-----------
  
  integer, parameter              :: nmaxscalingpoints = 750
  real(dp), allocatable, target   :: tabulated_scaling(:,:,:)
  integer, allocatable, target    :: tabulated_scaling_numpoint(:)


  integer, parameter :: nmaxcomponents = 20
  
  real*8 :: h_over_k = 0.0479924d0   ! = 6.626068e-34/1.3806503e-23*1e9              ! K/GHz
  real*8 :: T_E      = 8000d0        ! = Electron temperature for free-free emission ! K


  !From Table 1 of FDS'98
  real*8 :: T1_FDSMODEL1      = 20.0d0
  real*8 :: alpha1_FDSMODEL1  = 1.5d0
  real*8 :: f1_FDSMODEL1      = 1d0
  real*8 :: q1perq2_FDSMODEL1 = 1d0
  
  real*8 :: T1_FDSMODEL2      = 19.2d0
  real*8 :: alpha1_FDSMODEL2  = 1.7d0
  real*8 :: f1_FDSMODEL2      = 1d0
  real*8 :: q1perq2_FDSMODEL2 = 1d0

  real*8 :: T1_FDSMODEL3      = 18.1d0
  real*8 :: alpha1_FDSMODEL3  = 2.0d0
  real*8 :: f1_FDSMODEL3      = 1d0
  real*8 :: q1perq2_FDSMODEL3 = 1d0

  real*8 :: T1_FDSMODEL4      = 17.4d0
  real*8 :: alpha1_FDSMODEL4  = 2.2d0
  real*8 :: f1_FDSMODEL4      = 1d0
  real*8 :: q1perq2_FDSMODEL4 = 1d0

  real*8 :: T1_FDSMODEL5      = 17.0d0
  real*8 :: T2_FDSMODEL5      = 17.0d0
  real*8 :: alpha1_FDSMODEL5  = 1.5d0
  real*8 :: alpha2_FDSMODEL5  = 2.6d0
  real*8 :: f1_FDSMODEL5      = 0.25d0
  real*8 :: q1perq2_FDSMODEL5 = 0.61d0

  real*8 :: T1_FDSMODEL6      = 4.9d0
  real*8 :: T2_FDSMODEL6      = 18.1d0
  real*8 :: alpha1_FDSMODEL6  = 2.0d0
  real*8 :: alpha2_FDSMODEL6  = 2.0d0
  real*8 :: f1_FDSMODEL6      = 0.00261d0
  real*8 :: q1perq2_FDSMODEL6 = 2480d0

  real*8 :: T1_FDSMODEL7      = 9.6d0
  real*8 :: T2_FDSMODEL7      = 16.4d0
  real*8 :: alpha1_FDSMODEL7  = 1.5d0
  real*8 :: alpha2_FDSMODEL7  = 2.6d0
  real*8 :: f1_FDSMODEL7      = 0.0309d0
  real*8 :: q1perq2_FDSMODEL7 = 11.2d0

  real*8 :: T1_FDSMODEL8      = 9.4d0
  real*8 :: T2_FDSMODEL8      = 16.2d0
  real*8 :: alpha1_FDSMODEL8  = 1.67d0
  real*8 :: alpha2_FDSMODEL8  = 2.70d0
  real*8 :: f1_FDSMODEL8      = 0.0363d0
  real*8 :: q1perq2_FDSMODEL8 = 13.0d0

  integer :: CO_nline = 2
  real*8  :: CO_line1 = 115. !GHz
  real*8  :: CO_line2 = 230. !GHz


contains

  function S_greybody(nu_ghz,nu0,beta,Temp) result (S)
    implicit none
    real*8 S
    real*8 nu_ghz, beta,nu0,Temp,nu_dust,f
    
    nu_dust = Temp/h_over_k
    f       = (exp(nu0/nu_dust)-1d0)/(exp(nu_ghz/nu_dust)-1d0)
    S       = f*(nu_ghz/nu0)**(1d0+beta)
    
    return
    
  end function S_greybody

  function S_greybody_dbeta(nu_ghz,nu0,beta,Temp) result (S)
    implicit none
    real*8 S
    real*8 nu_ghz, beta,nu0,Temp,nu_dust,f
    
    nu_dust = Temp/h_over_k
    f       = (exp(nu0/nu_dust)-1d0)/(exp(nu_ghz/nu_dust)-1d0)
    S       = f*(nu_ghz/nu0)**(1d0+beta)*log(nu_ghz/nu0)
    
    return
    
  end function S_greybody_dbeta

  function S_greybody_dtemp(nu_ghz,nu0,beta,Temp) result (S)
    implicit none
    real*8 S
    real*8 nu_ghz, beta,nu0,Temp,nu_dust,x,x0,expx,expx0

    nu_dust = Temp/h_over_k
    x       = nu_ghz/nu_dust
    x0      = nu0/nu_dust
    expx    = exp(x)
    expx0   = exp(x0)

    S = ( x*expx*(expx0 - 1d0) - x0*expx0*(expx-1d0) )/( expx - 1d0 )**2.
    S = S/Temp*(nu_ghz/nu0)**(beta+1d0)
    
    return
    
  end function S_greybody_dtemp
  
  function S_dust_fds(nu_ghz,nu0,alpha1,Temp1,alpha2,Temp2,f1,q1perq2) result (S)
    implicit none
    real*8 S
    real*8 nu_ghz, d,d0,beta,nu0,alpha1,alpha2,Temp1,Temp2,f1,q1perq2

    real*8 nu1_dust,x1,x10,expx1,expx10,nu2_dust,x2,x20,expx2,expx20
    real*8 nu0dust
    real*8 nu(2),SS(2)
    integer ii

    nu(1)    = nu_ghz
    nu(2)    = nu0
    nu0dust  = 3000.

    nu1_dust = Temp1/h_over_k
    nu2_dust = Temp2/h_over_k
    x10      = nu0dust/nu1_dust
    x20      = nu0dust/nu2_dust
    expx10   = exp(x10)
    expx20   = exp(x20)

    do ii = 1, 2
       x1       = nu(ii)/nu1_dust
       x2       = nu(ii)/nu2_dust
       expx1    = exp(x1)
       expx2    = exp(x2)       
       
       d  = q1perq2*f1*(nu(ii)/nu0dust)**(1.0+alpha1)/(expx1-1d0) + &
            (1d0-f1)*(nu(ii)/nu0dust)**(1.0+alpha2)/(expx2-1d0)
       
       d0 = q1perq2*f1/(expx10-1d0) + (1d0-f1)/(expx20-1d0)
       
       beta      = dlog(d/d0)/dlog(nu(ii)/nu0dust)
       SS(ii)    = (nu(ii)/nu0dust)**beta
       
    end do
    
    S = SS(1)/SS(2)

    return
    
  end function S_dust_fds
  
  function S_cmb(nu_ghz) result (S)
    implicit none
    
    real*8 S,x,expx
    real*8 nu_ghz
    
    x    = nu_ghz/56.78d0!  = T k / h / 1e9 = (2.725/6.62607554d-34*1.38065812d-23)/1e9
    expx = exp(x)
    S    = x**2*expx/(expx-1d0)**2
    
    return
  end function S_cmb

  function S_dq(nu_ghz) result (S)
    implicit none
    
    real*8 S,x,expx
    real*8 nu_ghz
    
    x    = nu_ghz/56.78d0!  = T k / h / 1e9 = (2.725/6.62607554d-34*1.38065812d-23)/1e9
    expx = exp(x)
    S    = x**2*expx/(expx-1d0)**2
    
    return
  end function S_dq

  function S_absolute(nu_ghz) result (S)
    implicit none
    
    real*8 S,x,expx
    real*8 nu_ghz
    
    x    = nu_ghz/56.78d0
    expx = exp(x)
    S    = x/(expx-1d0)
    
    return
  end function S_absolute

  
  function powerlawslope_cmb(nu_ghz) result (b)
    implicit none
   
    real*8 b,x
    real*8 nu_ghz

    x = nu_ghz/56.68d0
    b = 2.0d0 + x - 2.0d0*x*exp(x)/(exp(x)-1d0) 

    return
  end function powerlawslope_cmb

  function S_thermalsz(nu_ghz) result (S)
    implicit none
   
    real*8 S,x,expx
    real*8 nu_ghz
    
    x    = nu_ghz/56.78d0
    expx = exp(x)
    S    = x**2*expx/(expx-1d0)**2*(x*(expx+1d0)/(expx-1d0)-4d0)
    
    return
  end function S_thermalsz

  function S_powerlaw(nu_ghz,nu0,beta) result (S)
    implicit none
    real*8 S
    real*8 nu_ghz, beta,nu0
    
    S = (nu_ghz/nu0)**beta

    return

  end function S_powerlaw

  function S_powerlawrunning(nu_ghz,nu0,beta,alpha) result (S)
    implicit none
    real*8 S,logS
    real*8 nu_ghz, beta,alpha,nu0
    
    logS = beta*log(nu_ghz/nu0)+0.5d0*alpha*log(nu_ghz/nu0)**2
    S    = exp(logS)
    
    return

  end function S_powerlawrunning

  function S_powerlaw3(nu_ghz,nu0,beta,alpha,gamma) result (S)
    implicit none
    real*8 S,logS
    real*8 nu_ghz, beta,alpha,gamma,nu0

    
    logS = beta*log(nu_ghz/nu0)+0.5d0*alpha*log(nu_ghz/nu0)**2+gamma/6d0*log(nu_ghz/nu0)**3
    S    = exp(logS)
    
    return

  end function S_powerlaw3

  function S_freefree(nu_ghz,nu0,temp) result (S)
    !Based on Dickinson, et al'03 Eq.(8) and (10)
    !NEEDS TESTING
    implicit none
    real*8 S,S0
    real*8 nu_ghz,nu0,temp
    real*8 const

    const = 1.5d0*dlog(temp)    
    S0    = nu0**(-2)*(dlog(4.995d-2/nu0)+const)
    S     = nu_ghz**(-2)*(dlog(4.995d-2/nu_ghz)+const)
    S     = S/S0
    
    return

  end function S_freefree
  
  function S_freemixingelements(freq,nfreq,fg_param,nchannel) result(S)
    real*8 freq(nfreq)
    integer nfreq,nchannel
    real*8 fg_param(nchannel-1)
    real*8 S(nfreq)
    
    integer ff,ii
    
    S(1)= 1d0    
    ii  = 0
    do ff= 2 ,nfreq
       if (freq(ff) .eq. freq(ff-1)) then
          S(ff) = S(ff-1)
       else
          ii    = ii+1
          S(ff) = fg_param(ii)
       endif
    enddo
  end function S_freemixingelements

  subroutine get_bands_with_lines(band,nband,nabscissa,line_ghz,nline,band_with_line,nband_with_line)
    real*8,        intent(in) :: band(nabscissa,nband)
    integer,       intent(in) :: nband
    integer,       intent(in) :: nabscissa
    real*8,        intent(in) :: line_ghz(nline)
    integer,       intent(in) :: nline
    integer                   :: band_with_line(nband)
    integer                   :: nband_with_line
    
    integer bb,ll
    
    nband_with_line = 0
    do ll = 1, nline
       do bb = 1, nband
          if( (line_ghz(ll) .gt. band(1,bb)) .and. (line_ghz(ll) .lt. band(nabscissa,bb)) ) then
             nband_with_line                 = nband_with_line + 1
             band_with_line(nband_with_line) = bb  !Need to rethink this for multiple bands at same freq
          endif          
       end do
    end do

  end subroutine get_bands_with_lines
  
  
  function effective_frequency(numin,numax,powerlawslope) result (f)
    implicit none
   
    real f,numin,numax,powerlawslope
    
    f = 1d0/(powerlawslope+1d0)*&
         (numax**(powerlawslope+1d0)-numax**(powerlawslope+1d0)/&
         (numax-numin))**(1d0/powerlawslope)

    return
  end function effective_frequency


function get_mixingmatrix(component_list,reference_frequency_list,ncomp,&
     fg_param,band,weight,nband,nabscissa) result(pmat)
  
  character(32), intent(in) :: component_list(ncomp)
  real*8,        intent(in) :: reference_frequency_list(ncomp)
  integer,       intent(in) :: ncomp
  real*8,        intent(in) :: fg_param(:)
  real*8,        intent(in) :: band(nabscissa,nband)
  real*8,        intent(in) :: weight(nabscissa,nband)
  integer,       intent(in) :: nband
  integer,       intent(in) :: nabscissa
  real(dp)                     pmat(nband,ncomp)

  integer bb,ff,cc,fg_param_index
  integer                   :: band_with_line(nband)
  integer                   :: nband_with_line
  real(dp)                  :: response(nabscissa)


  pmat(:,:)      = 0.0
  fg_param_index = 1
  do cc = 1, ncomp
     select case(component_list(cc))
     case ('tabulated')
        do bb = 1, nband
           call interp_linear(1,tabulated_scaling_numpoint(cc),&
                tabulated_scaling(cc,:,1),tabulated_scaling(cc,:,2),&
                nabscissa,band(:,bb),response)
           do ff = 1, nabscissa
              pmat(bb,cc)= pmat(bb,cc) + response(ff)*weight(ff,bb)
           end do
        end do
     case ('absolute')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)= pmat(bb,cc) + S_absolute(band(ff,bb))*weight(ff,bb)
           end do
        end do
     case ('cmb')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)= pmat(bb,cc) + S_cmb(band(ff,bb))*weight(ff,bb)
           end do
        end do
     case ('dq')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)= pmat(bb,cc) + S_dq(band(ff,bb))*weight(ff,bb)
           end do
        end do
     case ('thermalsz')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)= pmat(bb,cc) + S_thermalsz(band(ff,bb))*weight(ff,bb)
           end do
        end do
     case ('fdsmodel1')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)= pmat(bb,cc) + S_greybody(band(ff,bb),reference_frequency_list(cc),&
                   alpha1_FDSMODEL1,T1_FDSMODEL1)*weight(ff,bb)
           end do
        end do
     case ('fdsmodel2')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)= pmat(bb,cc) + S_greybody(band(ff,bb),reference_frequency_list(cc),&
                   alpha1_FDSMODEL2,T1_FDSMODEL2)*weight(ff,bb)
           end do
        end do
     case ('fdsmodel3')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)= pmat(bb,cc) + S_greybody(band(ff,bb),reference_frequency_list(cc),&
                   alpha1_FDSMODEL3,T1_FDSMODEL3)*weight(ff,bb)
           end do
        end do
     case ('fdsmodel4')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)= pmat(bb,cc) + S_greybody(band(ff,bb),reference_frequency_list(cc),&
                   alpha1_FDSMODEL4,T1_FDSMODEL4)*weight(ff,bb)
           end do
        end do
     case ('fdsmodel5')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)= pmat(bb,cc) + S_dust_fds(band(ff,bb),reference_frequency_list(cc),&
                   alpha1_FDSMODEL5,t1_FDSMODEL5,&
                   alpha2_FDSMODEL5,t2_FDSMODEL5,f1_FDSMODEL5,q1perq2_FDSMODEL5)*weight(ff,bb)
           end do
        end do
     case ('fdsmodel6')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)= pmat(bb,cc) + S_dust_fds(band(ff,bb),reference_frequency_list(cc),&
                   alpha1_FDSMODEL6,t1_FDSMODEL6,&
                   alpha2_FDSMODEL6,t2_FDSMODEL6,f1_FDSMODEL6,q1perq2_FDSMODEL6)*weight(ff,bb)
           end do
        end do
     case ('fdsmodel7')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)= pmat(bb,cc) + S_dust_fds(band(ff,bb),reference_frequency_list(cc),&
                   alpha1_FDSMODEL7,t1_FDSMODEL7,&
                   alpha2_FDSMODEL7,t2_FDSMODEL7,f1_FDSMODEL7,q1perq2_FDSMODEL7)*weight(ff,bb)
           end do
        end do
     case ('fdsmodel8')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)= pmat(bb,cc) + S_dust_fds(band(ff,bb),reference_frequency_list(cc),&
                   alpha1_FDSMODEL8,t1_FDSMODEL8,&
                alpha2_FDSMODEL8,t2_FDSMODEL8,f1_FDSMODEL8,q1perq2_FDSMODEL8)*weight(ff,bb)
           end do
        end do
     case ('twocomponent')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc) = pmat(bb,cc) + S_greybody(band(ff,bb),reference_frequency_list(cc),&
                   fg_param(fg_param_index),fg_param(fg_param_index+1))*weight(ff,bb) + &
                   S_greybody(band(ff,bb),reference_frequency_list(cc),&
                   fg_param(fg_param_index+2),fg_param(fg_param_index+3))*weight(ff,bb)*&
                   fg_param(fg_param_index+4)
           end do
        end do
        fg_param_index=fg_param_index+2
     case ('greybody')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)= pmat(bb,cc) + S_greybody(band(ff,bb),reference_frequency_list(cc),&
                   fg_param(fg_param_index),fg_param(fg_param_index+1))*weight(ff,bb)
           end do
        end do
        fg_param_index=fg_param_index+2
     case ('greybody18K')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)= pmat(bb,cc) + S_greybody(band(ff,bb),reference_frequency_list(cc),&
                   fg_param(fg_param_index),18d0)*weight(ff,bb)
           end do
        end do
        fg_param_index=fg_param_index+1
     case ('greybody21K')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)= pmat(bb,cc) + S_greybody(band(ff,bb),reference_frequency_list(cc),&
                   fg_param(fg_param_index),21d0)*weight(ff,bb)
           end do
        end do
        fg_param_index=fg_param_index+1
     case ('greybody_dbeta')
        !Use in combination with the greybody - requires the same dust temperature.
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)= pmat(bb,cc) + S_greybody_dbeta(band(ff,bb),reference_frequency_list(cc),&
                   fg_param(fg_param_index-2),fg_param(fg_param_index-1))*weight(ff,bb)
           end do
        end do
     case ('greybody_dtemp')
        !Use in combination with the greybody - requires the same dust spectral index.
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)= pmat(bb,cc) + S_greybody_dtemp(band(ff,bb),reference_frequency_list(cc),&
                   fg_param(fg_param_index-2),fg_param(fg_param_index-1))*weight(ff,bb)
           end do
        end do
     case ('sync')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)= pmat(bb,cc) + S_powerlaw(band(ff,bb),reference_frequency_list(cc),-3d0)*weight(ff,bb)
           end do
        end do
     case ('freefree')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)= pmat(bb,cc) + S_freefree(band(ff,bb),reference_frequency_list(cc),T_E)*weight(ff,bb)
           end do
        end do
        !May need to implement beta_FF as per Eq.3 of MAMD et al '08        
     case ('freefree2.14')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)= pmat(bb,cc) + S_powerlaw(band(ff,bb),reference_frequency_list(cc),-2.14d0)*weight(ff,bb)
           end do
        end do
        !May need to implement beta_FF as per Eq.3 of MAMD et al '08        
     case ('sd')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)= pmat(bb,cc) + S_greybody(band(ff,bb),reference_frequency_list(cc),1.6d0,&
                   fg_param(fg_param_index)*h_over_k)*weight(ff,bb)
           end do
        end do
        fg_param_index=fg_param_index+1
     case ('powerlaw')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)= pmat(bb,cc) + S_powerlaw(band(ff,bb),reference_frequency_list(cc),&
                   fg_param(fg_param_index))*weight(ff,bb)
           end do
        end do
        fg_param_index=fg_param_index+1
     case ('powerlaw_dbeta')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)=  pmat(bb,cc) + log(band(ff,bb)/reference_frequency_list(cc))*&
                   S_powerlaw(band(ff,bb),reference_frequency_list(cc),fg_param(fg_param_index-1))*weight(ff,bb)
           end do
        end do
     case ('powerlaw_dbeta2')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)=  pmat(bb,cc) + log(band(ff,bb)/reference_frequency_list(cc))**2*&
                   S_powerlaw(band(ff,bb),reference_frequency_list(cc),fg_param(fg_param_index-1))**2*weight(ff,bb)
           end do
        end do
     case ('powerlaw_dbeta3')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)=  pmat(bb,cc) + log(band(ff,bb)/reference_frequency_list(cc))**3*&
                   S_powerlaw(band(ff,bb),reference_frequency_list(cc),&
                   fg_param(fg_param_index-1))**2*weight(ff,bb)
           end do
        end do
     case ('powerlawrunning')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)= pmat(bb,cc) + S_powerlawrunning(band(ff,bb),reference_frequency_list(cc),&
                   fg_param(fg_param_index),fg_param(fg_param_index+1))*weight(ff,bb)
           end do
        end do
        fg_param_index=fg_param_index+2
     case ('powerlaw3')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)= pmat(bb,cc) + S_powerlaw3(band(ff,bb),reference_frequency_list(cc), &
                   fg_param(fg_param_index),fg_param(fg_param_index+1),&
                   fg_param(fg_param_index+2))*weight(ff,bb)
           end do
        end do
        fg_param_index=fg_param_index+3
     case ('constant')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)= pmat(bb,cc) +  1d0*weight(ff,bb)
           end do
        end do
     case ('linear')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)= pmat(bb,cc) +  dlog(band(ff,bb)/reference_frequency_list(cc))*weight(ff,bb)
           end do
        end do
     case ('quadratic')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)= pmat(bb,cc) +  dlog(band(ff,bb)/reference_frequency_list(cc))**2*weight(ff,bb)
           end do
        end do
     case ('quartic')
        do bb=1, nband
           do ff=1, nabscissa
              pmat(bb,cc)= pmat(bb,cc) +  dlog(band(ff,bb)/reference_frequency_list(cc))**4*weight(ff,bb)
           end do
        end do
     case ('co')
        call get_bands_with_lines(band,nband,nabscissa,[CO_line1,CO_line2],CO_nline,&
             band_with_line,nband_with_line)
        if( (nband_with_line .ge. 1)) then
           pmat(band_with_line(1),cc) = pmat(band_with_line(1),cc) + 1d0
           do bb = 2, nband_with_line
              pmat(band_with_line(bb),cc) = pmat(band_with_line(bb),cc) + fg_param(fg_param_index+bb-2)           
           end do
           fg_param_index = fg_param_index + nband_with_line - 1
        endif
     case default
     end select
  end do

 end function get_mixingmatrix

 function get_totalforegroundparameters(component_list,ncomp) result(nparam)
   character(32), intent(in) :: component_list(ncomp)
   integer, intent(in)       :: ncomp
   integer cc,nparam

   nparam = 0
   do cc=1,ncomp
      select case(component_list(cc))
      case ('tabulated')
         continue
      case ('absolute')
         continue
      case ('cmb')
         continue
      case ('dq')
         continue
      case ('thermalsz')
         continue
      case ('fdsmodel1')
         continue
      case ('fdsmodel2')
         continue
      case ('fdsmodel3')
         continue
      case ('fdsmodel4')
         continue
      case ('fdsmodel5')
         continue
      case ('fdsmodel6')
         continue
      case ('fdsmodel7')
         continue
      case ('fdsmodel8')
         continue
      case ('greybody')
         nparam = nparam+2
      case ('greybody18K')
         nparam = nparam+1
      case ('greybody21K')
         nparam = nparam+1
      case ('twocomponent')
         nparam = nparam+5
      case ('greybody_dbeta')
         continue
      case ('greybody_dtemp')
         continue
      case ('sync')
         continue
      case ('freefree')
         continue
      case ('freefree2.14')
         continue
      case ('sd')
         nparam = nparam+1
      case ('powerlaw')
         nparam = nparam+1
      case ('powerlaw_dbeta')
         continue
      case ('powerlaw_dbeta2')
         continue
      case ('powerlaw_dbeta3')
         continue
      case ('powerlawrunning')
         nparam = nparam+2
      case ('powerlaw3')
         nparam = nparam+3
      case ('constant')
         continue
      case ('linear')
         continue
      case ('quadratic')
         continue
      case ('quartic')
         continue
      case ('co')
         nparam = nparam + CO_nline - 1
      case default
         print *,'Component named:', component_list(cc),' not implemented.'
         stop
      end select
   end do

 end function get_totalforegroundparameters


 function get_mixingmatrix_partialderivative(component_list,reference_frequency_list,&
      ncomp,fg_param,freq,nfreq,fg_param_ii)  result(pprimemat)
   character(32), intent(in) :: component_list(ncomp)
   real*8,        intent(in) :: reference_frequency_list(ncomp)
   integer,       intent(in) :: ncomp
   real(dp),      intent(in) :: fg_param(:)
   real(dp),      intent(in) :: freq(nfreq)
   integer,       intent(in) :: nfreq,fg_param_ii
   real(dp) pprimemat(nfreq,ncomp)

   real(dp) pprimecol(nfreq),delta,deltap(nfreq,1)
   integer cc,ii,nparam_thiscomponent,index,index2
   real(dp), dimension(:), allocatable     :: local_fg_param(:)
   real(dp), dimension(:,:), allocatable   :: weight
   integer nabscissa

   allocate(weight(1,nfreq))

   !------------------------------------------------------------------------------------
   !Calculates partial derivative of mixing matrix with respect to parameter fg_param_ii
   !------------------------------------------------------------------------------------
   
   pprimemat(:,:) = 0.0d0
   delta          = 0.01d0
   
   index       = 0
   cc          = 0
   nabscissa   = 1
   weight(:,:) = 1.
   do while (index .lt. fg_param_ii ) 
      cc                   = cc + 1
      nparam_thiscomponent = get_totalforegroundparameters(component_list(cc),1)
      if (nparam_thiscomponent .eq. 0) then
         !-------------------------------------------------------------------------------
         !Components with no parameters (like CMB) have partial derivatives equal to zero
         !-------------------------------------------------------------------------------
         continue
      else
         index = index + nparam_thiscomponent
         if (index .lt. fg_param_ii) then
            continue
         else
            allocate(local_fg_param(nparam_thiscomponent))
            !Get mixing matrix column
            pprimecol(:) = 0.0d0
            index2       = nparam_thiscomponent - (index - fg_param_ii)
            do ii=1,2
               local_fg_param(1:nparam_thiscomponent) = &
                    fg_param(index-nparam_thiscomponent+1:index)
               !Perturb
               local_fg_param(index2) = local_fg_param(index2) + &
                    delta*(-1.d0)**ii*fg_param(fg_param_ii)
               
               deltap = get_mixingmatrix(component_list(cc),&
                    reference_frequency_list(cc),&
                    1,local_fg_param,freq,weight,nfreq,nabscissa)
               pprimecol(1:nfreq) = pprimecol(1:nfreq) + (-1d0)**ii*deltap(1:nfreq,1)
            end do
            !----------------------------------------
            !Simple approximation of first derivative
            !----------------------------------------
            pprimecol(1:nfreq) = pprimecol(1:nfreq)/2./delta/fg_param(fg_param_ii)         
         endif
      endif      
   end do

   pprimemat(1:nfreq,cc) = pprimecol(1:nfreq) 
   deallocate(local_fg_param)

 end function get_mixingmatrix_partialderivative
  
end module foregroundscalings


