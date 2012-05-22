module CalcLike
 use Random
 use settings
 use ParamDef
 use Likelihoodcmbforeground
 use runparameters
 implicit none

 real :: Temperature  = 1

contains
  
  function GetLogPrior(Params) !Get -Ln(Prior)
    type(ParamSet)  Params 
    real GetLogPrior
    integer ii

    GetLogPrior = 0.
    do ii=1,num_params_used
      if(scales%pvar(ii) .gt. 0.0) GetLogprior=GetLogprior+ &
           (params%p(ii)-scales%pmean(ii))**2/scales%pvar(ii)
    enddo
    GetLogPrior = Getlogprior*0.5
    
  end function GetLogPrior

  function GuessLnLikeOffset(Params)
    type(ParamSet)  Params 
    real GuessLnLikeOffset
    
      rp%lnlike_offset  = 0.
      GuessLnLikeOffset = -GetLogLike(Params)
      LastParams        = Params
    
    return
  end function


  function GetLogLike(Params) !Get -Ln(Likelihood)
    type(ParamSet)  Params 
    real GetLogLike
 
    if (any(Params%P > Scales%PMax) .or. any(Params%P < Scales%PMin)) then
       GetLogLike = logZero
       return
    end if

    select case(rp%likelihood_label)
       case(4)
          ! Spectral index likelihood.
          GetLogLike =  cmbforeground_LnLike4(params) 
       case(5)
          ! Spectral index likelihood with assumption of homogeneous noise
          GetLogLike =  cmbforeground_LnLike5(params) 
       case(7)
          ! Spectral index likelihood with calibration errors.
          GetLogLike =  cmbforeground_LnLike7(params) 
       case(8)
          ! Spectral index likelihood with calibration errors and assumption of homogeneous noise.
          GetLogLike =  cmbforeground_LnLike8(params) 
       case(9)
          ! Spectral index likelihood with calibration errors and offset marginalisation.
          GetLogLike =  cmbforeground_LnLike10(params) 
       case(11)
          ! 3x3 or 2x2 noise matrix
          GetLogLike =  cmbforeground_LnLike11(params) 
    end select

    LastParams=Params

   end function GetLogLike

    
  function GetLogLikePost() 
    real GetLogLikePost

    if ( GetLogLikePost >= logZero) then
       GetLogLikePost = logZero
    else 
    endif

     if (GetLogLikePost /= LogZero) GetLogLikePost = GetLogLikePost/Temperature

  end function GetLogLikePost

end module CalcLike
