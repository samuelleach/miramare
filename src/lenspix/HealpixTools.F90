module HealpixTools
!Extra tools by SML built around HealpixObj

  use HealpixObj
  use mpi_stop
  use Matrixutils

contains

  subroutine HealpixMap_invert(map,covmat_firstcol)
    
    type(HealpixMap), intent(inout) ::  map
    integer, intent(in), optional   :: covmat_firstcol
    real(dp), dimension(3,3) :: mat,testmat
    integer nelement
    integer offset
    logical mpd
    integer*8 pp
    integer*8 npd
    integer mykind 

    if (.not. present(covmat_firstcol)) then 
       offset = 0
    else
       offset = covmat_firstcol-1
    endif

    mykind   = kind(map%tqu(0,1+offset))
    nelement = map%nmaps - offset
    npd      = 0
    do pp = 0, map%npix-1
       if (map%tqu(pp,1+offset) .ne. 0.) then    ! Quick check on first diagonal component.
          select case(nelement)
          case(1)
             !Single map
             map%tqu(pp,1+offset) = 1./map%tqu(pp,1+offset)
          case(2)
             !Assume two maps are hit map and variance map
             map%tqu(pp,1+offset) = 1./map%tqu(pp,1+offset)
             map%tqu(pp,2+offset) = 1./map%tqu(pp,2+offset)
          case(3)
             !Three element map (assumed diagonal covariance)
             map%tqu(pp,1+offset) = 1./map%tqu(pp,1+offset)
             map%tqu(pp,2+offset) = 1./map%tqu(pp,2+offset)
             map%tqu(pp,3+offset) = 1./map%tqu(pp,3+offset)
          case(4)
             !Assume WMAP ordering TT,QQ,QU,UU (nhits)
             stop
          case(6:7)
             !Assume "Madam" and "Springtide" ordering TT,TQ,TU,QQ,QU,UU (nhits)
             mat(1,1) = map%tqu(pp,1+offset)
             mat(2,2) = map%tqu(pp,4+offset)
             mat(3,3) = map%tqu(pp,6+offset)
             mat(1,2) = map%tqu(pp,2+offset)
             mat(2,1) = mat(1,2)
             mat(1,3) = map%tqu(pp,3+offset)
             mat(3,1) = mat(1,3)
             mat(2,3) = map%tqu(pp,5+offset)
             mat(3,2) = mat(2,3)
             testmat  = mat
             mpd      = matrix_is_positive_definite(testmat)
             if (mpd) then
                npd = npd +1
                call Matrix_inverse(mat)
                map%tqu(pp,1+offset) = real(mat(1,1))
                map%tqu(pp,2+offset) = real(mat(1,2))
                map%tqu(pp,3+offset) = real(mat(1,3))
                map%tqu(pp,4+offset) = real(mat(2,2))
                map%tqu(pp,5+offset) = real(mat(2,3))
                map%tqu(pp,6+offset) = real(mat(3,3))
             else
                print *, 'Warning: pixel is not positive definite, pixel = ',pp
                map%tqu(pp,1+offset) = real(1d0/map%tqu(pp,1+offset)) !Kluge
                map%tqu(pp,2+offset) = 0.
                map%tqu(pp,3+offset) = 0.
                map%tqu(pp,4+offset) = 0.
                map%tqu(pp,5+offset) = 0.
                map%tqu(pp,6+offset) = 0.
             endif
          case default
             print *,'Nelement = ',nelement
             call doStop('This case is not handled in HealpixTools.F90:HealpixMap_invert')
          end select
       endif
    end do
    
  end subroutine HealpixMap_invert

  subroutine HealpixMap_gausslike(map,ninv,chi2,ninv_firstcol)
    
    !Purpose: Calculate map . ninv . map

    type(HealpixMap), intent(inout) ::  map
    type(HealpixMap), intent(in)    ::  ninv
    type(HealpixMap), intent(inout) ::  chi2
    integer, intent(in), optional   ::  ninv_firstcol

    integer pp
    integer offset
    real(dp), dimension(3,3) :: mat
    real(dp), dimension(3) :: tmpmap

    if (.not. present(ninv_firstcol)) then 
       offset = 0
    else
       offset = ninv_firstcol-1
    endif

    do pp = 0, map%npix -1        
       mat(1,1) = ninv%tqu(pp,1+offset)
       mat(2,2) = ninv%tqu(pp,4+offset)
       mat(3,3) = ninv%tqu(pp,6+offset)
       mat(1,2) = ninv%tqu(pp,2+offset)
       mat(2,1) = mat(1,2)
       mat(1,3) = ninv%tqu(pp,3+offset)
       mat(3,1) = mat(1,3)
       mat(2,3) = ninv%tqu(pp,5+offset)
       mat(3,2) = mat(2,3)

       tmpmap           = matmul(mat,map%tqu(pp,1:3))
       chi2%tqu(pp,1)   = tmpmap(1)*map%tqu(pp,1)
       chi2%tqu(pp,2)   = tmpmap(2)*map%tqu(pp,2)
       chi2%tqu(pp,3)   = tmpmap(3)*map%tqu(pp,3)

       chi2%tqu(pp,4) = chi2%tqu(pp,1) + chi2%tqu(pp,2) + chi2%tqu(pp,3)
       
    end do

  end subroutine HealpixMap_gausslike
  
  function HealpixMap_getmedian(map,missval,mapindex)
    !Purpose: Get the median value of a Healpix Map.
    
    
    type(HealpixMap), intent(in)  :: map
    real, intent(in), optional    :: missval
    integer, intent(in), optional :: mapindex
    real HealpixMap_getmedian
    
    real, dimension(:), allocatable :: x
    real                            :: mymissval
    integer                         :: mymapindex, pp, npix

    if(present(missval)) then
       mymissval = missval
    else
       mymissval = fmissval
    endif

    if(present(mapindex)) then
       mymapindex = mapindex
    else
       mymapindex = 1
    endif

    npix = 0
    do pp = 0, map%npix - 1
!       if( abs((map%tqu(pp,mymapindex) - mymissval)/mymissval) .ge. 1e-5) then
       if( map%tqu(pp,mymapindex) .ne. mymissval) then
          npix = npix + 1
       end if
    end do


    allocate(x(npix))
    npix = 0
    do pp = 0,map%npix-1
!       if( abs((map%tqu(pp,mymapindex) - mymissval)/mymissval) .ge. 1e-5) then
       if( map%tqu(pp,mymapindex) .ne. mymissval) then
          npix    = npix + 1
          x(npix) = map%tqu(pp,mymapindex)
       end if
    end do


    call median(x,npix,HealpixMap_getmedian)

    deallocate(x)


  end function HealpixMap_getmedian


  


end module HealpixTools
