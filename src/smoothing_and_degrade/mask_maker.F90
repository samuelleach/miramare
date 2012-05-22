module mask_maker

  use HealpixObj
  
contains
  
  
  subroutine get_mask_from_map(Map,Mask,mask_npix,nside_out)
    
    Type(HealpixMap), intent(in)   :: Map
    Type(HealpixMap), intent(out)  :: Mask
    integer*8, intent(out)         :: mask_npix
    integer, intent(in), optional  :: nside_out

    integer npix
    integer*8 pp
    Type(HealpixMap) :: Mask_temp
    
    character(len=*), parameter :: CODE = "GET_MASK_FROM_MAP"


    mask%TQU(:,1) = 1.
    mask%ordering = map%ordering

    do pp = 0, map%npix-1
       if ( abs(map%TQU(pp,1)/fmissval -1.) .lt. 1.e-5 ) mask%TQU(pp,1) = 0.
    enddo

    if(present(nside_out)) then
       if(nside_out .lt. map%nside) then
          npix = nside2npix(nside_out)
          call HealpixMap_init(mask_temp,npix,nmaps=1)
          call HealpixMap_udgrade(mask,mask_temp,nside_out)
          mask = mask_temp          
          where( mask%TQU(:,1).lt. 1. )  mask%TQU(:,1) = 0.
       end if
    endif
    
    mask_npix = 0
    do pp = 0, mask%npix-1
       if ( mask%TQU(pp,1) .lt. 1.e-5 ) mask_npix = mask_npix + 1
    enddo
    

  end subroutine get_mask_from_map
  

  subroutine smooth_mask(H,mask,beam,nlmax,nside_out)
    
    Type(HealpixInfo)               :: H
    Type(HealpixMap)                :: Mask
    real(dp), pointer               :: beam(:)
    integer nlmax
    integer, intent(in), optional   :: nside_out

    Type(HealpixMap)                :: Mask_Temp
    Type(HealpixAlm)                :: A
    integer npix
    logical verbose

    character(len=*), parameter :: CODE = "SMOOTH_MASK"
    verbose = .true.

    if(verbose) print *,trim(CODE),': map2alm'
    call HealpixMap2alm(H,mask,A,nlmax)!,dopol=.false.)
    
    if(verbose) print *,trim(CODE),': smoothing alm'
    call HealpixAlm_smooth_beam(A,beam)
    
    if(verbose) print *,trim(CODE),': alm2map '
    call HealpixAlm2Map(H,A, mask, mask%npix) 

    if(present(nside_out)) then
       if(nside_out .lt. mask%nside) then
          npix = nside2npix(nside_out)
          call HealpixMap_init(mask_temp,npix,nmaps=1)
          call HealpixMap_udgrade(mask,mask_temp,nside_out)
          mask = mask_temp          
       end if
    endif

    where( mask%TQU(:,1).lt. 1. )  mask%TQU(:,1) = 0.


  
  end subroutine smooth_mask
  
end module mask_maker






