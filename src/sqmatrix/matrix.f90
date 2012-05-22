MODULE matrix
!
! Matrix algebra
!
 implicit none

 CONTAINS

 Subroutine invert3_eig(cc,cdet,rcond,sing)
!
! Inversion of a symmetric positive-definite 3x3 matrix
!   through eigenvalue analysis.
!
   double precision,intent(inout) :: cc(6)
   double precision,intent(out)   :: cdet   ! determinant
   double precision,intent(out)   :: rcond  ! reciprocal condition number
   logical, intent(out),optional  :: sing   ! flag for singular matrix
   integer                        :: nrot
   double precision               :: a(3,3),v(3,3),d(3),di(3),dmin,dmax

   a(1,1) = cc(1)
   a(1,2) = cc(2)
   a(1,3) = cc(3)
   a(2,1) = cc(2)
   a(2,2) = cc(4)
   a(2,3) = cc(5)
   a(3,1) = cc(3)
   a(3,2) = cc(5)
   a(3,3) = cc(6)

   call Jacobi(a,3,d,v,nrot)

   dmin = minval(abs(d))
   dmax = maxval(abs(d))

   if (dmax.lt.1.e-30) then
      if (present(sing)) sing=.true.
      cc = 0.0
      cdet = 0.0
      rcond = 0.0
      return
   elseif (dmin.le.1.e-12*dmax) then
      if (present(sing)) sing=.true.
      cc = 0.0
      cdet = d(1)*d(2)*d(3)
      rcond = dmin/dmax
      return
   else
      if (present(sing)) sing=.false.
      di = 1/d
      cdet = d(1)*d(2)*d(3)
      rcond = dmin/dmax

      cc(1) = di(1)*v(1,1)*v(1,1)+di(2)*v(1,2)*v(1,2)+di(3)*v(1,3)*v(1,3)
      cc(2) = di(1)*v(1,1)*v(2,1)+di(2)*v(1,2)*v(2,2)+di(3)*v(1,3)*v(2,3)
      cc(3) = di(1)*v(1,1)*v(3,1)+di(2)*v(1,2)*v(3,2)+di(3)*v(1,3)*v(3,3)
      cc(4) = di(1)*v(2,1)*v(2,1)+di(2)*v(2,2)*v(2,2)+di(3)*v(2,3)*v(2,3)
      cc(5) = di(1)*v(2,1)*v(3,1)+di(2)*v(2,2)*v(3,2)+di(3)*v(2,3)*v(3,3)
      cc(6) = di(1)*v(3,1)*v(3,1)+di(2)*v(3,2)*v(3,2)+di(3)*v(3,3)*v(3,3)
   endif

 end subroutine


!-----------------------------------------------------------------------


 Subroutine invert3_eig_s(cc,cdet,rcond,noeig,plim)
!
! Inversion of a symmetric positive-definite 3x3 matrix
!   through eigenvalue analysis.
! If matrix is singular, invert the non-singular part.
!
   double precision,intent(inout) :: cc(6)
   double precision,intent(out)   :: cdet   ! determinant
   double precision,intent(out)   :: rcond  ! reciprocal condition number
   integer,intent(out)            :: noeig  ! number of non-zero eigenvalues
   double precision,intent(in)    :: plim   ! lowest accepted eigenvalue
   integer                        :: k,nrot
   double precision               :: a(3,3),v(3,3),d(3),di(3)
   double precision               :: dlim,dmin,dmax

   a(1,1) = cc(1)
   a(1,2) = cc(2)
   a(1,3) = cc(3)
   a(2,1) = cc(2)
   a(2,2) = cc(4)
   a(2,3) = cc(5)
   a(3,1) = cc(3)
   a(3,2) = cc(5)
   a(3,3) = cc(6)

   call Jacobi(a,3,d,v,nrot)

   dmin = minval(abs(d))
   dmax = maxval(abs(d))

   if (dmax.lt.1.e-30) then
      noeig = 0
      cc = 0.0
      cdet = 0.0
      rcond = 0.0
      return
   else
      noeig = 3
      dlim = max(plim*dmax,1.e-12*dmax)
      cdet = d(1)*d(2)*d(3)
      rcond = dmin/dmax
      do k = 1,3
         if (d(k).le.dlim) then
            di(k)=0.0
            noeig=noeig-1
         else
            di(k)=1/d(k)
         endif
      enddo

      cc(1) = di(1)*v(1,1)*v(1,1)+di(2)*v(1,2)*v(1,2)+di(3)*v(1,3)*v(1,3)
      cc(2) = di(1)*v(1,1)*v(2,1)+di(2)*v(1,2)*v(2,2)+di(3)*v(1,3)*v(2,3)
      cc(3) = di(1)*v(1,1)*v(3,1)+di(2)*v(1,2)*v(3,2)+di(3)*v(1,3)*v(3,3)
      cc(4) = di(1)*v(2,1)*v(2,1)+di(2)*v(2,2)*v(2,2)+di(3)*v(2,3)*v(2,3)
      cc(5) = di(1)*v(2,1)*v(3,1)+di(2)*v(2,2)*v(3,2)+di(3)*v(2,3)*v(3,3)
      cc(6) = di(1)*v(3,1)*v(3,1)+di(2)*v(3,2)*v(3,2)+di(3)*v(3,3)*v(3,3)
   endif

 end subroutine


!-----------------------------------------------------------------------


 Subroutine Jacobi(a,n,d,v,nrot)
 !
 ! Eigenvalues d and eigenvectors v of a matrix a.
 ! Jacobi algorithm from Numerical Recipes
 !
   integer,         intent(in)    :: n
   integer,         intent(out)   :: nrot
   double precision,intent(inout) :: a(n,n)
   double precision,intent(out)   :: d(n),v(n,n)
   integer                        :: i,j,ip,iq
   double precision               :: c,g,h,s,sm,t,tau,theta,tresh
   double precision               :: b(n),z(n)

   v=0.0
   do ip=1,n
      v(ip,ip)=1.d0
      d(ip)=a(ip,ip)
   enddo
   b=d
   z=0.0

   nrot=0
   do i=1,50

      sm=0.0
      do ip=1,n-1
         do iq=ip+1,n
            sm=sm+abs(a(ip,iq))
         enddo
      enddo
      if (sm.eq.0) return
      if (i.lt.4) then
         tresh=0.2*sm/(n*n)
      else
         tresh=0.0
      endif

      do ip=1,n-1
         do iq=ip+1,n
            g=100.*abs(a(ip,iq))
            if ((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.  &
                            (abs(d(iq))+g.eq.abs(d(iq)))) then
               a(ip,iq)=0.
            elseif (abs(a(ip,iq)).gt.tresh) then
               h=d(iq)-d(ip)
               if (abs(h)+g.eq.abs(h)) then
                  t=a(ip,iq)/h
               else
                  theta=0.5d0*h/a(ip,iq)
                  t=1/(abs(theta)+sqrt(1+theta*theta))
                  if (theta.lt.0) t=-t
               endif
               c=1/sqrt(1+t*t)
               s=t*c
               tau=s/(1+c)
               h=t*a(ip,iq)
               z(ip)=z(ip)-h
               z(iq)=z(iq)+h
               d(ip)=d(ip)-h
               d(iq)=d(iq)+h
               a(ip,iq)=0.
               do j=1,ip-1
                  g=a(j,ip)
                  h=a(j,iq)
                  a(j,ip)=g-s*(h+g*tau)
                  a(j,iq)=h+s*(g-h*tau)
               enddo
               do j=ip+1,iq-1
                  g=a(ip,j)
                  h=a(j,iq)
                  a(ip,j)=g-s*(h+g*tau)
                  a(j,iq)=h+s*(g-h*tau)
               enddo
               do j=iq+1,n
                  g=a(ip,j)
                  h=a(iq,j)
                  a(ip,j)=g-s*(h+g*tau)
                  a(iq,j)=h+s*(g-h*tau)
               enddo
               do j=1,n
                  g=v(j,ip)
                  h=v(j,iq)
                  v(j,ip)=g-s*(h+g*tau)
                  v(j,iq)=h+s*(g-h*tau)
               enddo
               nrot=nrot+1
            endif
         enddo
      enddo
      b=b+z
      d=b
      z=0.
   enddo
   write(*,*) 'Too many iterations'
   stop

 END SUBROUTINE

!------------------------------------------------------------------------


 SUBROUTINE invert3_LU(cc,cdet,sing)
! Invert a symmetric (3,3) matrix by LU decomposition
   double precision,intent(inout) :: cc(6)
   double precision,intent(out)   :: cdet
   logical, intent(out),optional  :: sing
   double precision               :: alu(3,3)

   alu(1,1) = cc(1)
   alu(1,2) = cc(2)
   alu(1,3) = cc(3)
   alu(2,1) = cc(2)
   alu(2,2) = cc(4)
   alu(2,3) = cc(5)
   alu(3,1) = cc(3)
   alu(3,2) = cc(5)
   alu(3,3) = cc(6)

   cc = 0.0
   cdet = 0.0
   if (present(sing)) sing = .true.

   if (abs(alu(1,1)).lt.1.e-30) return
   alu(2,1) = alu(2,1)/alu(1,1)
   alu(3,1) = alu(3,1)/alu(1,1)
   alu(2,2) = alu(2,2)-alu(2,1)*alu(1,2)
   if (abs(alu(2,2)).lt.1.e-30) return
   alu(3,2) = (alu(3,2)-alu(3,1)*alu(1,2))/alu(2,2)
   alu(2,3) = alu(2,3)-alu(2,1)*alu(1,3)
   alu(3,3) = alu(3,3)-alu(3,1)*alu(1,3)-alu(3,2)*alu(2,3)
   if (abs(alu(3,3)).lt.1.e-30) return

   cdet = alu(1,1)*alu(2,2)*alu(3,3)
   if (present(sing)) sing = .false.

   cc(3) = (-alu(3,1)+alu(3,2)*alu(2,1))/alu(3,3)
   cc(2) = (-alu(2,1)-alu(2,3)*cc(3))/alu(2,2)
   cc(1) = (1.d0-alu(1,2)*cc(2)-alu(1,3)*cc(3))/alu(1,1)
   cc(5) = -alu(3,2)/alu(3,3)
   cc(4) = (1.d0-alu(2,3)*cc(5))/alu(2,2)
   cc(6) = 1.d0/alu(3,3)

 END SUBROUTINE


!-----------------------------------------------------------------------

END MODULE
