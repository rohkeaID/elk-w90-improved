
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine zgerci(m,n,alpha,x,y,ld,a)
implicit none
! arguments
integer, intent(in) :: m,n
complex(8), intent(in) :: alpha
complex(8), intent(in) :: x(m),y(n)
integer, intent(in) :: ld
complex(8), intent(inout) :: a(ld,*)
! local variables
integer j
! numbers less than eps are considered to be zero
real(8), parameter :: eps=1.d-10
real(8) a1,b1
complex(8) z1
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(z1,a1,b1)
!$OMP DO
do j=1,n
  z1=alpha*y(j)
  if (abs(dble(z1)).gt.eps) then
    if (abs(aimag(z1)).gt.eps) then
! complex prefactor
      call zaxpy(m,z1,x,1,a(:,j),1)
    else
! real prefactor
      a1=dble(z1)
      call daxpy(2*m,a1,x,1,a(:,j),1)
    end if
  else if (abs(aimag(z1)).gt.eps) then
! imaginary prefactor
    b1=aimag(z1)
    a(1:m,j)=a(1:m,j)+b1*cmplx(-aimag(x(1:m)),dble(x(1:m)),8)
  end if
end do
!$OMP END DO
!$OMP END PARALLEL
return
end subroutine

