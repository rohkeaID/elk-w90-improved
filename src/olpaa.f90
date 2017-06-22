
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine olpaa(ias,ngp,apwalm,ld,o)
use modmain
implicit none
! arguments
integer, intent(in) :: ias
integer, intent(in) :: ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
integer, intent(in) :: ld
complex(8), intent(inout) :: o(*)
! local variables
integer is,l,m,lm,io
! automatic arrays
complex(8) x(ngp)
is=idxis(ias)
lm=0
do l=0,lmaxmat
  do m=-l,l
    lm=lm+1
    do io=1,apword(l,is)
      x(1:ngp)=conjg(apwalm(1:ngp,io,lm,ias))
      call zheri(ngp,1.d0,x,ld,o)
    end do
  end do
end do
return

contains

subroutine zheri(n,alpha,x,ld,a)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: alpha
complex(8), intent(in) :: x(n)
integer, intent(in) :: ld
complex(8), intent(inout) :: a(*)
! local variables
integer j,k
! numbers less than eps are considered to be zero
real(8), parameter :: eps=1.d-10
real(8) a1,b1
complex(8) z1
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,z1,a1,b1)
!$OMP DO
do j=1,n
  k=(j-1)*ld
  z1=alpha*conjg(x(j))
  if (abs(dble(z1)).gt.eps) then
    if (abs(aimag(z1)).gt.eps) then
! complex prefactor
      call zaxpy(j-1,z1,x,1,a(k+1),1)
      a(k+j)=dble(a(k+j))+dble(z1*x(j))
    else
! real prefactor
      a1=dble(z1)
      call daxpy(2*(j-1),a1,x,1,a(k+1),1)
      a(k+j)=dble(a(k+j))+a1*dble(x(j))
    end if
  else if (abs(aimag(z1)).gt.eps) then
! imaginary prefactor
    b1=aimag(z1)
    a(k+1:k+j-1)=a(k+1:k+j-1)+b1*cmplx(-aimag(x(1:j-1)),dble(x(1:j-1)),8)
    a(k+j)=dble(a(k+j))-b1*aimag(x(j))
  end if
end do
!$OMP END DO
!$OMP END PARALLEL
return
end subroutine

end subroutine

