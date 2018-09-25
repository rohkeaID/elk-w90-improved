
! Copyright (C) 2017 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine acpole(zfm,zwr,zr)
use modmain
use modgw
implicit none
! arguments
complex(8), intent(in) :: zfm(0:nwfm),zwr(nwplot)
complex(8), intent(out) :: zr(nwplot)
! local variables
integer, parameter :: maxit=1000
integer iter,iw
integer n,n2,i,j
real(8), parameter :: eps=1.d-5
! allocatable arrays
real(8), allocatable :: x(:,:)
! external functions
complex(8) zfpole
external zfpole
n=2*npole+1
n2=2*n
allocate(x(n2,n2+1))
! intialise simplex
x(:,:)=0.d0
! fit the constant
x(1,2)=0.5d0
x(2,3)=0.5d0
call minf_nm(1,zfm,n2,x,maxit,iter,eps)
! fit the constant and the first pole
x(3,1)=1.d0
do i=1,6
  x(i,i+1)=x(i,1)+0.1d0
end do
call minf_nm(1,zfm,n2,x,maxit,iter,eps)
! fit the remaining poles one-by-one
i=7
do j=2,npole
  x(i,1)=1.d0
  x(i,i+1)=x(i,1)+0.1d0
  i=i+1
  x(i,i+1)=0.1d0
  i=i+1
  x(i,i+1)=0.1d0
  i=i+1
  x(i,i+1)=0.1d0
  i=i+1
  call minf_nm(1,zfm,n2,x,maxit,iter,eps)
end do
! fit the constant and the first pole again
do i=1,6
  x(i,i+1)=x(i,1)+0.1d0
end do
call minf_nm(1,zfm,n2,x,maxit,iter,eps)
! fit everything together
if (npole.gt.1) then
  do i=1,n2
    x(i,i+1)=x(i,1)+0.1d0
  end do
  call minf_nm(1,zfm,n2,x,maxit,iter,eps)
end if
do iw=1,nwplot
  zr(iw)=zfpole(x(:,1),zwr(iw))
end do
deallocate(x)
return
end subroutine

