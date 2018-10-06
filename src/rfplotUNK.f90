
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rfplotUNK(np,vpl,rfmt,rfir,fp)
use modmain
implicit none
! arguments
integer, intent(in) :: np
real(8), intent(in) :: vpl(3,np)
complex(8), intent(in) :: rfmt(npmtmax,natmtot),rfir(ngtot)
complex(8), intent(out) :: fp(np)
! local variables
integer ia,is,ias
integer nr,nri,ir0,ir
integer lmax,l,m,lm
integer ig,ifg,ip
integer i1,i2,i3,i,j
real(8) rmt2,r,tp(2),ya(4),t1,t2
real(8) v1(3),v2(3),v3(3),v4(3),v5(3)
complex(8) sum
! automatic arrays
real(8) rlm(lmmaxo)
! allocatable arrays
complex(8), allocatable :: rfmt1(:,:,:)
complex(8), allocatable :: zfft(:)
! unpack the muffin-tin function
allocate(rfmt1(lmmaxo,nrmtmax,natmtot))
do ias=1,natmtot
  is=idxis(ias)
  call rfmtpack(.false.,nrmt(is),nrmti(is),rfmt(:,ias),rfmt1(:,:,ias))
end do
! Fourier transform rfir to G-space
allocate(zfft(ngtot))
zfft(:)=rfir(:)
call zfftifc(3,ngridg,-1,zfft)

deallocate(rfmt1,zfft)
return

contains

real(8) function poly4(xa,ya,x)
implicit none
! arguments
real(8), intent(in) :: xa(4),ya(4),x
! local variables
real(8) x0,x1,x2,x3,y0,y1,y2,y3
real(8) c1,c2,c3,t0,t1,t2,t3,t4,t5,t6
! evaluate the polynomial coefficients
x0=xa(1)
x1=xa(2)-x0
x2=xa(3)-x0
x3=xa(4)-x0
y0=ya(1)
y1=ya(2)-y0
y2=ya(3)-y0
y3=ya(4)-y0
t0=1.d0/(x1*x2*x3*(x1-x2)*(x1-x3)*(x2-x3))
t1=x1*x2*y3
t2=x2*x3*y1
t3=x3*x1*y2
c3=t1*(x1-x2)+t2*(x2-x3)+t3*(x3-x1)
t6=x3**2
t5=x2**2
t4=x1**2
c2=t1*(t5-t4)+t2*(t6-t5)+t3*(t4-t6)
c1=t1*(x2*t4-x1*t5)+t2*(x3*t5-x2*t6)+t3*(x1*t6-x3*t4)
t1=x-x0
! evaluate the polynomial
poly4=y0+t0*t1*(c1+t1*(c2+c3*t1))
return
end function

end subroutine
