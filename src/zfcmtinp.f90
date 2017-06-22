
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

complex(8) function zfcmtinp(nr,nri,r,r2,zfmt1,zfmt2)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: r(nr),r2(nr)
complex(8), intent(in) :: zfmt1(lmmaxvr,nr),zfmt2(lmmaxvr,nr)
! local variables
integer ir
real(8) a,b,t1
complex(8) z1
! automatic arrays
real(8) fr1(nr),fr2(nr)
! external functions
real(8) fintgt
complex(8) zdotc
external fintgt,zdotc
do ir=1,nri
  z1=zdotc(lmmaxinr,zfmt1(:,ir),1,zfmt2(:,ir),1)*r2(ir)
  fr1(ir)=dble(z1)
  fr2(ir)=aimag(z1)
end do
t1=dble(lmmaxinr)/dble(lmmaxvr)
do ir=nri+1,nr
  z1=zdotc(lmmaxvr,zfmt1(:,ir),1,zfmt2(:,ir),1)*(t1*r2(ir))
  fr1(ir)=dble(z1)
  fr2(ir)=aimag(z1)
end do
! integrate
a=fintgt(-1,nr,r,fr1)
b=fintgt(-1,nr,r,fr2)
zfcmtinp=(fourpi/lmmaxinr)*cmplx(a,b,8)
return
end function

