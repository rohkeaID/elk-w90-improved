
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

complex(8) function zfcmtinp(nr,nri,r,r2,zfmt1,zfmt2)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: r(nr),r2(nr)
complex(8), intent(in) :: zfmt1(*),zfmt2(*)
! local variables
integer ir,i0,i1
real(8) a,b,t1
complex(8) z1
! automatic arrays
real(8) fr1(nr),fr2(nr)
! external functions
real(8) fintgt
external fintgt
i1=0
do ir=1,nri
  i0=i1+1
  i1=i1+lmmaxi
  z1=dot_product(zfmt1(i0:i1),zfmt2(i0:i1))*r2(ir)
  fr1(ir)=dble(z1)
  fr2(ir)=aimag(z1)
end do
t1=dble(lmmaxi)/dble(lmmaxo)
do ir=nri+1,nr
  i0=i1+1
  i1=i1+lmmaxo
  z1=dot_product(zfmt1(i0:i1),zfmt2(i0:i1))*(t1*r2(ir))
  fr1(ir)=dble(z1)
  fr2(ir)=aimag(z1)
end do
! integrate
a=fintgt(-1,nr,r,fr1)
b=fintgt(-1,nr,r,fr2)
zfcmtinp=(fourpi/dble(lmmaxi))*cmplx(a,b,8)
return
end function

