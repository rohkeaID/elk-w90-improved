
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dmatch(ias,ip,ngp,vgpc,apwalm,dapwalm)
use modmain
implicit none
! arguments
integer, intent(in) :: ias,ip,ngp
real(8), intent(in) :: vgpc(3,ngkmax)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(out) :: dapwalm(ngkmax,apwordmax,lmmaxapw)
! local variables
integer is,l,m,lm,io,igp
complex(8) z1
! take derivative with respect to atomic displacement
is=idxis(ias)
lm=0
do l=0,lmaxapw
  do m=-l,l
    lm=lm+1
    do io=1,apword(l,is)
      do igp=1,ngp
        z1=apwalm(igp,io,lm,ias)
        dapwalm(igp,io,lm)=vgpc(ip,igp)*cmplx(-aimag(z1),dble(z1),8)
      end do
    end do
  end do
end do
return
end subroutine

