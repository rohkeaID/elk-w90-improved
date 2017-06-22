
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dolpaa(ias,ngp,ngpq,apwalm,apwalmq,dapwalm,dapwalmq,ld,od)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: ias
integer, intent(in) :: ngp,ngpq
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(in) :: apwalmq(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(in) :: dapwalm(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: dapwalmq(ngkmax,apwordmax,lmmaxapw)
integer, intent(in) :: ld
complex(8), intent(inout) :: od(ld,*)
! local variables
integer is,l,m,lm,io
! automatic arrays
complex(8) x(ngpq)
if (ias.ne.iasph) return
is=idxis(ias)
lm=0
do l=0,lmaxmat
  do m=-l,l
    lm=lm+1
    do io=1,apword(l,is)
      x(1:ngpq)=conjg(apwalmq(1:ngpq,io,lm,ias))
      call zgerci(ngpq,ngp,zone,x,dapwalm(:,io,lm),ld,od)
      x(1:ngpq)=conjg(dapwalmq(1:ngpq,io,lm))
      call zgerci(ngpq,ngp,zone,x,apwalm(:,io,lm,ias),ld,od)
    end do
  end do
end do
return
end subroutine

