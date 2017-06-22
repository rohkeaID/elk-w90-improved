
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dolpalo(ias,ngp,ngpq,dapwalm,dapwalmq,ld,od)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: ias
integer, intent(in) :: ngp,ngpq
complex(8), intent(in) :: dapwalm(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: dapwalmq(ngkmax,apwordmax,lmmaxapw)
integer, intent(in) :: ld
complex(8), intent(inout) :: od(ld,*)
! local variables
integer is,ilo,io
integer l,m,lm,i,j
if (ias.ne.iasph) return
is=idxis(ias)
do ilo=1,nlorb(is)
  l=lorbl(ilo,is)
  do m=-l,l
    lm=idxlm(l,m)
    j=ngp+idxlo(lm,ilo,ias)
    do i=1,ngpq
      do io=1,apword(l,is)
        od(i,j)=od(i,j)+conjg(dapwalmq(i,io,lm))*oalo(io,ilo,ias)
      end do
    end do
    i=ngpq+idxlo(lm,ilo,ias)
    do j=1,ngp
      do io=1,apword(l,is)
        od(i,j)=od(i,j)+oalo(io,ilo,ias)*dapwalm(j,io,lm)
      end do
    end do
  end do
end do
return
end subroutine

