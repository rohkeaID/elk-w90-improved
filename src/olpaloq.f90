
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine olpaloq(ias,ngp,ngpq,apwalm,apwalmq,ld,oq)
use modmain
implicit none
! arguments
integer, intent(in) :: ias
integer, intent(in) :: ngp,ngpq
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(in) :: apwalmq(ngkmax,apwordmax,lmmaxapw,natmtot)
integer, intent(in) :: ld
complex(8), intent(inout) :: oq(ld,*)
! local variables
integer is,ilo,io
integer l,m,lm,i,j
is=idxis(ias)
do ilo=1,nlorb(is)
  l=lorbl(ilo,is)
  do m=-l,l
    lm=idxlm(l,m)
    j=ngp+idxlo(lm,ilo,ias)
    do i=1,ngpq
      do io=1,apword(l,is)
        oq(i,j)=oq(i,j)+conjg(apwalmq(i,io,lm,ias))*oalo(io,ilo,ias)
      end do
    end do
    i=ngpq+idxlo(lm,ilo,ias)
    do j=1,ngp
      do io=1,apword(l,is)
        oq(i,j)=oq(i,j)+oalo(io,ilo,ias)*apwalm(j,io,lm,ias)
      end do
    end do
  end do
end do
return
end subroutine





