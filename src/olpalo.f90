
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine olpalo(ias,ngp,apwalm,ld,o)
use modmain
implicit none
! arguments
integer, intent(in) :: ias
integer, intent(in) :: ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
integer, intent(in) :: ld
complex(8), intent(inout) :: o(*)
! local variables
integer is,ilo,io
integer l,m,lm,i,j,k
is=idxis(ias)
do ilo=1,nlorb(is)
  l=lorbl(ilo,is)
  do m=-l,l
    lm=idxlm(l,m)
    j=ngp+idxlo(lm,ilo,ias)
    k=(j-1)*ld
    do i=1,ngp
      k=k+1
      do io=1,apword(l,is)
        o(k)=o(k)+conjg(apwalm(i,io,lm,ias))*oalo(io,ilo,ias)
      end do
    end do
  end do
end do
return
end subroutine

