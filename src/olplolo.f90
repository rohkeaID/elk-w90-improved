
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine olplolo(ias,ngp,ld,o)
use modmain
implicit none
! arguments
integer, intent(in) :: ias
integer, intent(in) :: ngp
integer, intent(in) :: ld
complex(8), intent(inout) :: o(ld,*)
! local variables
integer is,ilo,jlo
integer l,m,lm,i,j
is=idxis(ias)
do ilo=1,nlorb(is)
  l=lorbl(ilo,is)
  do jlo=1,nlorb(is)
    if (lorbl(jlo,is).eq.l) then
      do m=-l,l
        lm=idxlm(l,m)
        i=ngp+idxlo(lm,ilo,ias)
        j=ngp+idxlo(lm,jlo,ias)
        if (i.le.j) then
          o(i,j)=o(i,j)+ololo(ilo,jlo,ias)
        end if
      end do
    end if
  end do
end do
return
end subroutine

