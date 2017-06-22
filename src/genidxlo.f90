
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genidxlo
! !INTERFACE:
subroutine genidxlo
! !USES:
use modmain
! !DESCRIPTION:
!   Generates an index array which maps the local-orbitals in each atom to their
!   locations in the overlap or Hamiltonian matrices. Also finds the total
!   number of local-orbitals.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ias,i,ilo,l,m,lm
! allocate global local-orbital index
if (allocated(idxlo)) deallocate(idxlo)
allocate(idxlo(lolmmax,nlomax,natmtot))
i=0
do ias=1,natmtot
  is=idxis(ias)
  do ilo=1,nlorb(is)
    l=lorbl(ilo,is)
    do m=-l,l
      i=i+1
      lm=idxlm(l,m)
      idxlo(lm,ilo,ias)=i
    end do
  end do
end do
nlotot=i
return
end subroutine
!EOC
