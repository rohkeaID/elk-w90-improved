

! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genidxbdg
use modmain
use modscdft
implicit none
! local variables
integer ik,jk,ist,i
i=0
! loop over non-reduced k-points
do ik=1,nkptnr
! equivalent reduced k-point
  jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
  do ist=1,nstsv
    if (abs(evalsv(ist,jk)-efermi).lt.ewbdg) i=i+1
  end do
end do
nbdg=i
if (nbdg.eq.0) then
  write(*,*)
  write(*,'("Error(genidxbdg): no states found in BdG energy window")')
  write(*,*)
  stop
end if
! allocate global index array
if (allocated(idxbdg)) deallocate(idxbdg)
allocate(idxbdg(2,nbdg))
i=0
do ik=1,nkptnr
  jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
  do ist=1,nstsv
    if (abs(evalsv(ist,jk)-efermi).lt.ewbdg) then
      i=i+1
      idxbdg(1,i)=ik
      idxbdg(2,i)=ist
    end if
  end do
end do
return
end subroutine

