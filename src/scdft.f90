
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine scdft
use modmain
use modscdft
implicit none
! local variables
integer ik
! allocatable arrays
complex(8), allocatable :: achi(:,:)
! initialise global variables
call init0
call init1
call init2
! read Fermi energy from file
call readfermi
! get the eigenvalues and occupancies from file
do ik=1,nkpt
  call getevalsv(filext,vkl(:,ik),evalsv(:,ik))
  call getoccsv(filext,vkl(:,ik),occsv(:,ik))
end do
! generate the index to the BdG states
call genidxbdg
! initialise the anomalous density
allocate(achi(nbdg,nbdg))
call achiinit(achi)
! begin self-consistent loop
do iscl=1,maxscl
! invert the BdG equations
  call bdginv(achi)
end do
deallocate(achi)
return
end subroutine

