
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writestress
use modmain
use modmpi
implicit none
! local variables
integer i,j,k
! initialise universal variables
call init0
! start from the atomic densities
trdstate=.false.
! generate the stress tensor
call genstress
! write the stress tensor to file
if (mp_mpi) then
  open(50,file='STRESS.OUT',action='WRITE',form='FORMATTED')
  write(50,*)
  write(50,'("Lattice vector matrix, A, changed by")')
  write(50,*)
  write(50,'("     A --> A + e_i dt,")')
  write(50,*)
  write(50,'("where dt is an infinitesimal scalar and e_i is a strain tensor")')
  write(50,*)
  write(50,'("Stress is given by the derivative of the total energy dE/dt")')
  do k=1,nstrain
    write(50,*); write(50,*)
    write(50,'("Strain matrix : ",I1)') k
    do i=1,3
      write(50,'(3G18.10)') (strain(i,j,k),j=1,3)
    end do
    write(50,*)
    write(50,'("Stress matrix component : ",G18.10)') stress(k)
  end do
  close(50)
  write(*,*)
  write(*,'("Info(writestress):")')
  write(*,'(" Stress matrix written to STRESS.OUT")')
end if
return
end subroutine

