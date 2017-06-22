
! Copyright (C) 2014 D. Ernsting, S. Dugdale and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine emdplot
use modmain
use modpw
implicit none
! local variables
integer ik
real(8) t1
! allocatable arrays
real(4), allocatable :: emds(:,:)
t1=sum(abs(vkloff(:)))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(emdplot): use vkloff = 0 for the ground-state run")')
  write(*,*)
  stop
end if
! initialise universal variables
call init0
call init1
call init4
! read density and potentials from file
call readstate
! read Fermi energy from file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! get the occupancies from file
do ik=1,nkpt
  call getoccsv(filext,vkl(:,ik),occsv(:,ik))
end do
! read in the electron momentum density
allocate(emds(nhkmax,nkpt))
call reademd(emds)
! write the density plot to file
select case(task)
case(171)
  call emdplot1d(emds)
  write(*,*)
  write(*,'("Info(emdplot): 1D electron momentum density written to &
   &EMD1D.OUT")')
case(172)
  call emdplot2d(emds)
  write(*,*)
  write(*,'("Info(emdplot): 2D electron momentum density written to &
   &EMD2D.OUT")')
case(173)
  call emdplot3d(emds)
  write(*,*)
  write(*,'("Info(emdplot): 3D electron momentum density written to &
   &EMD3D.OUT")')
end select
deallocate(emds)
return
end subroutine

