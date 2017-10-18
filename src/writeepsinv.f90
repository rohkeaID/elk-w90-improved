
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeepsinv
use modmain
use modmpi
implicit none
! local variables
integer ik
! initialise global variables
call init0
call init1
call init2
call init3
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
! get the eigenvalues and occupancies from file
do ik=1,nkpt
  call getevalsv(filext,ik,vkl(:,ik),evalsv(:,ik))
  call getoccsv(filext,ik,vkl(:,ik),occsv(:,ik))
end do
! generate the inverse dielectric function and write to file
call epsinv
if (mp_mpi) then
  write(*,*)
  write(*,'("Info(writeepsinv):")')
  write(*,'(" inverse RPA dielectric function, eps^(-1)(G,G'',q,w), written to &
   &EPSINV.OUT")')
end if
return
end subroutine

