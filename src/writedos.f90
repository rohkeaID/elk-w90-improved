
! Copyright (C) 2015 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writedos
use modmain
implicit none
! local variables
integer ik
! initialise universal variables
call init0
call init1
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
  if (dosocc) call getoccsv(filext,ik,vkl(:,ik),occsv(:,ik))
end do
! generate the partial and total DOS and write to file
call dos('.OUT',dosocc,occsv)
write(*,*)
write(*,'("Info(writedos):")')
write(*,'(" Total density of states written to TDOS.OUT")')
write(*,*)
write(*,'(" Partial density of states written to PDOS_Sss_Aaaaa.OUT")')
write(*,'(" for all species and atoms")')
if (dosmsum) then
  write(*,'(" PDOS summed over m")')
end if
if (dosssum) then
  write(*,'(" PDOS summed over spin")')
end if
write(*,*)
write(*,'(" Spin-quantisation axis : ",3G18.10)') sqados(:)
if (lmirep) then
  write(*,*)
  write(*,'(" Eigenvalues of a random matrix in the (l,m) basis symmetrised")')
  write(*,'(" with the site symmetries written to ELMIREP.OUT for all")')
  write(*,'(" species and atoms. Degenerate eigenvalues correspond to")')
  write(*,'(" irreducible representations of each site symmetry group")')
end if
write(*,*)
write(*,'(" Interstitial density of states written to IDOS.OUT")')
write(*,*)
write(*,'(" Fermi energy is at zero in plots")')
write(*,*)
write(*,'(" DOS units are states/Hartree/unit cell")')
return
end subroutine

