! Copyright (C) 2015 Manh Duc Le
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writew90win
! !INTERFACE:
subroutine writew90win
! !USES:
use modmain
use modw90
! !DESCRIPTION:
!   Writes a template WIN input file for Wannier90 with the crystal structure
!
! !REVISION HISTORY:
!   Created January 2015 (Manh Duc Le)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ik,i
integer wnkptnr,wnkpt
! allocatable arrays
integer, allocatable :: wikmap(:,:,:),wikmapnr(:,:,:),wivk(:,:)
real(8), allocatable :: wvkl(:,:),wvkc(:,:),wwkpt(:)
! automatic arrays
real(8) v2(3)
real(8) wkptboxl(3,4)
character(256) filename
! fixed values

! Checks that user has defined bands and projections
if(wann_projlines.eq.-1) then
  write(*,*)
  write(*,'("Error(writew90win): No projections specified - &
        &please input using the wann_projections block.")')
  write(*,*)
  stop
end if
if(wann_nband.eq.-1) then
  write(*,*)
  write(*,'("Error(writew90win): No bands specified - &
        &please input using the wann_bands block.")')
  write(*,*)
  stop
end if

! initialise universal variables
call init0

! Code to generate the k-point mesh for w90 (from init1.f90)
! set up the default k-point box
wkptboxl=0.d0
wkptboxl(1,2)=1.d0
wkptboxl(2,3)=1.d0
wkptboxl(3,4)=1.d0
! allocate the k-point set arrays
allocate(wikmap(0:ngridk(1)-1,0:ngridk(2)-1,0:ngridk(3)-1))
allocate(wikmapnr(0:ngridk(1)-1,0:ngridk(2)-1,0:ngridk(3)-1))
wnkptnr=ngridk(1)*ngridk(2)*ngridk(3)
allocate(wivk(3,wnkptnr))
allocate(wvkl(3,wnkptnr))
allocate(wvkc(3,wnkptnr))
allocate(wwkpt(wnkptnr))
! generate the k-point set
call genppts(.false.,1,symlat,ngridk,wnkptnr,epslat,bvec,wkptboxl,wnkpt, &
    wikmap,wikmapnr,wivk,wvkl,wvkc,wwkpt,wkptnr)

! reads the projections and checks that the number of wannier functions
! given is equal to the number of projections
call getw90proj
! checks that the number of bands is greater than or equal to the number projections
if(wann_nband.lt.wann_nwf) then
  write(*,*)
  write(*,'("Error(writew90win): Number of bands less than number of&
        &wannier functions. Please check wann_nwf/wann_projections and wann_bands")')
  write(*,*)
  stop
end if

filename = trim(wann_seedname)//'.win'
open(50,file=filename,action='WRITE',form='FORMATTED')
! writes the number of Wannier functions and number of bands to calculate
write(50,'(" length_unit=Bohr")')
write(50,'(" num_wann =  ",I8)') wann_nwf
write(50,'(" num_bands = ",I8)') wann_nband
write(50,'(" num_iter = ",I8)') wann_numiter

! writes spinors if necessary
if (nspinor.eq.2) then
write(50,*)
write(50,'(" spinors = true")')
write(50,'(" spn_formatted = true")')
end if

! writes the unit cell
write(50,*)
write(50,'(" begin unit_cell_cart")')
write(50,'(" bohr")')
write(50,'(3G18.10)') avec(:,1)
write(50,'(3G18.10)') avec(:,2)
write(50,'(3G18.10)') avec(:,3)
write(50,'(" end unit_cell_cart")')
! writes the atomic positions
write(50,*)
write(50,'(" begin atoms_frac")')
do is=1,nspecies
  do ia=1,natoms(is)
! code taken from writegeom.f90
    if (molecule) then
! map lattice coordinates to [-0.5,0.5)
      v2(:)=atposl(:,ia,is)
      do i=1,3
        if (v2(i).gt.0.5d0) v2(i)=v2(i)-1.d0
      end do
    else
! otherwise write lattice coordinates
      v2(:)=atposl(:,ia,is)
    end if
    write(50,'(A5,3G18.10)') trim(spsymb(is)),v2(:)
  end do
end do
write(50,'(" end atoms_frac")')
! writes the projections block
write(50,*)
write(50,'(" begin projections")')
do i=1,wann_projlines
  write(50,'(A)') trim(wann_projstr(i))
end do
write(50,'(" end projections")')
! writes the list of k-points
write(50,*)
write(50,'(" mp_grid = ",3I8)') ngridk
write(50,*)
write(50,'(" begin kpoints")')
do ik=1,wnkpt
  write(50,'(3G18.10)') wvkl(:,ik)
end do
write(50,'(" end kpoints")')
if(wann_inputlines.gt.0) then
  do i=1,wann_inputlines
    write(50,'(A)') trim(wann_input(i))
  end do
end if
close(50)
!end do

deallocate(wikmap)
deallocate(wikmapnr)
deallocate(wivk)
deallocate(wvkl)
deallocate(wvkc)
deallocate(wwkpt)

end subroutine
!EOC
