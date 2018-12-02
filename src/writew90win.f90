
! Copyright (C) 2015 Manh Duc Le, 2017-18 Arsenii Gerasimov, Yaroslav Kvashnin
! and Lars Nordstrom.
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
!   Writes out a template {\tt seedname.win} file with the Wannier90 input
!   parameters. Uses {\tt wannier} and {\tt wannierExtra} blocks.
!
! !REVISION HISTORY:
!   Created January 2015 (Manh Duc Le)
!   Modified August 2018 (Arsenii Gerasimov)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ik,i
! automatic arrays
real(8)        v2(3)
character(256) filename
!-------------------------------------------------------------------------------

! Check that user has defined bands and projections
if( wann_projlines .eq. -1 ) then
  write(*,*)
  write(*,'("Error(writew90win): No projections specified - &
            &please input using the wann_projections block.")')
  write(*,*)
  stop
end if
if( wann_nband .eq. -1 ) then
  write(*,*)
  write(*,'("Error(writew90win): No bands specified - &
            &please input using the wann_bands block.")')
  write(*,*)
  stop
end if

write(*,*)
write(*,'("  Info(Wannier): Number of bands to be used for &
                                            &wannierization :",i3)'), wann_nband
write(*,'(17x,"Indices of the corresponding bands/eigenvalues :",99i4)') &
                                                        wann_bands(1:wann_nband)
write(*,*)

! Initialise universal variables
call init0

! Code to generate the k-point mesh for w90 (from init1.f90)
! Set up the k-point box (based on k-points set, specified in Wannier-block)
ngridk   = wann_ngridk
reducek0 = reducek ! If reducek=1 was used in ground state calculations,
                   ! need to regenerate the eigenvectors set for the full BZ.
reducek  = 0
call init1

! Check that the number of bands is greater than or equal to the number projections
if( wann_nband .lt. wann_nwf ) then
  write(*,*)
  write(*,'("Error(writew90win): Number of bands less than number of wannier &
                      &functions. Please check wann_nwf/wann_projections and &
                                                                &wann_bands")')
  write(*,*)
  stop
end if

filename = trim(wann_seedname)//'.win'
open(500,file=filename,action='WRITE',form='FORMATTED')
! Write the number of Wannier functions and number of bands to calculate
write(500,'(" length_unit=Bohr")')
write(500,'(" num_wann =  ",I8)')  wann_nwf
write(500,'(" num_bands = ",I8)')  wann_nband
write(500,'(" num_iter = ",I8)')   wann_numiter

! Write spinors, if necessary
if ( nspinor .eq. 2 ) then
  write(500,*)
  write(500,'(" spinors = true")')
  write(500,'(" spn_formatted = true")')
end if

! Write unit cell
write(500,*)
write(500,'(" begin unit_cell_cart")')
write(500,'(" bohr")')
write(500,'(3G18.10)') avec(:,1)
write(500,'(3G18.10)') avec(:,2)
write(500,'(3G18.10)') avec(:,3)
write(500,'(" end unit_cell_cart")')

! Writes atomic positions
write(500,*)
write(500,'(" begin atoms_frac")')
do is = 1,nspecies
  do ia = 1,natoms(is)
    ! Code taken from writegeom.f90
    if ( molecule ) then
      ! Map lattice coordinates to [-0.5,0.5)
      v2(:) = atposl(:,ia,is)
      do i = 1,3
        if ( v2(i) .gt. 0.5d0 ) v2(i) = v2(i) - 1.d0
      end do
    else
      ! Otherwise write lattice coordinates
      v2(:) = atposl(:,ia,is)
    end if
    write(500,'(A5,3G18.10)') trim(spsymb(is)),v2(:)
  end do
end do
write(500,'(" end atoms_frac")')

! Write the projections block
write(500,*)
write(500,'(" begin projections")')
do i = 1,wann_projlines
  write(500,'(A)') trim(wann_projstr(i))
end do
write(500,'(" end projections")')

! Write the list of k-points
write(500,*)
write(500,'(" mp_grid = ",3I8)') ngridk
write(500,*)
write(500,'(" begin kpoints")')
do ik = 1,nkpt
  write(500,'(3G18.10)') vkl(:,ik)
end do
write(500,'(" end kpoints")')
if( wann_inputlines .gt. 0 ) then
  do i = 1,wann_inputlines
    write(500,'(A)') trim(wann_input(i))
  end do
end if

close(500)

reducek = reducek0

write(*,*)
write(*,*) " Info(Wannier): <seedname>.win file has been created"
write(*,*)

end subroutine
!EOC
