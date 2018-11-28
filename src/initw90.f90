
! Copyright (C) 2017-18 Arsenii Gerasimov, Yaroslav Kvashnin and Lars Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: initw90
! !INTERFACE:
subroutine initw90
! !USES:
use modmain
use modw90
! !DESCRIPTION:
!
!
! !REVISION HISTORY:
!   Created November 2018 (Arsenii Gerasimov)
!EOP
!BOC
implicit none
! local variables
logical exists
integer redkfil,nstsv_,recl
integer ikp
real(8) t1
! allocatable arrays
real(8),    allocatable :: evalsv_(:)
! automatic arrays
real(8) vkl_(3)
!-------------------------------------------------------------------------------
reducek0 = reducek ! If reducek = 1 was used in ground state calculations,
                   ! need to regenerate the eigenvectors set for the full BZ.
reducek  = 0

!
call setupw90lib

! Read density and potentials from file
call readstate
! Fourier transform Kohn-Sham potential to G-space
call genvsig
! Read Fermi energy from a file
call readfermi
! Find the new linearisation energies
call linengy
! Generate the APW radial functions
call genapwfr
! Generate the local-orbital radial functions
call genlofr

! Check that EVECSV.OUT has all necessary k-points
allocate(evalsv_(nstsv))
redkfil = 0
inquire(iolength=recl) vkl_,nstsv_,evalsv_
do ikp = 1,nkpt
  exists = .false.
  t1 = 9.d99
  inquire(file='EVALSV'//trim(filext),exist=exists)
  if(exists) then
    open(70,file='EVALSV'//trim(filext),action='READ',form='UNFORMATTED', &
        access='DIRECT',recl=recl,err=101)
    read(70,rec=ikp,err=101) vkl_,nstsv_,evalsv_
    close(70)
    t1 = abs( vkl(1,ikp) - vkl_(1) )&
       + abs( vkl(2,ikp) - vkl_(2) )&
       + abs( vkl(3,ikp) - vkl_(3) )
  end if
101 continue
  if ( .not.exists .or. t1.gt.epslat .or. nstsv.ne.nstsv_ ) then
    redkfil = 1
    exit
  end if
end do
! If kpoint not found in saved eigen-values/vectors, then need to recompute
! EVEC*OUT.
if ( redkfil .ne. 0 ) then
  ! Compute the overlap radial integrals
  call olprad
  ! Compute the Hamiltonian radial integrals
  call hmlrad
  ! Generate the spin-orbit coupling radial functions
  call gensocfr
  ! Generate the first- and second-variational eigenvectors and eigenvalues
  call genevfsv
end if

deallocate(evalsv_)

end subroutine initw90
!EOC
