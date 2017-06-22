
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: eveqnbdg
! !INTERFACE:
subroutine eveqnbdg(evecbdg)
! !USES:
use modmain
use modscdft
! !INPUT/OUTPUT PARAMETERS:
!   evecbdg : the BdG Hamiltonian matrix on input, the eigenvectors on output
!             (inout,complex(nmbdg,nmbdg))
! !DESCRIPTION:
!   Finds the eigenvalues and eigenvectors of the Bogoliubov-de Gennes
!   Hamiltonian. See {\tt hmlbdg}.
!
! !REVISION HISTORY:
!   Created January 2012 (JKD)
!EOP
!BOC
implicit none
! arguments
complex(8), intent(inout) :: evecbdg(nmbdg,nmbdg)
! local variables
integer lwork,info
! allocatable arrays
real(8), allocatable :: rwork(:)
complex(8), allocatable :: work(:)
allocate(rwork(3*nmbdg))
lwork=2*nmbdg
allocate(work(lwork))
call zheev('V','U',nmbdg,evecbdg,nmbdg,evalbdg,work,lwork,rwork,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(eveqnbdg): diagonalisation of the BdG Hamiltonian failed")')
  write(*,'(" ZHEEV returned INFO = ",I8)') info
  write(*,*)
  stop
end if
deallocate(rwork,work)
return
end subroutine
!EOC

