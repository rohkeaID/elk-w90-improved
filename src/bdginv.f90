
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine bdginv(achi0,adelta)
use modmain
use modscdft
implicit none
! arguments
complex(8), intent(in) :: achi0(nbdg,nbdg)
complex(8), intent(out) :: adelta(nbdg,nbdg)
! local variables
integer, parameter :: maxit=1000
integer it
real(8) d,dp
! allocatable arrays
complex(8), allocatable :: achi(:,:)
complex(8), allocatable :: evecbdg(:,:)
! external functions
real(8) dznrm2
allocate(achi(nbdg,nbdg))
allocate(evecbdg(nmbdg,nmbdg))
! zero the anomalous potential
adelta(:,:)=0.d0
dp=0.d0
do it=1,maxit
! set up the BdG Hamiltonian
  call hmlbdg(adelta,evecbdg)
! solve the BdG eigenvalue equations
  call eveqnbdg(evecbdg)
! generate the anomalous density Chi
  call genachi(evecbdg,achi)
! add the residual to the potential
  adelta(:,:)=adelta(:,:)+taubdg*(achi(:,:)-achi0(:,:))
! compute the sum of diagonal elements squared
  d=dznrm2(nbdg,adelta,nbdg)
  if (it.ge.2) then
    d=sqrt(abs(d)/dble(nbdg))
    if (abs(d-dp).lt.epspot) return
  end if
  dp=d
end do
write(*,*)
write(*,'("Warning(bdginv): BdG equation inverter failed to converge")')
deallocate(achi,evecbdg)
return
end subroutine

