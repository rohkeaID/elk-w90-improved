
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rdmeval
! !INTERFACE:
subroutine rdmeval
! !USES:
use modmain
use modrdm
! !DESCRIPTION:
!   RDMFT eigenvalues are determined by calculating the derivative of the total
!   energy with respect to the occupation number at half the maximum occupancy
!   ($n_{\rm max}/2$).
!
! !REVISION HISTORY:
!   Created 2009 (Sharma)
!EOP
!BOC
implicit none
! local variables
integer ik,ist
real(8) t1
! allocatable arrays
real(8), allocatable :: dedn(:,:)
allocate(dedn(nstsv,nkpt))
do ik=1,nkpt
  do ist=1,nstsv
    t1=occsv(ist,ik)
    occsv(ist,ik)=occmax/2.d0
    call rdmdedn(dedn)
    evalsv(ist,ik)=-dedn(ist,ik)
    occsv(ist,ik)=t1
  end do
  call putevalsv(filext,ik,evalsv(:,ik))
end do
deallocate(dedn)
return
end subroutine
!EOC

