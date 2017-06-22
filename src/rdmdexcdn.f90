
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rdmdexcdn
! !INTERFACE:
subroutine rdmdexcdn(dedn)
! !USES:
use modmain
use modrdm
! !INPUT/OUTPUT PARAMETERS:
!   dedn : energy derivative (inout,real(nstsv,nkpt))
! !DESCRIPTION:
!   Calculates the derivative of the exchange-correlation energy w.r.t.
!   occupation numbers and adds the result to the total.
!
! !REVISION HISTORY:
!   Created 2008 (Sharma)
!EOP
!BOC
implicit none
! arguments
real(8), intent(inout) :: dedn(nstsv,nkpt)
! local variables
integer ik1,ik2,jk,iv(3)
integer ist1,ist2
! parameter for calculating the functional derivatives
real(8), parameter :: eps=1.d-12
real(8) t1,t2,t3,t4
! allocatable arays
real(8), allocatable :: vclijji(:,:,:)
if (rdmxctype.eq.0) return
! calculate the prefactor
if (rdmxctype.eq.1) then
  t1=1.d0/occmax
! power functional
else if (rdmxctype.eq.2) then
  if (spinpol) then
    t1=rdmalpha
  else
    t1=2.d0*rdmalpha*(0.25d0)**rdmalpha
  end if
else
  write(*,*)
  write(*,'("Error(rdmdexcdn): rdmxctype not defined : ",I8)') rdmxctype
  write(*,*)
  stop
end if
allocate(vclijji(nstsv,nstsv,nkpt))
! start loop over non-reduced k-points
do ik1=1,nkptnr
! get the Coulomb matrix elements
  call getvclijji(ik1,vclijji)
! find the equivalent reduced k-point
  iv(:)=ivk(:,ik1)
  jk=ivkik(iv(1),iv(2),iv(3))
!  loop over reduced k-points
  do ik2=1,nkpt
    do ist1=1,nstsv
      do ist2=1,nstsv
! Hartree-Fock functional
        if (rdmxctype.eq.1) then
          t2=t1*occsv(ist1,jk)
! power functional
        else if (rdmxctype.eq.2) then
          t3=sum(abs(vkl(:,ik2)-vkl(:,ik1)))
          if ((ist2.eq.ist1).and.(t3.lt.epslat)) then
            t2=(1.d0/occmax)*occsv(ist1,jk)
          else
            t3=max(occsv(ist2,ik2),eps)
            t4=max(occsv(ist1,jk),eps)
            t2=t1*(t4**rdmalpha)/(t3**(1.d0-rdmalpha))
          end if
        end if
        dedn(ist2,ik2)=dedn(ist2,ik2)+t2*vclijji(ist2,ist1,ik2)
      end do
    end do
  end do
end do
deallocate(vclijji)
return
end subroutine
!EOC

