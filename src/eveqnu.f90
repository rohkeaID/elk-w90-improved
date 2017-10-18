
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine eveqnu(ik0,ist)
use modmain
use modultra
implicit none
! arguments
integer, intent(in) :: ik0,ist
! local variables
integer i,j,ik
integer iv(3),jv(3)
integer iq
real(8) ca
real(8) t1,t2
real(8) s1,s2,s3
complex(8) z1
!***************


! allocatable arrays
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: evecu(:,:)
! external functions
real(8) dznrm2
complex(8) zdotc
external dznrm2,zdotc


allocate(evecu(nkpa,nkpa))

! coupling constant of the A-field (1/c)
ca=1.d0/solsc

! spin-polarised case
if (spinpol) then
! read in the second-variational eigenvector for kappa = 0
  allocate(evecsv(nstsv,nstsv))
  ik=nkpa*(ik0-1)+ikpa0
  call getevecsv(filext,ik,vkl(:,ik),evecsv)
  i=nstfv+1
  t1=dznrm2(nstfv,evecsv(:,ist),1)**2
  t2=dznrm2(nstfv,evecsv(i,ist),1)**2
  z1=zdotc(nstfv,evecsv(:,ist),1,evecsv(i,ist),1)
  s1=2.d0*dble(z1)
  s2=2.d0*aimag(z1)
  s3=t1-t2
  deallocate(evecsv)
end if
do j=1,nkpa
  jv(:)=ivkpa(:,j)
  do i=1,j
    iv(:)=ivkpa(:,i)-jv(:)
    iq=ivqiq(iv(1),iv(2),iv(3))
! scalar potential
    z1=vsuq(iq)
    if (spinpol) z1=z1+s1*bsuq(iq,1)+s2*bsuq(iq,2)+s3*bsuq(iq,3)
    evecu(i,j)=vsuq(iq)
  end do
end do
! add the second-variational eigenvalues along the diagonal
do i=1,nkpa
! index of the k+kappa-point
  ik=nkpa*(ik0-1)+i
  evecu(i,i)=evecu(i,i)+evalsv(ist,ik)
end do


deallocate(evecu)

return
end subroutine

