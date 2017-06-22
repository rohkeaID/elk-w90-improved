
! Copyright (C) 2014 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeafpdt
use modmain
use modtddft
implicit none
! local variables
integer its,i
real(8) ed,t1
! conversion factor of power density to J/cm^2
real(8), parameter :: ced=ha_si/(100.d0*br_si)**2
! allocatable arrays
real(8), allocatable :: f(:),g(:),pd(:)
! allocate local arrays
allocate(f(ntimes),g(ntimes),pd(ntimes))
! compute the power density at each time step
pd(:)=0.d0
do i=1,3
  f(:)=afieldt(i,:)
  call fderiv(1,ntimes,times,f,g)
  pd(:)=pd(:)+g(:)**2
end do
t1=fourpi/solsc
pd(:)=t1*pd(:)
! write the power density to file
open(50,file='AFPDT.OUT',action='WRITE',form='FORMATTED')
do its=1,ntimes
  write(50,'(2G18.10)') times(its),pd(its)
end do
close(50)
! integrate power density to find the total energy density
call fderiv(-1,ntimes,times,pd,g)
ed=g(ntimes)
open(50,file='AFTED.OUT',action='WRITE',form='FORMATTED')
write(50,*)
write(50,'("Total energy density : ",G18.10)') ed
write(50,'(" in J/cm^2           : ",G18.10)') ed*ced
close(50)
write(*,*)
write(*,'("Info(writeafpdt):")')
write(*,'(" Power density of A-field written to AFPDT.OUT")')
write(*,'(" Total energy density of A-field written to AFTED.OUT")')
deallocate(f,g,pd)
return
end subroutine

