
! Copyright (C) 2015 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine doccupy
use modmain
use modphonon
implicit none
! local variables
integer ik,jk,ist
real(8) sum1,sum2
real(8) x,dx,t1,t2
! allocatable arrays
real(8), allocatable :: devalsv(:)
! external functions
real(8) sdelta
external sdelta
if (.not.tphiq0) return
allocate(devalsv(nstsv))
t1=1.d0/swidth
! compute the derivative of the Fermi energy
sum1=0.d0
sum2=0.d0
do ik=1,nkptnr
  jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
  call getdevalsv(ik,iqph,isph,iaph,ipph,devalsv)
  do ist=1,nstsv
    x=(efermi-evalsv(ist,jk))*t1
    t2=wkptnr*occmax*sdelta(stype,x)*t1
    sum1=sum1+t2
    sum2=sum2+t2*devalsv(ist)
  end do
end do
if (abs(sum1).gt.1.d-6) then
  defermi=sum2/sum1
else
  defermi=0.d0
end if
! determine the occupation number derivatives
do ik=1,nkptnr
  jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
  call getdevalsv(ik,iqph,isph,iaph,ipph,devalsv)
  do ist=1,nstsv
    x=(efermi-evalsv(ist,jk))*t1
    dx=(defermi-devalsv(ist))*t1
    doccsv(ist,ik)=occmax*sdelta(stype,x)*dx
  end do
end do
deallocate(devalsv)
return
end subroutine

