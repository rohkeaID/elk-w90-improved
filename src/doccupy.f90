
! Copyright (C) 2015 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine doccupy
use modmain
use modphonon
implicit none
! local variables
integer, parameter :: maxit=1000
integer ik,jk,ist,it
real(8) de0,de1,de
real(8) dchg,x,dx,t1
! external functions
real(8) sdelta
external sdelta
if (.not.tphq0) return
de0=1.d6
de1=-1.d6
do ik=1,nkptnr
  do ist=1,nstsv
    de=devalsv(ist,ik)
    if (de.lt.de0) de0=de
    if (de.gt.de1) de1=de
  end do
end do
t1=1.d0/swidth
do it=1,maxit
  defermi=0.5d0*(de0+de1)
  dchg=0.d0
  do ik=1,nkptnr
    jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
    do ist=1,nstsv
      x=(efermi-evalsv(ist,jk))*t1
      dx=(defermi-devalsv(ist,ik))*t1
      doccsv(ist,ik)=occmax*sdelta(stype,x)*dx
      dchg=dchg+wkptnr*doccsv(ist,ik)
    end do
  end do
  if (dchg.lt.0.d0) then
    de0=defermi
  else
    de1=defermi
  end if
  if ((de1-de0).lt.1.d-12) goto 10
end do
write(*,*)
write(*,'("Warning(doccupy): could not find Fermi energy derivative")')
10 continue
return
end subroutine

