
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rhonorm
! !INTERFACE:
subroutine rhonorm
! !USES:
use modmain
! !DESCRIPTION:
!   Loss of precision of the calculated total charge can result because the
!   muffin-tin density is computed on a set of $(\theta,\phi)$ points and then
!   transformed to a spherical harmonic representation. This routine adds a
!   constant to the density so that the total charge is correct. If the error in
!   total charge exceeds a certain tolerance then a warning is issued.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!   Changed from rescaling to adding, September 2006 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias,ir
real(8) t1,t2
if (.not.trhonorm) return
! check error in total charge
t1=chgcalc/chgtot-1.d0
if (abs(t1).gt.epschg) then
  write(*,*)
  write(*,'("Warning(rhonorm): total charge density incorrect for s.c. &
   &loop ",I5)') iscl
  write(*,'(" Calculated : ",G18.10)') chgcalc
  write(*,'(" Required   : ",G18.10)') chgtot
end if
! error in average density
t1=(chgtot-chgcalc)/omega
! add the constant difference to the density
t2=t1/y00
do ias=1,natmtot
  is=idxis(ias)
  do ir=1,nrmt(is)
    rhomt(1,ir,ias)=rhomt(1,ir,ias)+t2
  end do
end do
rhoir(:)=rhoir(:)+t1
! add the difference to the charges
do is=1,nspecies
  t2=t1*(4.d0*pi/3.d0)*rmt(is)**3
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    chgmt(ias)=chgmt(ias)+t2
    chgmttot=chgmttot+t2
  end do
end do
chgir=chgtot-chgmttot
chgcalc=chgtot
return
end subroutine
!EOC
