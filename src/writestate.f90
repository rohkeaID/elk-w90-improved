
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writestate
! !INTERFACE:
subroutine writestate
! !USES:
use modmain
use moddftu
! !DESCRIPTION:
!   Writes the charge density, potentials and other relevant variables to the
!   file {\tt STATE.OUT}. Note to developers: changes to the way the variables
!   are written should be mirrored in {\tt readstate}.
!
! !REVISION HISTORY:
!   Created May 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is
open(50,file='STATE'//trim(filext),action='WRITE',form='UNFORMATTED')
write(50) version
write(50) spinpol
write(50) nspecies
write(50) lmmaxvr
write(50) nrmtmax
write(50) nrcmtmax
do is=1,nspecies
  write(50) natoms(is)
  write(50) nrmt(is)
  write(50) rsp(1:nrmt(is),is)
  write(50) nrcmt(is)
  write(50) rcmt(1:nrcmt(is),is)
end do
write(50) ngridg
write(50) ngvec
write(50) ndmag
write(50) nspinor
write(50) fsmtype
write(50) ftmtype
write(50) dftu
write(50) lmmaxdm
! write the density
write(50) rhomt,rhoir
! write the Coulomb potential
write(50) vclmt,vclir
! write the exchange-correlation potential
write(50) vxcmt,vxcir
! write the Kohn-Sham effective potential
write(50) vsmt,vsir
if (spinpol) then
! write the magnetisation, exchange-correlation and effective magnetic fields
  write(50) magmt,magir
  write(50) bxcmt,bxcir
  write(50) bsmt,bsir
! write fixed spin moment magnetic fields
  if (fsmtype.ne.0) then
    write(50) bfsmc
    write(50) bfsmcmt
  end if
end if
! write the potential matrix in each muffin-tin
if ((dftu.ne.0).or.(ftmtype.ne.0)) then
  write(50) vmatmt
end if
! write the fixed tensor moment potential matrix
if (ftmtype.ne.0) then
  write(50) vmftm
end if
close(50)
return
end subroutine
!EOC

