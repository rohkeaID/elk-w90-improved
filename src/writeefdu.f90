
! Copyright (C) 2008  F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writeefdu
! !INTERFACE:
subroutine writeefdu
! !USES:
use modmain
use moddftu
! !DESCRIPTION:
!   Writes to file the linearisation energies for all radial functions used to
!   calculate the Slater integrals through a Yukawa potential.
!
! !REVISION HISTORY:
!   Created July 2008 (Francesco Cricchio)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias,i,l
open(50,file='ENGYFDU'//trim(filext),form='FORMATTED')
do i=1,ndftu
  is=idftu(1,i)
  l=idftu(2,i)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(50,*)
    write(50,'("Species : ",I4," (",A,"), atom : ",I4)') is,trim(spsymb(is)),ia
    write(50,'(" Radial functions used for Slater parameters :")')
    write(50,'("  l = ",I2," : ",G18.10)') l,fdue(l,ias)
  end do
end do
close(50)
return
end subroutine
!EOC

