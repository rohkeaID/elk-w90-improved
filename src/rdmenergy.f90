
! Copyright (C) 2005-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rdmenergy
! !INTERFACE:
subroutine rdmenergy
! !USES:
use modmain
use modrdm
use modmpi
use modtest
! !DESCRIPTION:
!   Calculates RDMFT total energy (free energy for finite temperatures).
!
! !REVISION HISTORY:
!   Created 2008 (Sharma)
!   Updated for free energy 2009 (Baldsiefen)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias
integer nr,ik,ist
real(8) vn,t1
complex(8) z1
! allocatable arrays
real(8), allocatable :: rfmt(:,:)
complex(8), allocatable :: evecsv(:,:)
! external functions
real(8) rfmtinp
complex(8) zdotc
external rfmtinp,zdotc
allocate(rfmt(lmmaxvr,nrmtmax))
! Coulomb energy from core states
engyvcl=0.d0
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  rfmt(:,:)=0.d0
  if (spincore) then
    rfmt(1,1:nr)=(rhocr(1:nr,ias,1)+rhocr(1:nr,ias,2))/y00
  else
    rfmt(1,1:nr)=rhocr(1:nr,ias,1)/y00
  end if
  engyvcl=engyvcl+rfmtinp(nr,nrmtinr(is),1,rsp(:,is),r2sp(:,is),rfmt, &
   vclmt(:,:,ias))
end do
deallocate(rfmt)
engykn=engykncr
allocate(evecsv(nstsv,nstsv))
do ik=1,nkpt
  call getevecsv(filext,vkl(:,ik),evecsv)
  do ist=1,nstsv
    t1=wkpt(ik)*occsv(ist,ik)
! Coulomb energy from valence states
    engyvcl=engyvcl+t1*dble(vclmat(ist,ist,ik))
! kinetic energy from valence states
    z1=zdotc(nstsv,evecsv(:,ist),1,dkdc(:,ist,ik),1)
    engykn=engykn+t1*dble(z1)
  end do
end do
deallocate(evecsv)
! Madelung term
engymad=0.d0
do is=1,nspecies
! compute the bare nucleus potential at the origin
  call potnucl(ptnucl,1,rsp(:,is),spzn(is),vn)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    engymad=engymad+0.5d0*spzn(is)*(vclmt(1,1,ias)*y00-vn)
  end do
end do
! exchange-correlation energy
call rdmengyxc
! total energy
engytot=0.5d0*engyvcl+engymad+engykn+engyx
if (rdmtemp.gt.0.d0) then
  call rdmentropy
  engytot=engytot-rdmtemp*rdmentrpy
end if
! write total energy to test file
call writetest(300,'RDMFT total energy',tol=1.d-6,rv=engytot)
return
end subroutine
!EOC
