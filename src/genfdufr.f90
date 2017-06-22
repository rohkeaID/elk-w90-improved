
! Copyright (C) 2008  F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genfdufr
! !INTERFACE:
subroutine genfdufr
! !USES:
use modmain
use moddftu
! !DESCRIPTION:
!   Generates the radial functions used to calculate the Slater integrals
!   through a Yukawa potential.
!
! !REVISION HISTORY:
!   Created April 2008 from genapwfr (Francesco Cricchio)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias
integer nr,ir,nn,i,l
real(8) t1
! automatic arrays
real(8) vr(nrmtmax),fr(nrmtmax)
real(8) p0(nrmtmax),p1(nrmtmax),q0(nrmtmax),q1(nrmtmax)
! external functions
real(8) fintgt
external fintgt
do i=1,ndftu
  is=idftu(1,i)
  l=idftu(2,i)
  nr=nrmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    vr(1:nr)=vsmt(1,1:nr,ias)*y00
! integrate the radial Schrodinger equation
    call rschrodint(solsc,l,fdue(l,ias),nr,rsp(:,is),vr,nn,p0,p1,q0,q1)
! normalise radial functions
    fr(1:nr)=p0(1:nr)**2
    t1=fintgt(-1,nr,rsp(:,is),fr)
    if (t1.lt.1.d-20) then
      write(*,*)
      write(*,'("Error(genfdufr): degenerate APW radial functions")')
      write(*,'(" for species ",I4)') is
      write(*,'(" atom ",I4)') ia
      write(*,'(" and angular momentum ",I4)') l
      write(*,*)
      stop
    end if
    t1=1.d0/sqrt(abs(t1))
    p0(1:nr)=t1*p0(1:nr)
! divide by r and store in global array
    do ir=1,nr
      fdufr(ir,l,ias)=p0(ir)/rsp(ir,is)
    end do
  end do
end do
return
end subroutine
!EOC

