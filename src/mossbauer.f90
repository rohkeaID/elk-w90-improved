
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: mossbauer
! !INTERFACE:
subroutine mossbauer
! !USES:
use modmain
! !DESCRIPTION:
!   Computes the contact charge density and contact magnetic hyperfine field for
!   each atom and outputs the data to the file {\tt MOSSBAUER.OUT}. The Thomson
!   radius, $R_{\rm T}=Z/c^2$, is used for determining the contact moment $m_c$,
!   the relation of which to field strength is via Fermi's formula
!   $$ B_c=\frac{8\pi}{3}\mu_B m_c, $$
!   where the orbital and dipolar contributions are neglected. See
!   S. Bl\"{u}gel, H. Akai, R. Zeller, and P. H. Dederichs, {\it Phys. Rev. B}
!   {\bf 35}, 3271 (1987). See also {\tt radnucl}.
!
! !REVISION HISTORY:
!   Created May 2004 (JKD)
!EOP
!BOC
implicit none
! local variables
integer idm,is,ia,ias
integer nr,nri,nrn,nrt,ir
real(8) rt,vt,mc,bc,t1
real(8) rho0,rhon,rhoa
! allocatable arrays
real(8), allocatable :: fr(:,:)
! external functions
real(8) fintgt
external fintgt
! initialise universal variables
call init0
! read density and potentials from file
call readstate
! allocate local arrays
allocate(fr(nrmtmax,3))
open(50,file='MOSSBAUER.OUT',form='FORMATTED')
! loop over species
do is=1,nspecies
  nr=nrmt(is)
  nri=nrmti(is)
  nrn=nrnucl(is)
! Thomson radius and volume
  rt=abs(spzn(is))/solsc**2
  nrt=nr
  do ir=1,nr
    if (rsp(ir,is).gt.rt) then
      nrt=ir
      exit
    end if
  end do
  rt=rsp(nrt,is)
  vt=(4.d0/3.d0)*pi*rt**3
! loop over atoms
  do ia=1,natoms(is)
    ias=idxas(ia,is)
!--------------------------------!
!     contact charge density     !
!--------------------------------!
    call rfmtlm(1,nr,nri,rhomt(:,ias),fr)
    rho0=fr(1,1)*y00
    rhon=fr(nrn,1)*y00
    fr(1:nrn,1)=fr(1:nrn,1)*r2sp(1:nrn,is)
    t1=fintgt(-1,nrn,rsp(:,is),fr)
    rhoa=fourpi*y00*t1/vnucl(is)
    write(50,*)
    write(50,*)
    write(50,'("Species : ",I4," (",A,"), atom : ",I4)') is,trim(spsymb(is)),ia
    write(50,*)
    write(50,'(" approximate nuclear radius : ",G18.10)') rnucl(is)
    write(50,'(" number of mesh points to nuclear radius : ",I6)') nrn
    write(50,'(" density at nuclear center      : ",G18.10)') rho0
    write(50,'(" density at nuclear surface     : ",G18.10)') rhon
    write(50,'(" average contact charge density : ",G18.10)') rhoa
!------------------------------------------!
!     contact magnetic hyperfine field     !
!------------------------------------------!
    if (spinpol) then
      do idm=1,ndmag
        call rfmtlm(1,nr,nri,magmt(:,ias,idm),fr(:,idm))
      end do
      if (ncmag) then
! non-collinear
        fr(1:nrt,1)=sqrt(fr(1:nrt,1)**2 &
                        +fr(1:nrt,2)**2 &
                        +fr(1:nrt,3)**2)*r2sp(1:nrt,is)
      else
! collinear
        fr(1:nrt,1)=abs(fr(1:nrt,1))*r2sp(1:nrt,is)
      end if
      t1=fintgt(-1,nrt,rsp(:,is),fr)
      mc=fourpi*y00*t1/vt
      write(50,*)
      write(50,'(" Thomson radius : ",G18.10)') rt
      write(50,'(" number of mesh points to Thomson radius : ",I6)') nrt
      write(50,'(" contact magnetic moment (mu_B) : ",G18.10)') mc
      bc=(8.d0*pi/3.d0)*mc/(2.d0*solsc)
      write(50,'(" contact hyperfine field : ",G18.10)') bc
      write(50,'("  tesla                  : ",G18.10)') bc*b_si/solsc
    end if
  end do
end do
close(50)
write(*,*)
write(*,'("Info(mossbauer):")')
write(*,'(" Mossbauer parameters written to MOSSBAUER.OUT")')
deallocate(fr)
return
end subroutine
!EOC

