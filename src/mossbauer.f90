
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
integer is,ia,ias,ir,nrt
real(8) rt,vt,rho0,mc,bc,t1
! allocatable arrays
real(8), allocatable :: fr(:),gr(:)
! initialise universal variables
call init0
! read density and potentials from file
call readstate
! allocate local arrays
allocate(fr(nrmtmax),gr(nrmtmax))
open(50,file='MOSSBAUER.OUT',action='WRITE',form='FORMATTED')
! loop over species
do is=1,nspecies
! Thomson radius and volume
  rt=abs(spzn(is))/solsc**2
  do ir=1,nrmt(is)-1
    if (rsp(ir,is).gt.rt) exit
  end do
  nrt=ir
  rt=rsp(nrt,is)
  vt=(4.d0/3.d0)*pi*rt**3
! loop over atoms
  do ia=1,natoms(is)
    ias=idxas(ia,is)
!--------------------------------!
!     contact charge density     !
!--------------------------------!
    do ir=1,nrnucl(is)
      fr(ir)=rhomt(1,ir,ias)*r2sp(ir,is)
    end do
    call fderiv(-1,nrnucl(is),rsp(:,is),fr,gr)
    rho0=fourpi*y00*gr(nrnucl(is))/vnucl(is)
    write(50,*)
    write(50,*)
    write(50,'("Species : ",I4," (",A,"), atom : ",I4)') is,trim(spsymb(is)),ia
    write(50,*)
    write(50,'(" approximate nuclear radius : ",G18.10)') rnucl(is)
    write(50,'(" number of mesh points to nuclear radius : ",I6)') nrnucl(is)
    write(50,'(" density at nuclear center      : ",G18.10)') rhomt(1,1,ias)*y00
    write(50,'(" density at nuclear surface     : ",G18.10)') &
     rhomt(1,nrnucl(is),ias)*y00
    write(50,'(" average contact charge density : ",G18.10)') rho0
!------------------------------------------!
!     contact magnetic hyperfine field     !
!------------------------------------------!
    if (spinpol) then
      do ir=1,nrt
        if (ncmag) then
! non-collinear
          t1=sqrt(magmt(1,ir,ias,1)**2 &
                 +magmt(1,ir,ias,2)**2 &
                 +magmt(1,ir,ias,3)**2)
        else
! collinear
          t1=magmt(1,ir,ias,1)
        end if
        fr(ir)=t1*r2sp(ir,is)
      end do
      call fderiv(-1,nrt,rsp(:,is),fr,gr)
      mc=fourpi*y00*gr(nrt)/vt
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
deallocate(fr,gr)
return
end subroutine
!EOC
