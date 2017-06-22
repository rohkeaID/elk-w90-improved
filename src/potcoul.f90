
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: potcoul
! !INTERFACE:
subroutine potcoul
! !USES:
use modmain
! !DESCRIPTION:
!   Calculates the Coulomb potential of the real charge density stored in the
!   global variables {\tt rhomt} and {\tt rhoir} by solving Poisson's equation.
!   These variables are coverted to complex representations and passed to the
!   routine {\tt zpotcoul}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias,nr
real(8) t1
complex(8) zrho0
! automatic arrays
real(8) vn(nrmtmax)
! allocatable arrays
complex(8), allocatable :: zrhomt(:,:,:),zrhoir(:)
complex(8), allocatable :: zvclmt(:,:,:),zvclir(:)
allocate(zrhomt(lmmaxvr,nrmtmax,natmtot))
! convert real muffin-tin charge density to complex spherical harmonic expansion
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(is)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  call rtozfmt(nrmt(is),nrmtinr(is),1,rhomt(:,:,ias),1,zrhomt(:,:,ias))
end do
!$OMP END DO
!$OMP END PARALLEL
! solve the complex Poisson's equation in the muffin-tins
allocate(zvclmt(lmmaxvr,nrmtmax,natmtot))
call genzvclmt(nrmt,nrmtinr,nrspmax,rsp,nrmtmax,zrhomt,zvclmt)
deallocate(zrhomt)
! add the nuclear monopole potentials
t1=1.d0/y00
do is=1,nspecies
  nr=nrmt(is)
  call potnucl(ptnucl,nr,rsp(:,is),spzn(is),vn)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    zvclmt(1,1:nr,ias)=zvclmt(1,1:nr,ias)+t1*vn(1:nr)
  end do
end do
! store real interstitial charge density in complex array
allocate(zrhoir(ngtot))
zrhoir(:)=rhoir(:)
! solve Poisson's equation in the entire unit cell
allocate(zvclir(ngtot))
call zpotcoul(nrmt,nrmtinr,nrspmax,rsp,1,gc,jlgr,ylmg,sfacg,zrhoir,nrmtmax, &
 zvclmt,zvclir,zrho0)
! convert complex muffin-tin potential to real spherical harmonic expansion
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(is)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  call ztorfmt(nrmt(is),nrmtinr(is),1,zvclmt(:,:,ias),1,vclmt(:,:,ias))
end do
!$OMP END DO
!$OMP END PARALLEL
! store complex interstitial potential in real array
call dcopy(ngtot,zvclir,2,vclir,1)
deallocate(zrhoir,zvclmt,zvclir)
! apply constant electric field if required
if (efieldpol) call potefield
return
end subroutine
!EOC

