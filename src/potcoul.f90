
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: potcoul
! !INTERFACE:
subroutine potcoul
! !USES:
use modmain
use modomp
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
integer is,ia,ias,nthd
integer nr,nri,ir,i
real(8) t1
! automatic arrays
real(8) vn(nrmtmax)
! allocatable arrays
complex(8), allocatable :: zrhomt(:,:),zrhoir(:)
complex(8), allocatable :: zvclmt(:,:),zvclir(:)
allocate(zrhomt(npmtmax,natmtot))
! convert real muffin-tin charge density to complex spherical harmonic expansion
call omp_hold(natmtot,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(is) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  call rtozfmt(nrmt(is),nrmti(is),rhomt(:,ias),zrhomt(:,ias))
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
! solve the complex Poisson's equation in the muffin-tins
allocate(zvclmt(npmtmax,natmtot))
call genzvclmt(nrmt,nrmti,nrspmax,rsp,npmtmax,zrhomt,zvclmt)
deallocate(zrhomt)
! add the nuclear monopole potentials
t1=1.d0/y00
do is=1,nspecies
  nr=nrmt(is)
  nri=nrmti(is)
  call potnucl(ptnucl,nr,rsp(:,is),spzn(is),vn)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    i=1
    do ir=1,nri
      zvclmt(i,ias)=zvclmt(i,ias)+t1*vn(ir)
      i=i+lmmaxi
    end do
    do ir=nri+1,nr
      zvclmt(i,ias)=zvclmt(i,ias)+t1*vn(ir)
      i=i+lmmaxo
    end do
  end do
end do
! store real interstitial charge density in complex array
allocate(zrhoir(ngtot))
zrhoir(:)=rhoir(:)
! solve Poisson's equation in the entire unit cell
allocate(zvclir(ngtot))
call zpotcoul(nrmt,nrmti,npmt,npmti,nrspmax,rsp,ngridg,igfft,ngvec,gc,gclg, &
 ngvec,jlgrmt,ylmg,sfacg,zrhoir,npmtmax,zvclmt,zvclir)
! convert complex muffin-tin potential to real spherical harmonic expansion
call omp_hold(natmtot,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(is) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  call ztorfmt(nrmt(is),nrmti(is),zvclmt(:,ias),vclmt(:,ias))
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
! store complex interstitial potential in real array
call dcopy(ngtot,zvclir,2,vclir,1)
deallocate(zrhoir,zvclmt,zvclir)
! apply constant electric field if required
if (tefield) call potefield
return
end subroutine
!EOC

