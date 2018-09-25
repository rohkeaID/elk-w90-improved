
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dforce(dyn)
use modmain
use modphonon
use modomp
implicit none
! arguments
complex(8), intent(out) :: dyn(3,natmtot)
! local variables
integer ik,is,ias,nthd
integer nr,nri,ir,np,i
complex(8) zsum,z1
! automatic arrays
real(8) vn(nrmtmax)
! allocatable arrays
complex(8), allocatable :: zrhomt(:,:),zrhoir(:)
complex(8), allocatable :: grhomt(:,:,:),grhoir(:,:)
complex(8), allocatable :: zvclmt(:,:),zvclir(:)
complex(8), allocatable :: gvclmt(:,:,:),gvclir(:,:)
complex(8), allocatable :: zfmt(:),gzfmt(:,:)
! external functions
real(8) fintgt
complex(8) zfmtinp
external fintgt,zfmtinp
allocate(zrhomt(npmtmax,natmtot),zrhoir(ngtot))
allocate(grhomt(npmtmax,natmtot,3),grhoir(ngtot,3))
allocate(zvclmt(npmtmax,natmtot),zvclir(ngtot))
allocate(gvclmt(npmtmax,natmtot,3),gvclir(ngtot,3))
allocate(zfmt(npmtmax),gzfmt(npmtmax,3))
! make complex copy of the density
do ias=1,natmtot
  is=idxis(ias)
  call rtozfmt(nrmt(is),nrmti(is),rhomt(:,ias),zrhomt(:,ias))
end do
zrhoir(:)=rhoir(:)
! compute the gradient of the density
call gradzf(zrhomt,zrhoir,grhomt,grhoir)
!--------------------------------------------------------------!
!     Hellmann-Feynman force derivative for displaced atom     !
!--------------------------------------------------------------!
nr=nrmt(isph)
nri=nrmti(isph)
np=npmt(isph)
! zero the interstitial density
zrhoir(:)=0.d0
! compute the gradient of the nuclear potential
call potnucl(.false.,nr,rsp(:,isph),spzn(isph),vn)
zfmt(1:np)=0.d0
i=1
do ir=1,nri
  zfmt(i)=vn(ir)/y00
  i=i+lmmaxi
end do
do ir=nri+1,nr
  zfmt(i)=vn(ir)/y00
  i=i+lmmaxo
end do
call gradzfmt(nr,nri,rsp(:,isph),zfmt,npmtmax,gzfmt)
! compute the q-dependent nuclear Coulomb potential derivative
zvclmt(:,:)=0.d0
zvclmt(1:np,iasph)=gzfmt(1:np,ipph)
tphdyn=.true.
call zpotcoul(nrmt,nrmti,npmt,npmti,nrspmax,rsp,ngridg,igfft,ngvec,gqc,gclgq, &
 ngvec,jlgqrmt,ylmgq,sfacgq,zrhoir,npmtmax,zvclmt,zvclir)
zfmt(1:np)=zvnmt(1:np)
! multiply with density derivative and integrate
zsum=0.d0
do ir=1,ngtot
  zsum=zsum+cfunir(ir)*conjg(zvclir(ir))*drhoir(ir)
end do
zsum=zsum*omega/dble(ngtot)
do ias=1,natmtot
  is=idxis(ias)
  zsum=zsum+zfmtinp(nrmt(is),nrmti(is),rsp(:,is),r2sp(:,is),zvclmt(:,ias), &
   drhomt(:,ias))
end do
dyn(ipph,iasph)=-zsum
! compute the lattice-periodic nuclear Coulomb potential derivative
zvclmt(:,:)=0.d0
zvclmt(1:np,iasph)=gzfmt(1:np,ipph)
call zpotcoul(nrmt,nrmti,npmt,npmti,nrspmax,rsp,ngridg,igfft,ngvec,gc,gclg, &
 ngvec,jlgrmt,ylmg,sfacg,zrhoir,npmtmax,zvclmt,zvclir)
tphdyn=.false.
! multiply with density gradient and integrate
zsum=0.d0
do ir=1,ngtot
  zsum=zsum+cfunir(ir)*zvclir(ir)*grhoir(ir,ipph)
end do
zsum=zsum*omega/dble(ngtot)
do ias=1,natmtot
  is=idxis(ias)
  zsum=zsum+zfmtinp(nrmt(is),nrmti(is),rsp(:,is),r2sp(:,is),zvclmt(:,ias), &
   grhomt(:,ias,ipph))
end do
dyn(ipph,iasph)=dyn(ipph,iasph)-zsum
! nuclear-nuclear term
zvclmt(1:np,iasph)=zvnmt(1:np)-zfmt(1:np)
call gradzf(zvclmt,zvclir,gvclmt,gvclir)
do ias=1,natmtot
  is=idxis(ias)
  z1=spzn(is)*gvclmt(1,ias,ipph)*y00
  dyn(ipph,iasph)=dyn(ipph,iasph)+z1
end do
!-------------------------------------------------------------------!
!     Hellmann-Feynman force derivative for non-displaced atoms     !
!-------------------------------------------------------------------!
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  np=npmt(is)
! remove the gradient part of the Coulomb potential for displaced muffin-tin
  if (ias.eq.iasph) then
    call rtozfmt(nr,nri,vclmt(:,iasph),zfmt)
    call gradzfmt(nr,nri,rsp(:,isph),zfmt,npmtmax,gzfmt)
    dvclmt(1:np,ias)=dvclmt(1:np,ias)+gzfmt(1:np,ipph)
  end if
! compute the gradient of the Coulomb potential derivative at the nucleus
  call gradzfmt(nr,nri,rsp(:,is),dvclmt(:,ias),npmtmax,gzfmt)
  do i=1,3
    if ((ias.eq.iasph).and.(i.eq.ipph)) cycle
    dyn(i,ias)=spzn(is)*gzfmt(1,i)*y00
  end do
end do
!--------------------------------------------!
!     IBS correction to force derivative     !
!--------------------------------------------!
! k-point dependent part
call omp_hold(nkptnr,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ik=1,nkptnr
  call dforcek(ik,dyn)
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
! k-point independent part
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  np=npmt(is)
  do i=1,3
    z1=zfmtinp(nr,nri,rsp(:,is),r2sp(:,is),grhomt(:,ias,i),dvsmt(:,ias))
    dyn(i,ias)=dyn(i,ias)-z1
  end do
! convert Kohn-Sham potential to complex spherical harmonics
  call rtozfmt(nr,nri,vsmt(:,ias),zfmt)
! remove the gradient part from the density derivative for displaced muffin-tin
  if (ias.eq.iasph) then
    drhomt(1:np,ias)=drhomt(1:np,ias)+grhomt(1:np,ias,ipph)
  end if
! compute the gradient of the density derivative
  call gradzfmt(nr,nri,rsp(:,is),drhomt(:,ias),npmtmax,gzfmt)
  do i=1,3
    z1=zfmtinp(nr,nri,rsp(:,is),r2sp(:,is),zfmt,gzfmt(:,i))
    dyn(i,ias)=dyn(i,ias)-z1
  end do
end do
deallocate(zrhomt,zrhoir,grhomt,grhoir)
deallocate(zvclmt,zvclir,gvclmt,gvclir,zfmt,gzfmt)
return
end subroutine

