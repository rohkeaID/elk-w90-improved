
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dforce(dyn)
use modmain
use modphonon
implicit none
! arguments
complex(8), intent(out) :: dyn(3,natmtot)
! local variables
integer ik,is,ias
integer nr,nri,ir
integer igq0,i
complex(8) zrho0,zsum,z1
! automatic arrays
real(8) vn(nrmtmax)
! allocatable arrays
complex(8), allocatable :: zrhomt(:,:,:),zrhoir(:)
complex(8), allocatable :: grhomt(:,:,:,:),grhoir(:,:)
complex(8), allocatable :: zvclmt(:,:,:),zvclir(:)
complex(8), allocatable :: gvclmt(:,:,:,:),gvclir(:,:)
complex(8), allocatable :: zfmt(:,:),gzfmt(:,:,:)
! external functions
complex(8) zfinp,zfmtinp
external zfinp,zfmtinp
allocate(zrhomt(lmmaxvr,nrmtmax,natmtot),zrhoir(ngtot))
allocate(grhomt(lmmaxvr,nrmtmax,natmtot,3),grhoir(ngtot,3))
allocate(zvclmt(lmmaxvr,nrmtmax,natmtot),zvclir(ngtot))
allocate(gvclmt(lmmaxvr,nrmtmax,natmtot,3),gvclir(ngtot,3))
allocate(zfmt(lmmaxvr,nrmtmax),gzfmt(lmmaxvr,nrmtmax,3))
! make complex copy of the density
do ias=1,natmtot
  is=idxis(ias)
  call rtozfmt(nrmt(is),nrmtinr(is),1,rhomt(:,:,ias),1,zrhomt(:,:,ias))
end do
zrhoir(:)=rhoir(:)
! compute the gradient of the density
call gradzf(zrhomt,zrhoir,grhomt,grhoir)
!--------------------------------------------------------------!
!     Hellmann-Feynman force derivative for displaced atom     !
!--------------------------------------------------------------!
if (iqph.eq.iq0) then
  igq0=1
else
  igq0=0
end if
nr=nrmt(isph)
! zero the interstitial density
zrhoir(:)=0.d0
! compute the gradient of the nuclear potential
call potnucl(ptnucl,nr,rsp(:,isph),spzn(isph),vn)
do ir=1,nr
  zfmt(:,ir)=0.d0
  zfmt(1,ir)=vn(ir)/y00
end do
call gradzfmt(nr,nrmtinr(isph),rsp(:,isph),zfmt,nrmtmax,gzfmt)
! compute the q-dependent nuclear Coulomb potential derivative
zvclmt(:,:,:)=0.d0
do ir=1,nr
  zvclmt(2:4,ir,iasph)=gzfmt(2:4,ir,ipph)
end do
tphdyn=.true.
call zpotcoul(nrmt,nrmtinr,nrspmax,rsp,igq0,gqc,jlgqr,ylmgq,sfacgq,zrhoir, &
 nrmtmax,zvclmt,zvclir,zrho0)
zfmt(:,:)=zvnmt(:,:)
! multiply with density derivative and integrate
zsum=0.d0
do ir=1,ngtot
  zsum=zsum+cfunir(ir)*conjg(zvclir(ir))*drhoir(ir)
end do
zsum=zsum*omega/dble(ngtot)
do ias=1,natmtot
  is=idxis(ias)
  zsum=zsum+zfmtinp(nrmt(is),nrmtinr(is),rsp(:,is),r2sp(:,is),zvclmt(:,:,ias), &
   drhomt(:,:,ias))
end do
dyn(ipph,iasph)=-zsum
! compute the lattice-periodic nuclear Coulomb potential derivative
zvclmt(:,:,:)=0.d0
do ir=1,nr
  zvclmt(2:4,ir,iasph)=gzfmt(2:4,ir,ipph)
end do
call zpotcoul(nrmt,nrmtinr,nrspmax,rsp,1,gc,jlgr,ylmg,sfacg,zrhoir,nrmtmax, &
 zvclmt,zvclir,zrho0)
tphdyn=.false.
! multiply with density gradient and integrate
zsum=0.d0
do ir=1,ngtot
  zsum=zsum+cfunir(ir)*zvclir(ir)*grhoir(ir,ipph)
end do
zsum=zsum*omega/dble(ngtot)
do ias=1,natmtot
  is=idxis(ias)
  zsum=zsum+zfmtinp(nrmt(is),nrmtinr(is),rsp(:,is),r2sp(:,is),zvclmt(:,:,ias), &
   grhomt(:,:,ias,ipph))
end do
dyn(ipph,iasph)=dyn(ipph,iasph)-zsum
! nuclear-nuclear term
zvclmt(:,:,iasph)=zvnmt(:,:)-zfmt(:,:)
call gradzf(zvclmt,zvclir,gvclmt,gvclir)
do ias=1,natmtot
  is=idxis(ias)
  z1=spzn(is)*gvclmt(1,nrnucl(is),ias,ipph)*y00
  dyn(ipph,iasph)=dyn(ipph,iasph)+z1
end do
!-------------------------------------------------------------------!
!     Hellmann-Feynman force derivative for non-displaced atoms     !
!-------------------------------------------------------------------!
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmtinr(is)
! remove the gradient part of the Coulomb potential for displaced muffin-tin
  if (ias.eq.iasph) then
    call rtozfmt(nr,nri,1,vclmt(:,:,iasph),1,zfmt)
    call gradzfmt(nr,nri,rsp(:,isph),zfmt,nrmtmax,gzfmt)
    dvclmt(:,1:nr,ias)=dvclmt(:,1:nr,ias)+gzfmt(:,1:nr,ipph)
  end if
! compute the gradient of the Coulomb potential derivative at the nucleus
  call gradzfmt(nr,nri,rsp(:,is),dvclmt(:,:,ias),nrmtmax,gzfmt)
  do i=1,3
    if ((ias.eq.iasph).and.(i.eq.ipph)) cycle
    dyn(i,ias)=spzn(is)*gzfmt(1,nrnucl(is),i)*y00
  end do
end do
!--------------------------------------------!
!     IBS correction to force derivative     !
!--------------------------------------------!
if (tfibs) then
! k-point dependent part
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
  do ik=1,nkptnr
    call dforcek(ik,dyn)
  end do
!$OMP END DO
!$OMP END PARALLEL
! k-point independent part
  do ias=1,natmtot
    is=idxis(ias)
    nr=nrmt(is)
    nri=nrmtinr(is)
    do i=1,3
      z1=zfmtinp(nr,nri,rsp(:,is),r2sp(:,is),grhomt(:,:,ias,i),dvsmt(:,:,ias))
      dyn(i,ias)=dyn(i,ias)-z1
    end do
! convert Kohn-Sham potential to complex spherical harmonics
    call rtozfmt(nr,nri,1,vsmt(:,:,ias),1,zfmt)
! remove the gradient part from the density derivative for displaced muffin-tin
    if (ias.eq.iasph) then
      drhomt(:,1:nr,ias)=drhomt(:,1:nr,ias)+grhomt(:,1:nr,ias,ipph)
    end if
! compute the gradient of the density derivative
    call gradzfmt(nr,nri,rsp(:,is),drhomt(:,:,ias),nrmtmax,gzfmt)
    do i=1,3
      z1=zfmtinp(nr,nri,rsp(:,is),r2sp(:,is),zfmt,gzfmt(:,:,i))
      dyn(i,ias)=dyn(i,ias)-z1
    end do
  end do
end if
deallocate(zrhomt,zrhoir,grhomt,grhoir)
deallocate(zvclmt,zvclir,gvclmt,gvclir,zfmt,gzfmt)
return
end subroutine

