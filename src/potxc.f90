
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: potxc
! !INTERFACE:
subroutine potxc
! !USES:
use modmain
use modxcifc
! !DESCRIPTION:
!   Computes the exchange-correlation potential and energy density. In the
!   muffin-tin, the density is transformed from spherical harmonic coefficients
!   $\rho_{lm}$ to spherical coordinates $(\theta,\phi)$ with a backward
!   spherical harmonic transformation (SHT). Once calculated, the
!   exchange-correlation potential and energy density are transformed with a
!   forward SHT.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer idm,is,ia,ias
integer nr,nri,ir,i,n
real(8) t1,t2,t3,t4
! allocatable arrays
real(8), allocatable :: rho(:),rhoup(:),rhodn(:)
real(8), allocatable :: gvrho(:),gvup(:),gvdn(:)
real(8), allocatable :: grho(:),gup(:),gdn(:)
real(8), allocatable :: g2rho(:),g2up(:),g2dn(:)
real(8), allocatable :: g3rho(:),g3up(:),g3dn(:)
real(8), allocatable :: grho2(:),gup2(:),gdn2(:),gupdn(:)
real(8), allocatable :: taumt(:,:,:,:),tauir(:,:)
real(8), allocatable :: ex(:),ec(:),vxc(:)
real(8), allocatable :: vx(:),vxup(:),vxdn(:)
real(8), allocatable :: vc(:),vcup(:),vcdn(:)
real(8), allocatable :: dxdg2(:),dxdgu2(:),dxdgd2(:),dxdgud(:)
real(8), allocatable :: dcdg2(:),dcdgu2(:),dcdgd2(:),dcdgud(:)
real(8), allocatable :: mag(:,:),bxc(:,:)
! meta-GGA variables if required
if (xcgrad.eq.3) then
! generate the kinetic energy density if required
  allocate(taumt(lmmaxvr,nrmtmax,natmtot,nspinor))
  allocate(tauir(ngtot,nspinor))
  call gentau(taumt,tauir)
! compute the Tran-Blaha '09 constant c
  call xc_c_tb09
end if
! allocate local arrays
n=lmmaxvr*nrmtmax
allocate(ex(n),ec(n),vxc(n))
if (spinpol) then
  allocate(mag(n,3),bxc(n,3))
end if
n=max(n,ngtot)
allocate(rho(n))
if (spinpol) then
  allocate(rhoup(n),rhodn(n))
  allocate(vxup(n),vxdn(n),vcup(n),vcdn(n))
  if (xcgrad.eq.1) then
    allocate(grho(n),gup(n),gdn(n))
    allocate(g2up(n),g2dn(n))
    allocate(g3rho(n),g3up(n),g3dn(n))
  else if (xcgrad.eq.2) then
    allocate(g2up(n),g2dn(n))
    allocate(gvup(3*n),gvdn(3*n))
    allocate(gup2(n),gdn2(n),gupdn(n))
    allocate(dxdgu2(n),dxdgd2(n),dxdgud(n))
    allocate(dcdgu2(n),dcdgd2(n),dcdgud(n))
  else if (xcgrad.eq.3) then
    allocate(g2up(n),g2dn(n))
    allocate(gvup(3*n),gvdn(3*n))
    allocate(gup2(n),gdn2(n),gupdn(n))
  end if
else
  allocate(vx(n),vc(n))
  if (xcgrad.eq.1) then
    allocate(grho(n),g2rho(n),g3rho(n))
  else if (xcgrad.eq.2) then
    allocate(g2rho(n),gvrho(3*n),grho2(n))
    allocate(dxdg2(n),dcdg2(n))
  else if (xcgrad.eq.3) then
    allocate(g2rho(n),gvrho(3*n),grho2(n))
  end if
end if
!---------------------------------------!
!     muffin-tin potential and field    !
!---------------------------------------!
do is=1,nspecies
  nr=nrmt(is)
  nri=nrmtinr(is)
  n=lmmaxvr*nr
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! compute the density in spherical coordinates
    call rbsht(nr,nri,1,rhomt(:,:,ias),1,rho)
    if (spinpol) then
!------------------------!
!     spin-polarised     !
!------------------------!
! magnetisation in spherical coordinates
      do idm=1,ndmag
        call rbsht(nr,nri,1,magmt(:,:,ias,idm),1,mag(:,idm))
      end do
      if (ncmag) then
! non-collinear (use Kubler's trick)
        do i=1,n
! compute rhoup=(rho+|m|)/2 and rhodn=(rho-|m|)/2
          t1=sqrt(mag(i,1)**2+mag(i,2)**2+mag(i,3)**2)
          rhoup(i)=0.5d0*(rho(i)+t1)
          rhodn(i)=0.5d0*(rho(i)-t1)
        end do
      else
! collinear
        do i=1,n
! compute rhoup=(rho+m_z)/2 and rhodn=(rho-m_z)/2
          rhoup(i)=0.5d0*(rho(i)+mag(i,1))
          rhodn(i)=0.5d0*(rho(i)-mag(i,1))
        end do
      end if
! call the exchange-correlation interface routine
      if (xcgrad.le.0) then
        call xcifc(xctype,n=n,rhoup=rhoup,rhodn=rhodn,ex=ex,ec=ec,vxup=vxup, &
         vxdn=vxdn,vcup=vcup,vcdn=vcdn)
      else if (xcgrad.eq.1) then
        if (ncgga) then
          rho(1:n)=0.5d0*rho(1:n)
          call ggamt_sp_1(is,rho,rho,grho,gup,gdn,g2up,g2dn,g3rho,g3up,g3dn)
        else
          call ggamt_sp_1(is,rhoup,rhodn,grho,gup,gdn,g2up,g2dn,g3rho,g3up,g3dn)
        end if
        call xcifc(xctype,n=n,rhoup=rhoup,rhodn=rhodn,grho=grho,gup=gup, &
         gdn=gdn,g2up=g2up,g2dn=g2dn,g3rho=g3rho,g3up=g3up,g3dn=g3dn,ex=ex, &
         ec=ec,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
      else if (xcgrad.eq.2) then
        call ggamt_sp_2a(is,rhoup,rhodn,g2up,g2dn,gvup,gvdn,gup2,gdn2,gupdn)
        call xcifc(xctype,n=n,rhoup=rhoup,rhodn=rhodn,gup2=gup2,gdn2=gdn2, &
         gupdn=gupdn,ex=ex,ec=ec,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn, &
         dxdgu2=dxdgu2,dxdgd2=dxdgd2,dxdgud=dxdgud,dcdgu2=dcdgu2, &
         dcdgd2=dcdgd2,dcdgud=dcdgud)
        call ggamt_sp_2b(is,g2up,g2dn,gvup,gvdn,vxup,vxdn,vcup,vcdn,dxdgu2, &
         dxdgd2,dxdgud,dcdgu2,dcdgd2,dcdgud)
      else if (xcgrad.eq.3) then
        call ggamt_sp_2a(is,rhoup,rhodn,g2up,g2dn,gvup,gvdn,gup2,gdn2,gupdn)
        call xcifc(xctype,n=n,c_tb09=c_tb09,rhoup=rhoup,rhodn=rhodn,g2up=g2up, &
         g2dn=g2dn,gup2=gup2,gdn2=gdn2,gupdn=gupdn,tauup=taumt(:,:,ias,1), &
         taudn=taumt(:,:,ias,2),vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
        ex(1:n)=0.d0; ec(1:n)=0.d0
      end if
! hybrid functionals
      if (hybrid) then
        t1=1.d0-hybridc
! scale exchange part of energy
        ex(1:n)=t1*ex(1:n)
! scale exchange part of potential
        vxup(1:n)=t1*vxup(1:n)
        vxdn(1:n)=t1*vxdn(1:n)
      end if
      if (ncmag) then
! non-collinear: locally spin rotate the exchange-correlation potential
        do i=1,n
          t1=vxup(i)+vcup(i)
          t2=vxdn(i)+vcdn(i)
          vxc(i)=0.5d0*(t1+t2)
! determine the exchange-correlation magnetic field
          t3=0.5d0*(t1-t2)
! |m| = rhoup - rhodn
          t4=rhoup(i)-rhodn(i)
          if (abs(t4).gt.1.d-8) t4=t3/t4
          bxc(i,1:3)=mag(i,1:3)*t4
        end do
      else
! collinear
        do i=1,n
          t1=vxup(i)+vcup(i)
          t2=vxdn(i)+vcdn(i)
          vxc(i)=0.5d0*(t1+t2)
          bxc(i,1)=0.5d0*(t1-t2)
        end do
      end if
! convert field to spherical harmonics
      do idm=1,ndmag
        call rfsht(nr,nri,1,bxc(:,idm),1,bxcmt(:,:,ias,idm))
      end do
    else
!--------------------------!
!     spin-unpolarised     !
!--------------------------!
      if (xcgrad.le.0) then
        call xcifc(xctype,n=n,rho=rho,ex=ex,ec=ec,vx=vx,vc=vc)
      else if (xcgrad.eq.1) then
        call ggamt_1(ias,grho,g2rho,g3rho)
        call xcifc(xctype,n=n,rho=rho,grho=grho,g2rho=g2rho,g3rho=g3rho,ex=ex, &
         ec=ec,vx=vx,vc=vc)
      else if (xcgrad.eq.2) then
        call ggamt_2a(ias,g2rho,gvrho,grho2)
        call xcifc(xctype,n=n,rho=rho,grho2=grho2,ex=ex,ec=ec,vx=vx,vc=vc, &
         dxdg2=dxdg2,dcdg2=dcdg2)
        call ggamt_2b(is,g2rho,gvrho,vx,vc,dxdg2,dcdg2)
      else if (xcgrad.eq.3) then
        call ggamt_2a(ias,g2rho,gvrho,grho2)
        call xcifc(xctype,n=n,c_tb09=c_tb09,rho=rho,g2rho=g2rho,grho2=grho2, &
         tau=taumt(:,:,ias,1),vx=vx,vc=vc)
        ex(1:n)=0.d0; ec(1:n)=0.d0
      end if
! hybrid functionals
      if (hybrid) then
        t1=1.d0-hybridc
! scale exchange part of energy
        ex(1:n)=t1*ex(1:n)
! scale exchange part of potential
        vxc(1:n)=t1*vx(1:n)+vc(1:n)
      else
        vxc(1:n)=vx(1:n)+vc(1:n)
      end if
    end if
! convert exchange and correlation energy densities to spherical harmonics
    call rfsht(nr,nri,1,ex,1,exmt(:,:,ias))
    call rfsht(nr,nri,1,ec,1,ecmt(:,:,ias))
! convert exchange-correlation potential to spherical harmonics
    call rfsht(nr,nri,1,vxc,1,vxcmt(:,:,ias))
  end do
end do
!------------------------------------------!
!     interstitial potential and field     !
!------------------------------------------!
if (spinpol) then
!------------------------!
!     spin-polarised     !
!------------------------!
  if (ncmag) then
! non-collinear
    do ir=1,ngtot
      t1=sqrt(magir(ir,1)**2+magir(ir,2)**2+magir(ir,3)**2)
      rhoup(ir)=0.5d0*(rhoir(ir)+t1)
      rhodn(ir)=0.5d0*(rhoir(ir)-t1)
    end do
  else
! collinear
    do ir=1,ngtot
      rhoup(ir)=0.5d0*(rhoir(ir)+magir(ir,1))
      rhodn(ir)=0.5d0*(rhoir(ir)-magir(ir,1))
    end do
  end if
  if (xcgrad.le.0) then
    call xcifc(xctype,n=ngtot,rhoup=rhoup,rhodn=rhodn,ex=exir,ec=ecir, &
     vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
  else if (xcgrad.eq.1) then
    if (ncgga) then
      rho(1:ngtot)=0.5d0*rhoir(1:ngtot)
      call ggair_sp_1(rho,rho,grho,gup,gdn,g2up,g2dn,g3rho,g3up,g3dn)
    else
      call ggair_sp_1(rhoup,rhodn,grho,gup,gdn,g2up,g2dn,g3rho,g3up,g3dn)
    end if
    call xcifc(xctype,n=ngtot,rhoup=rhoup,rhodn=rhodn,grho=grho,gup=gup, &
     gdn=gdn,g2up=g2up,g2dn=g2dn,g3rho=g3rho,g3up=g3up,g3dn=g3dn,ex=exir, &
     ec=ecir,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
  else if (xcgrad.eq.2) then
    call ggair_sp_2a(rhoup,rhodn,g2up,g2dn,gvup,gvdn,gup2,gdn2,gupdn)
    if (ncgga) then
      g2up(1:ngtot)=0.5d0*(g2up(1:ngtot)+g2dn(1:ngtot))
      g2dn(1:ngtot)=g2up(1:ngtot)
    end if
    call xcifc(xctype,n=ngtot,rhoup=rhoup,rhodn=rhodn,gup2=gup2,gdn2=gdn2, &
     gupdn=gupdn,ex=exir,ec=ecir,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn, &
     dxdgu2=dxdgu2,dxdgd2=dxdgd2,dxdgud=dxdgud,dcdgu2=dcdgu2,dcdgd2=dcdgd2, &
     dcdgud=dcdgud)
    call ggair_sp_2b(g2up,g2dn,gvup,gvdn,vxup,vxdn,vcup,vcdn,dxdgu2,dxdgd2, &
     dxdgud,dcdgu2,dcdgd2,dcdgud)
  else if (xcgrad.eq.3) then
    call ggair_sp_2a(rhoup,rhodn,g2up,g2dn,gvup,gvdn,gup2,gdn2,gupdn)
    call xcifc(xctype,n=ngtot,c_tb09=c_tb09,rhoup=rhoup,rhodn=rhodn,g2up=g2up, &
     g2dn=g2dn,gup2=gup2,gdn2=gdn2,gupdn=gupdn,tauup=tauir(:,1), &
     taudn=tauir(:,2),vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
    exir(1:ngtot)=0.d0; ecir(1:ngtot)=0.d0
  end if
! hybrid functionals
  if (hybrid) then
    t1=1.d0-hybridc
! scale exchange part of energy
    exir(1:ngtot)=t1*exir(1:ngtot)
! scale exchange part of potential
    vxup(1:ngtot)=t1*vxup(1:ngtot)
    vxdn(1:ngtot)=t1*vxdn(1:ngtot)
  end if
  if (ncmag) then
! non-collinear: spin rotate the local exchange potential
    do ir=1,ngtot
      t1=vxup(ir)+vcup(ir)
      t2=vxdn(ir)+vcdn(ir)
      vxcir(ir)=0.5d0*(t1+t2)
! determine the exchange-correlation magnetic field
      t3=0.5d0*(t1-t2)
      t4=rhoup(ir)-rhodn(ir)
      if (abs(t4).gt.1.d-8) t4=t3/t4
      bxcir(ir,:)=magir(ir,:)*t4
    end do
  else
! collinear
    do ir=1,ngtot
      t1=vxup(ir)+vcup(ir)
      t2=vxdn(ir)+vcdn(ir)
      vxcir(ir)=0.5d0*(t1+t2)
      bxcir(ir,1)=0.5d0*(t1-t2)
    end do
  end if
else
!--------------------------!
!     spin-unpolarised     !
!--------------------------!
  if (xcgrad.le.0) then
    call xcifc(xctype,n=ngtot,rho=rhoir,ex=exir,ec=ecir,vx=vx,vc=vc)
  else if (xcgrad.eq.1) then
    call ggair_1(grho,g2rho,g3rho)
    call xcifc(xctype,n=ngtot,rho=rhoir,grho=grho,g2rho=g2rho,g3rho=g3rho, &
     ex=exir,ec=ecir,vx=vx,vc=vc)
  else if (xcgrad.eq.2) then
    call ggair_2a(g2rho,gvrho,grho2)
    call xcifc(xctype,n=ngtot,rho=rhoir,grho2=grho2,ex=exir,ec=ecir,vx=vx, &
     vc=vc,dxdg2=dxdg2,dcdg2=dcdg2)
    call ggair_2b(g2rho,gvrho,vx,vc,dxdg2,dcdg2)
  else if (xcgrad.eq.3) then
    call ggair_2a(g2rho,gvrho,grho2)
    call xcifc(xctype,n=ngtot,c_tb09=c_tb09,rho=rhoir,g2rho=g2rho,grho2=grho2, &
     tau=tauir(:,1),vx=vx,vc=vc)
    exir(1:ngtot)=0.d0; ecir(1:ngtot)=0.d0
  end if
! hybrid functionals
  if (hybrid) then
    t1=1.d0-hybridc
! scale exchange part of energy
    exir(1:ngtot)=t1*exir(1:ngtot)
! scale exchange part of potential
    vxcir(1:ngtot)=t1*vx(1:ngtot)+vc(1:ngtot)
  else
    vxcir(1:ngtot)=vx(1:ngtot)+vc(1:ngtot)
  end if
end if
! optimised effective potential
if (xctype(1).lt.0) call oepmain
! symmetrise the exchange-correlation potential
call symrf(1,vxcmt,vxcir)
if (spinpol) then
! symmetrise the exchange-correlation effective field
  call symrvf(1,bxcmt,bxcir)
! remove the source contribution if required
  if (nosource) call projsbf
end if
if (xcgrad.eq.3) deallocate(taumt,tauir)
deallocate(rho,ex,ec,vxc)
if (spinpol) then
  deallocate(mag,bxc)
  deallocate(rhoup,rhodn,vxup,vxdn,vcup,vcdn)
  if (xcgrad.eq.1) then
    deallocate(grho,gup,gdn,g2up,g2dn,g3rho,g3up,g3dn)
  else if (xcgrad.eq.2) then
    deallocate(g2up,g2dn)
    deallocate(gvup,gvdn)
    deallocate(gup2,gdn2,gupdn)
    deallocate(dxdgu2,dxdgd2,dxdgud)
    deallocate(dcdgu2,dcdgd2,dcdgud)
  else if (xcgrad.eq.3) then
    deallocate(g2up,g2dn)
    deallocate(gvup,gvdn)
    deallocate(gup2,gdn2,gupdn)
  end if
else
  deallocate(vx,vc)
  if (xcgrad.eq.1) then
    deallocate(grho,g2rho,g3rho)
  else if (xcgrad.eq.2) then
    deallocate(g2rho,gvrho,grho2)
    deallocate(dxdg2,dcdg2)
  else if (xcgrad.eq.3) then
    deallocate(g2rho,gvrho,grho2)
  end if
end if
return
end subroutine
!EOC

