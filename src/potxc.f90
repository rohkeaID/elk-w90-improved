
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
integer ispn,idm,is,ia,ias
integer nr,nri,ir,np,i,n
real(8) t0,t1,t2,t3,t4
! allocatable arrays
real(8), allocatable :: rho(:),rhoup(:),rhodn(:)
real(8), allocatable :: gvrho(:),gvup(:),gvdn(:)
real(8), allocatable :: grho(:),gup(:),gdn(:)
real(8), allocatable :: g2rho(:),g2up(:),g2dn(:)
real(8), allocatable :: g3rho(:),g3up(:),g3dn(:)
real(8), allocatable :: grho2(:),gup2(:),gdn2(:),gupdn(:)
real(8), allocatable :: ex(:),ec(:),vxc(:)
real(8), allocatable :: vx(:),vxup(:),vxdn(:)
real(8), allocatable :: vc(:),vcup(:),vcdn(:)
real(8), allocatable :: mag(:,:),bxc(:,:),tau(:,:)
real(8), allocatable :: dxdgr2(:),dxdgu2(:),dxdgd2(:),dxdgud(:)
real(8), allocatable :: dcdgr2(:),dcdgu2(:),dcdgd2(:),dcdgud(:)
real(8), allocatable :: dxdg2r(:),dxdg2u(:),dxdg2d(:)
real(8), allocatable :: dcdg2r(:),dcdg2u(:),dcdg2d(:)
real(8), allocatable :: wx(:),wxup(:),wxdn(:)
real(8), allocatable :: wc(:),wcup(:),wcdn(:)
n=max(npmtmax,ngtot)
! meta-GGA variables if required
if ((xcgrad.eq.3).or.(xcgrad.eq.4)) then
! generate the kinetic energy density
  call gentau
! compute the Tran-Blaha '09 constant c
  if (xcgrad.eq.3) call xc_c_tb09
  allocate(tau(npmtmax,nspinor))
end if
! allocate local arrays
allocate(rho(npmtmax),ex(npmtmax),ec(npmtmax),vxc(npmtmax))
if (spinpol) then
  allocate(mag(npmtmax,3),bxc(npmtmax,3))
end if
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
  else if (xcgrad.eq.4) then
    allocate(g2up(n),g2dn(n))
    allocate(gvup(3*n),gvdn(3*n))
    allocate(gup2(n),gdn2(n),gupdn(n))
    allocate(dxdgu2(n),dxdgd2(n),dxdgud(n))
    allocate(dcdgu2(n),dcdgd2(n),dcdgud(n))
    allocate(dxdg2u(n),dxdg2d(n))
    allocate(dcdg2u(n),dcdg2d(n))
    allocate(wxup(n),wxdn(n),wcup(n),wcdn(n))
  end if
else
  allocate(vx(n),vc(n))
  if (xcgrad.eq.1) then
    allocate(grho(n),g2rho(n),g3rho(n))
  else if (xcgrad.eq.2) then
    allocate(g2rho(n),gvrho(3*n),grho2(n))
    allocate(dxdgr2(n),dcdgr2(n))
  else if (xcgrad.eq.3) then
    allocate(g2rho(n),gvrho(3*n),grho2(n))
  else if (xcgrad.eq.4) then
    allocate(g2rho(n),gvrho(3*n),grho2(n))
    allocate(dxdgr2(n),dcdgr2(n))
    allocate(dxdg2r(n),dcdg2r(n))
    allocate(wx(n),wc(n))
  end if
end if
!---------------------------------------!
!     muffin-tin potential and field    !
!---------------------------------------!
do is=1,nspecies
  nr=nrmt(is)
  nri=nrmti(is)
  np=npmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! convert the density to spherical coordinates
    call rbsht(nr,nri,rhomt(:,ias),rho)
! convert tau to spherical coordinates if required
    if ((xcgrad.eq.3).or.(xcgrad.eq.4)) then
      do ispn=1,nspinor
        call rbsht(nr,nri,taumt(:,ias,ispn),tau(:,ispn))
      end do
    end if
    if (spinpol) then
!------------------------!
!     spin-polarised     !
!------------------------!
! magnetisation in spherical coordinates
      do idm=1,ndmag
        call rbsht(nr,nri,magmt(:,ias,idm),mag(:,idm))
      end do
! use scaled spin exchange-correlation (SSXC) if required
      if (tssxc) mag(1:np,1:ndmag)=mag(1:np,1:ndmag)*ssxc
      if (ncmag) then
! non-collinear (use Kubler's trick)
        if (xcgrad.eq.0) then
! LSDA
          do i=1,np
! compute rhoup=(rho+|m|)/2 and rhodn=(rho-|m|)/2
            t0=rho(i)
            t1=sqrt(mag(i,1)**2+mag(i,2)**2+mag(i,3)**2)
            rhoup(i)=0.5d0*(t0+t1)
            rhodn(i)=0.5d0*(t0-t1)
          end do
        else
! functionals which require gradients
          do i=1,np
            t0=rho(i)
            t1=sqrt(mag(i,1)**2+mag(i,2)**2+mag(i,3)**2+dncgga)
            rhoup(i)=0.5d0*(t0+t1)
            rhodn(i)=0.5d0*(t0-t1)
          end do
        end if
      else
! collinear
        do i=1,np
! compute rhoup=(rho+m_z)/2 and rhodn=(rho-m_z)/2
          t0=rho(i)
          t1=mag(i,1)
          rhoup(i)=0.5d0*(t0+t1)
          rhodn(i)=0.5d0*(t0-t1)
        end do
      end if
! call the exchange-correlation interface routine
      if (xcgrad.le.0) then
        call xcifc(xctype,n=np,rhoup=rhoup,rhodn=rhodn,ex=ex,ec=ec,vxup=vxup, &
         vxdn=vxdn,vcup=vcup,vcdn=vcdn)
      else if (xcgrad.eq.1) then
        call ggamt_sp_1(is,rhoup,rhodn,grho,gup,gdn,g2up,g2dn,g3rho,g3up,g3dn)
        call xcifc(xctype,n=np,rhoup=rhoup,rhodn=rhodn,grho=grho,gup=gup, &
         gdn=gdn,g2up=g2up,g2dn=g2dn,g3rho=g3rho,g3up=g3up,g3dn=g3dn,ex=ex, &
         ec=ec,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
      else if (xcgrad.eq.2) then
        call ggamt_sp_2a(is,rhoup,rhodn,g2up,g2dn,gvup,gvdn,gup2,gdn2,gupdn)
        call xcifc(xctype,n=np,rhoup=rhoup,rhodn=rhodn,gup2=gup2,gdn2=gdn2, &
         gupdn=gupdn,ex=ex,ec=ec,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn, &
         dxdgu2=dxdgu2,dxdgd2=dxdgd2,dxdgud=dxdgud,dcdgu2=dcdgu2, &
         dcdgd2=dcdgd2,dcdgud=dcdgud)
        call ggamt_sp_2b(is,g2up,g2dn,gvup,gvdn,vxup,vxdn,vcup,vcdn,dxdgu2, &
         dxdgd2,dxdgud,dcdgu2,dcdgd2,dcdgud)
      else if (xcgrad.eq.3) then
        call ggamt_sp_2a(is,rhoup,rhodn,g2up,g2dn,gvup,gvdn,gup2,gdn2,gupdn)
        call xcifc(xctype,n=np,c_tb09=c_tb09,rhoup=rhoup,rhodn=rhodn, &
         g2up=g2up,g2dn=g2dn,gup2=gup2,gdn2=gdn2,gupdn=gupdn,tauup=tau(:,1), &
         taudn=tau(:,2),vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
        ex(1:np)=0.d0; ec(1:np)=0.d0
      else if (xcgrad.eq.4) then
        call ggamt_sp_2a(is,rhoup,rhodn,g2up,g2dn,gvup,gvdn,gup2,gdn2,gupdn)
        call xcifc(xctype,n=np,rhoup=rhoup,rhodn=rhodn,g2up=g2up,g2dn=g2dn, &
         gup2=gup2,gdn2=gdn2,gupdn=gupdn,tauup=tau(:,1),taudn=tau(:,2),ex=ex, &
         ec=ec,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn,dxdgu2=dxdgu2, &
         dxdgd2=dxdgd2,dxdgud=dxdgud,dcdgu2=dcdgu2,dcdgd2=dcdgd2, &
         dcdgud=dcdgud,dxdg2u=dxdg2u,dxdg2d=dxdg2d,dcdg2u=dcdg2u, &
         dcdg2d=dcdg2d,wxup=wxup,wxdn=wxdn,wcup=wcup,wcdn=wcdn)
        call ggamt_sp_2b(is,g2up,g2dn,gvup,gvdn,vxup,vxdn,vcup,vcdn,dxdgu2, &
         dxdgd2,dxdgud,dcdgu2,dcdgd2,dcdgud)
        wxup(1:np)=0.5d0*(wxup(1:np)+wxdn(1:np)+wcup(1:np)+wcdn(1:np))
        call rfsht(nr,nri,wxup,wxcmt(:,ias))
      end if
! hybrid functionals
      if (hybrid) then
        t1=1.d0-hybridc
! scale exchange part of energy
        ex(1:np)=t1*ex(1:np)
! scale exchange part of potential
        vxup(1:np)=t1*vxup(1:np)
        vxdn(1:np)=t1*vxdn(1:np)
      end if
      if (ncmag) then
! non-collinear: locally spin rotate the exchange-correlation potential
        do i=1,np
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
        do i=1,np
          t1=vxup(i)+vcup(i)
          t2=vxdn(i)+vcdn(i)
          vxc(i)=0.5d0*(t1+t2)
          bxc(i,1)=0.5d0*(t1-t2)
        end do
      end if
! scale B_xc for SSXC if required
      if (tssxc) bxc(1:np,1:ndmag)=bxc(1:np,1:ndmag)*ssxc
! convert field to spherical harmonics
      do idm=1,ndmag
        call rfsht(nr,nri,bxc(:,idm),bxcmt(:,ias,idm))
      end do
    else
!--------------------------!
!     spin-unpolarised     !
!--------------------------!
      if (xcgrad.le.0) then
        call xcifc(xctype,n=np,rho=rho,ex=ex,ec=ec,vx=vx,vc=vc)
      else if (xcgrad.eq.1) then
        call ggamt_1(ias,grho,g2rho,g3rho)
        call xcifc(xctype,n=np,rho=rho,grho=grho,g2rho=g2rho,g3rho=g3rho, &
         ex=ex,ec=ec,vx=vx,vc=vc)
      else if (xcgrad.eq.2) then
        call ggamt_2a(ias,g2rho,gvrho,grho2)
        call xcifc(xctype,n=np,rho=rho,grho2=grho2,ex=ex,ec=ec,vx=vx,vc=vc, &
         dxdgr2=dxdgr2,dcdgr2=dcdgr2)
        call ggamt_2b(is,g2rho,gvrho,vx,vc,dxdgr2,dcdgr2)
      else if (xcgrad.eq.3) then
        call ggamt_2a(ias,g2rho,gvrho,grho2)
        call xcifc(xctype,n=np,c_tb09=c_tb09,rho=rho,g2rho=g2rho,grho2=grho2, &
         tau=tau,vx=vx,vc=vc)
        ex(1:np)=0.d0; ec(1:np)=0.d0
      else if (xcgrad.eq.4) then
        call ggamt_2a(ias,g2rho,gvrho,grho2)
        call xcifc(xctype,n=np,rho=rho,g2rho=g2rho,grho2=grho2,tau=tau,ex=ex, &
         ec=ec,vx=vx,vc=vc,dxdgr2=dxdgr2,dcdgr2=dcdgr2,dxdg2r=dxdg2r, &
         dcdg2r=dcdg2r,wx=wx,wc=wc)
        call ggamt_2b(is,g2rho,gvrho,vx,vc,dxdgr2,dcdgr2)
        wx(1:np)=wx(1:np)+wc(1:np)
        call rfsht(nr,nri,wx,wxcmt(:,ias))
      end if
! hybrid functionals
      if (hybrid) then
        t1=1.d0-hybridc
! scale exchange part of energy
        ex(1:np)=t1*ex(1:np)
! scale exchange part of potential
        vxc(1:np)=t1*vx(1:np)+vc(1:np)
      else
        vxc(1:np)=vx(1:np)+vc(1:np)
      end if
    end if
! convert exchange and correlation energy densities to spherical harmonics
    call rfsht(nr,nri,ex,exmt(:,ias))
    call rfsht(nr,nri,ec,ecmt(:,ias))
! convert exchange-correlation potential to spherical harmonics
    call rfsht(nr,nri,vxc,vxcmt(:,ias))
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
    if (xcgrad.eq.0) then
! LSDA
      do ir=1,ngtot
        t0=rhoir(ir)
        t1=sqrt(magir(ir,1)**2+magir(ir,2)**2+magir(ir,3)**2)*ssxc
        rhoup(ir)=0.5d0*(t0+t1)
        rhodn(ir)=0.5d0*(t0-t1)
      end do
    else
! functionals which require gradients
      do ir=1,ngtot
        t0=rhoir(ir)
        t1=sqrt(magir(ir,1)**2+magir(ir,2)**2+magir(ir,3)**2+dncgga)*ssxc
        rhoup(ir)=0.5d0*(t0+t1)
        rhodn(ir)=0.5d0*(t0-t1)
      end do
    end if
  else
! collinear
    do ir=1,ngtot
      t0=rhoir(ir)
      t1=magir(ir,1)*ssxc
      rhoup(ir)=0.5d0*(t0+t1)
      rhodn(ir)=0.5d0*(t0-t1)
    end do
  end if
  if (xcgrad.le.0) then
    call xcifc(xctype,n=ngtot,rhoup=rhoup,rhodn=rhodn,ex=exir,ec=ecir, &
     vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
  else if (xcgrad.eq.1) then
    call ggair_sp_1(rhoup,rhodn,grho,gup,gdn,g2up,g2dn,g3rho,g3up,g3dn)
    call xcifc(xctype,n=ngtot,rhoup=rhoup,rhodn=rhodn,grho=grho,gup=gup, &
     gdn=gdn,g2up=g2up,g2dn=g2dn,g3rho=g3rho,g3up=g3up,g3dn=g3dn,ex=exir, &
     ec=ecir,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
  else if (xcgrad.eq.2) then
    call ggair_sp_2a(rhoup,rhodn,g2up,g2dn,gvup,gvdn,gup2,gdn2,gupdn)
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
  else if (xcgrad.eq.4) then
    call ggair_sp_2a(rhoup,rhodn,g2up,g2dn,gvup,gvdn,gup2,gdn2,gupdn)
    call xcifc(xctype,n=ngtot,rhoup=rhoup,rhodn=rhodn,g2up=g2up,g2dn=g2dn, &
     gup2=gup2,gdn2=gdn2,gupdn=gupdn,tauup=tauir(:,1),taudn=tauir(:,2), &
     ex=exir,ec=ecir,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn,dxdgu2=dxdgu2, &
     dxdgd2=dxdgd2,dxdgud=dxdgud,dcdgu2=dcdgu2,dcdgd2=dcdgd2,dcdgud=dcdgud, &
     dxdg2u=dxdg2u,dxdg2d=dxdg2d,dcdg2u=dcdg2u,dcdg2d=dcdg2d,wxup=wxup, &
     wxdn=wxdn,wcup=wcup,wcdn=wcdn)
    call ggair_sp_2b(g2up,g2dn,gvup,gvdn,vxup,vxdn,vcup,vcdn,dxdgu2,dxdgd2, &
     dxdgud,dcdgu2,dcdgd2,dcdgud)
    wxcir(1:ngtot)=0.5d0*(wxup(1:ngtot)+wxdn(1:ngtot) &
                         +wcup(1:ngtot)+wcdn(1:ngtot))
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
! scale field if required
  if (tssxc) bxcir(:,1:ndmag)=bxcir(:,1:ndmag)*ssxc
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
     vc=vc,dxdgr2=dxdgr2,dcdgr2=dcdgr2)
    call ggair_2b(g2rho,gvrho,vx,vc,dxdgr2,dcdgr2)
  else if (xcgrad.eq.3) then
    call ggair_2a(g2rho,gvrho,grho2)
    call xcifc(xctype,n=ngtot,c_tb09=c_tb09,rho=rhoir,g2rho=g2rho,grho2=grho2, &
     tau=tauir,vx=vx,vc=vc)
    exir(1:ngtot)=0.d0; ecir(1:ngtot)=0.d0
  else if (xcgrad.eq.4) then
    call ggair_2a(g2rho,gvrho,grho2)
    call xcifc(xctype,n=ngtot,rho=rhoir,g2rho=g2rho,grho2=grho2,tau=tauir, &
     ex=exir,ec=ecir,vx=vx,vc=vc,dxdgr2=dxdgr2,dcdgr2=dcdgr2,dxdg2r=dxdg2r, &
     dcdg2r=dcdg2r,wx=wx,wc=wc)
    call ggair_2b(g2rho,gvrho,vx,vc,dxdgr2,dcdgr2)
    wxcir(1:ngtot)=wx(1:ngtot)+wc(1:ngtot)
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
call symrf(nrmt,nrmti,npmt,npmtmax,vxcmt,vxcir)
if (spinpol) then
! symmetrise the exchange-correlation effective field
  call symrvf(.true.,ncmag,nrmt,nrmti,npmt,npmtmax,bxcmt,bxcir)
! remove the source contribution if required
  if (nosource) call projsbf
end if
deallocate(rho,ex,ec,vxc)
if ((xcgrad.eq.3).or.(xcgrad.eq.4)) deallocate(tau)
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
    deallocate(dxdgr2,dcdgr2)
  else if (xcgrad.eq.3) then
    deallocate(g2rho,gvrho,grho2)
  else if (xcgrad.eq.4) then
    deallocate(g2rho,gvrho,grho2)
    deallocate(dxdgr2,dcdgr2,dxdg2r,dcdg2r)
    deallocate(wx,wc)
  end if
end if
return
end subroutine
!EOC

