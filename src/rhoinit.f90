
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rhoinit
! !INTERFACE:
subroutine rhoinit
! !USES:
use modmain
! !DESCRIPTION:
!   Initialises the crystal charge density. Inside the muffin-tins it is set to
!   the spherical atomic density. In the interstitial region it is taken to be
!   constant such that the total charge is correct. Requires that the atomic
!   densities have already been calculated.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer, parameter :: lmax=1
integer is,ia,ias
integer nr,nrc,nro,ir,irc
integer l,m,lm,ig,ifg
real(8) x,t1,t2
complex(8) z1,z2,z3
! allocatable arrays
real(8), allocatable :: jl(:,:),ffg(:),fr(:)
complex(8), allocatable :: zfmt(:,:),zfft(:)
! external functions
real(8) fintgt
external fintgt
allocate(zfft(ngtot))
! zero the charge density and magnetisation arrays
rhomt(:,:,:)=0.d0
rhoir(:)=0.d0
if (spinpol) then
  magmt(:,:,:,:)=0.d0
  magir(:,:)=0.d0
end if
! compute the superposition of all the atomic density tails
zfft(:)=0.d0
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ffg,fr,nr,nro,ig) &
!$OMP PRIVATE(ir,x,t1,ia,ias,ifg)
!$OMP DO
do is=1,nspecies
  allocate(ffg(ngvec),fr(nrspmax))
  nr=nrmt(is)
  nro=nrsp(is)-nrmt(is)+1
  do ig=1,ngvec
    do ir=nr,nrsp(is)
! spherical bessel function j_0(x)
      x=gc(ig)*rsp(ir,is)
      if (x.gt.1.d-8) then
        t1=sin(x)/x
      else
        t1=1.d0
      end if
      fr(ir)=t1*rhosp(ir,is)*r2sp(ir,is)
    end do
    t1=fintgt(-1,nro,rsp(nr,is),fr(nr))
    ffg(ig)=(fourpi/omega)*t1
  end do
  do ia=1,natoms(is)
    ias=idxas(ia,is)
!$OMP CRITICAL
    do ig=1,ngvec
      ifg=igfft(ig)
      zfft(ifg)=zfft(ifg)+ffg(ig)*conjg(sfacg(ig,ias))
    end do
!$OMP END CRITICAL
  end do
  deallocate(fr,ffg)
end do
!$OMP END DO
!$OMP END PARALLEL
! compute the tails in each muffin-tin
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(jl,zfmt,is,nrc,ig,ifg) &
!$OMP PRIVATE(irc,x,z1,z2,z3,lm,l,m)
!$OMP DO
do ias=1,natmtot
  allocate(jl(0:lmax,nrcmtmax))
  allocate(zfmt(lmmaxvr,nrcmtmax))
  is=idxis(ias)
  nrc=nrcmt(is)
  zfmt(:,:)=0.d0
  do ig=1,ngvec
    ifg=igfft(ig)
    do irc=1,nrc
      x=gc(ig)*rcmt(irc,is)
      call sbessel(lmax,x,jl(:,irc))
    end do
    z1=fourpi*zfft(ifg)*sfacg(ig,ias)
    lm=0
    do l=0,lmax
      z2=z1*zil(l)
      do m=-l,l
        lm=lm+1
        z3=z2*conjg(ylmg(lm,ig))
        do irc=1,nrc
          zfmt(lm,irc)=zfmt(lm,irc)+jl(l,irc)*z3
        end do
      end do
    end do
  end do
  call ztorfmt(nrc,nrcmtinr(is),1,zfmt,lradstp,rhomt(:,:,ias))
  deallocate(jl,zfmt)
end do
!$OMP END DO
!$OMP END PARALLEL
! convert the density from a coarse to a fine radial mesh
call rfmtctof(rhomt)
! add the atomic charge density and the excess charge in each muffin-tin
t1=chgexs/omega
do ias=1,natmtot
  is=idxis(ias)
  do ir=1,nrmt(is)
    t2=(t1+rhosp(ir,is))/y00
    rhomt(1,ir,ias)=rhomt(1,ir,ias)+t2
  end do
end do
! interstitial density determined from the atomic tails and excess charge
call zfftifc(3,ngridg,1,zfft)
do ir=1,ngtot
  rhoir(ir)=dble(zfft(ir))+t1
end do
deallocate(zfft)
return
end subroutine
!EOC

