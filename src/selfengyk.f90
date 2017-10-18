
! Copyright (C) 2016 A. Davydov, A. Sanna, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine selfengyk(ikp,vmt,vir,bmt,bir,swfm)
use modmain
use modgw
implicit none
! arguments
integer, intent(in) :: ikp
real(8), intent(in) :: vmt(npcmtmax,natmtot),vir(ngtot)
real(8), intent(in) :: bmt(npcmtmax,natmtot,ndmag),bir(ngtot,ndmag)
complex(8), intent(out) :: swfm(nstsv,nstsv,0:nwfm)
! local variables
integer ik,jk,ist1,ist2,ist3
integer iv(3),iq,igq0,i
integer ig,jg,iw,jw,it
real(8) vl(3),vc(3),tp(2)
real(8) cfq,t1,t2
complex(8) zgq01,zgq02,z1,z2
! automatic arrays
integer idx(nstsv)
real(8) vcl(ngrf)
complex(8) sfacgq0(natmtot)
! allocatable arrays
real(8), allocatable :: vgqc(:,:),gqc(:)
real(8), allocatable :: jlgqr(:,:,:),jlgqrmt(:,:,:),jlgq0r(:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:),evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:)
complex(8), allocatable :: wfmt1(:,:,:,:),wfir1(:,:,:)
complex(8), allocatable :: wfmt2(:,:,:,:),wfir2(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:),zrhoir(:,:)
complex(8), allocatable :: zvclmt(:,:),zvclir(:)
complex(8), allocatable :: zfgq(:),zrho(:,:,:),epsi(:,:,:)
complex(8), allocatable :: v(:,:),stau(:,:,:)
complex(8), allocatable :: gs(:,:),wc(:,:),zv(:)
! external functions
complex(8) zfinp
external zfinp
! allocate local arrays
allocate(vgqc(3,ng2gk),gqc(ng2gk))
allocate(jlgqr(njcmax,nspecies,ngrf),jlgqrmt(0:lnpsd,ng2gk,nspecies))
allocate(jlgq0r(0:lmaxo,nrcmtmax,nspecies))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(ylmgq(lmmaxo,ng2gk),sfacgq(ng2gk,natmtot))
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
allocate(wfmt1(npcmtmax,natmtot,nspinor,nstsv),wfir1(ngtot,nspinor,nstsv))
allocate(wfmt2(npcmtmax,natmtot,nspinor,nstsv),wfir2(ngtot,nspinor,nstsv))
allocate(zrhomt(npcmtmax,natmtot,nstsv),zrhoir(ngtot,nstsv))
allocate(zrho(nstsv,ngrf,nstsv),epsi(ngrf,ngrf,nwrf))
allocate(v(nstsv,nstsv),stau(nstsv,nstsv,nwgw))
allocate(gs(nwgw,nstsv),wc(nwgw,ngrf),zv(max(nstsv,nwgw)))
! coefficient of long-range term
cfq=0.5d0*(omega/pi)**2
! get the eigenvectors from file for input reduced k-point
call getevecfv(filext,ikp,vkl(:,ikp),vgkl(:,:,:,ikp),evecfv)
call getevecsv(filext,ikp,vkl(:,ikp),evecsv)
! find the matching coefficients
call match(ngk(1,ikp),gkc(:,1,ikp),tpgkc(:,:,1,ikp),sfacgk(:,:,1,ikp),apwalm)
! index to all states
do ist1=1,nstsv
  idx(ist1)=ist1
end do
! calculate the wavefunctions for all states of the input k-point
call genwfsv(.false.,.true.,nstsv,idx,ngk(1,ikp),igkig(:,1,ikp),apwalm, &
 evecfv,evecsv,wfmt1,ngtot,wfir1)
! local -V_xc and -B_xc matrix elements
if (spinpol) then
  call genvbmatk(vmt,vir,bmt,bir,ngk(1,ikp),igkig(:,1,ikp),wfmt1,ngtot,wfir1,v)
else
  call genvmatk(vmt,vir,ngk(1,ikp),igkig(:,1,ikp),wfmt1,ngtot,wfir1,v)
end if
! Fourier transform wavefunctions to real-space
call zftwfir(ngk(1,ikp),igkig(:,1,ikp),wfir1)
! add the core Fock matrix elements
call vclcore(wfmt1,v)
! zero the self-energy matrix elements in tau-space
stau(:,:,:)=0.d0
! loop over non-reduced k-point set
do ik=1,nkptnr
! find the equivalent reduced k-point
  jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
! determine the q-vector
  iv(:)=ivk(:,ikp)-ivk(:,ik)
  iv(:)=modulo(iv(:),ngridk(:))
! check if the q-point is in user-defined set
  iv(:)=iv(:)*ngridq(:)
  do i=1,3
    if (modulo(iv(i),ngridk(i)).ne.0) goto 10
  end do
  iv(:)=iv(:)/ngridk(:)
  iq=iqmap(iv(1),iv(2),iv(3))
  vl(:)=vkl(:,ikp)-vkl(:,ik)
  vc(:)=vkc(:,ikp)-vkc(:,ik)
  do ig=1,ng2gk
! determine G+q-vectors
    vgqc(:,ig)=vgc(:,ig)+vc(:)
! G+q-vector length and (theta, phi) coordinates
    call sphcrd(vgqc(:,ig),gqc(ig),tp)
! spherical harmonics for G+q-vector
    call genylm(lmaxo,tp,ylmgq(:,ig))
  end do
! structure factors for G+q
  call gensfacgp(ng2gk,vgqc,ng2gk,sfacgq)
! find the shortest G+q-vector
  call findigp0(ng2gk,gqc,igq0)
  sfacgq0(:)=sfacgq(igq0,:)
! compute the required spherical Bessel functions
  call genjlgprmt(lnpsd,ng2gk,gqc,ng2gk,jlgqrmt)
  call genjlgq0r(gqc(igq0),jlgq0r)
  call genjlgqr(gqc,jlgqr)
! find the matching coefficients
  call match(ngk(1,ik),gkc(:,1,ik),tpgkc(:,:,1,ik),sfacgk(:,:,1,ik),apwalm)
! get the eigenvectors from file for non-reduced k-point
  call getevecfv(filext,0,vkl(:,ik),vgkl(:,:,1,ik),evecfv)
  call getevecsv(filext,0,vkl(:,ik),evecsv)
! calculate the wavefunctions for all states
  call genwfsv(.false.,.false.,nstsv,idx,ngk(1,ik),igkig(:,1,ik),apwalm, &
   evecfv,evecsv,wfmt2,ngtot,wfir2)
! determine the complex densities and transform to G+q-space
  do ist3=1,nstsv
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(zfgq)
!$OMP DO
    do ist1=1,nstsv
      allocate(zfgq(ngrf))
      call genzrho(.true.,.true.,wfmt2(:,:,:,ist3),wfir2(:,:,ist3), &
       wfmt1(:,:,:,ist1),wfir1(:,:,ist1),zrhomt(:,:,ist1),zrhoir(:,ist1))
      call zftzf(ngrf,jlgqr,ylmgq,ng2gk,sfacgq,zrhomt(:,:,ist1), &
       zrhoir(:,ist1),zfgq)
      zrho(ist1,:,ist3)=conjg(zfgq(:))
      deallocate(zfgq)
    end do
!$OMP END DO
!$OMP END PARALLEL
!--------------------------------------!
!     valence Fock matrix elements     !
!--------------------------------------!
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zvclmt,zvclir,zgq01,zgq02) &
!$OMP PRIVATE(ist1,z1,z2,t1)
!$OMP DO
    do ist2=1,nstsv
      allocate(zvclmt(npcmtmax,natmtot),zvclir(ngtot))
! calculate the Coulomb potential
      call genzvclmt(nrcmt,nrcmti,nrcmtmax,rcmt,npcmtmax,zrhomt(:,:,ist2), &
       zvclmt)
      call zpotcoul(nrcmt,nrcmti,npcmt,npcmti,nrcmtmax,rcmt,ng2gk,igq0,gqc, &
       ng2gk,jlgqrmt,ylmgq,sfacgq,zrhoir(:,ist2),npcmtmax,zvclmt,zvclir,zgq02)
      do ist1=1,ist2
        z1=zfinp(zrhomt(:,:,ist1),zrhoir(:,ist1),zvclmt,zvclir)
! compute the density coefficient of the smallest G+q-vector
        call zftgq0(jlgq0r,ylmgq(:,igq0),sfacgq0,zrhomt(:,:,ist1), &
         zrhoir(:,ist1),zgq01)
        z2=cfq*wiq2(iq)*conjg(zgq01)*zgq02
        t1=occsv(ist3,jk)/occmax
        v(ist1,ist2)=v(ist1,ist2)-t1*(wqptnr*z1+z2)
      end do
      deallocate(zvclmt,zvclir)
    end do
!$OMP END DO
!$OMP END PARALLEL
  end do
!-------------------------------------!
!     correlation matrix elements     !
!-------------------------------------!
! compute the regularised Coulomb interaction
  do ig=1,ngrf
    if (ig.eq.igq0) then
! volume of small parallelepiped around q-point (see genwiq2)
      t2=omegabz*wqptnr
! average symmetrised interaction over volume
      vcl(ig)=sqrt(fourpi*wiq2(iq)/t2)
    else
! G+q-vector is outside FBZ so use symmetrised 4 pi/(G+q)^2 interaction
      vcl(ig)=sqrt(fourpi)/gqc(ig)
    end if
  end do
! generate G_s in tau-space
  do ist1=1,nstsv
    t1=efermi-evalsv(ist1,jk)
    gs(:,ist1)=0.d0
    do iw=-nwfm,nwfm,2
      gs(iwfft(iw),ist1)=1.d0/cmplx(t1,wgw(iw),8)
    end do
    call zfftifc(1,nwgw,1,gs(:,ist1))
  end do
! get RPA inverse epsilon from file
  call getcf2pt('EPSINV.OUT',vl,ngrf,nwrf,epsi)
  do jg=1,ngrf
! subtract one from inverse epsilon to leave just the correlation part
    epsi(jg,jg,:)=epsi(jg,jg,:)-1.d0
! compute the correlation part of the screened interaction W_c
    t1=vcl(jg)
    do ig=1,ngrf
      t2=t1*vcl(ig)
      wc(:,ig)=0.d0
      do iw=-nwbs,nwbs,2
        jw=(iw+nwbs)/2+1
        wc(iwfft(iw),ig)=t2*epsi(ig,jg,jw)
      end do
! Fourier transform W_c to tau-space
      call zfftifc(1,nwgw,1,wc(:,ig))
    end do
    do ist3=1,nstsv
      do it=1,nwgw
        z1=gs(it,ist3)
        zv(1:nstsv)=0.d0
        if (gwdiag.eq.2) then
! use only the diagonal elements of W_c
          z2=z1*wc(it,jg)
          call zaxpy(nstsv,z2,zrho(:,jg,ist3),1,zv,1)
        else
! use the full W_c
          do ig=1,ngrf
            z2=z1*wc(it,ig)
            call zaxpy(nstsv,z2,zrho(:,ig,ist3),1,zv,1)
          end do
        end if
        do ist2=1,nstsv
          z2=conjg(zrho(ist2,jg,ist3))
          if (gwdiag.eq.0) then
! compute the full self-energy matrix
            call zaxpy(nstsv,z2,zv,1,stau(:,ist2,it),1)
          else
! compute only the diagonal elements of the self-energy
            stau(ist2,ist2,it)=stau(ist2,ist2,it)+z2*zv(ist2)
          end if
        end do
      end do
    end do
  end do
10 continue
! end loop over k-points
end do
! Fourier transform the self-energy to frequency space, multiply by GW diagram
! prefactor and store in output array
t1=-wqptnr*omega*kboltz*tempk
do ist1=1,nstsv
  do ist2=1,nstsv
    zv(1:nwgw)=stau(ist1,ist2,1:nwgw)
    call zfftifc(1,nwgw,-1,zv)
    do iw=-nwfm,nwfm,2
      jw=(iw+nwfm)/2
      swfm(ist1,ist2,jw)=t1*zv(iwfft(iw))
    end do
  end do
end do
! add the local potential and Fock matrix elements to the self-energy for each
! Matsubara frequency
do iw=0,nwfm
  do ist2=1,nstsv
    do ist1=1,ist2
      swfm(ist1,ist2,iw)=swfm(ist1,ist2,iw)+v(ist1,ist2)
    end do
    do ist1=ist2+1,nstsv
      swfm(ist1,ist2,iw)=swfm(ist1,ist2,iw)+conjg(v(ist2,ist1))
    end do
  end do
end do
deallocate(vgqc,gqc,jlgqr,jlgqrmt,jlgq0r)
deallocate(ylmgq,sfacgq,apwalm,evecfv,evecsv)
deallocate(wfmt1,wfir1,wfmt2,wfir2)
deallocate(zrhomt,zrhoir,zrho)
deallocate(epsi,v,gs,wc,zv)
return
end subroutine

