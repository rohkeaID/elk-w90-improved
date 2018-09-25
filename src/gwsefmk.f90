
! Copyright (C) 2016 A. Davydov, A. Sanna, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine gwsefmk(ikp,vmt,vir,bmt,bir,swfm)
use modmain
use modgw
use modomp
implicit none
! arguments
integer, intent(in) :: ikp
real(8), intent(in) :: vmt(npcmtmax,natmtot),vir(ngtot)
real(8), intent(in) :: bmt(npcmtmax,natmtot,ndmag),bir(ngtot,ndmag)
complex(8), intent(out) :: swfm(nstsv,nstsv,0:nwfm)
! local variables
integer ik,jk,ist1,ist2,ist3
integer iv(3),iq,ig,jg,i
integer iw,jw,it,nthd
real(8) vl(3),vc(3),tp(2)
real(8) t1,t2
complex(8) z1,z2
! automatic arrays
integer idx(nstsv)
! allocatable arrays
integer(8), allocatable :: lock(:)
real(8), allocatable :: vgqc(:,:),gqc(:),gclgq(:)
real(8), allocatable :: jlgqr(:,:,:),jlgqrmt(:,:,:)
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
allocate(vgqc(3,ngvc),gqc(ngvc),gclgq(ngvc))
allocate(jlgqr(njcmax,nspecies,ngrf),jlgqrmt(0:lnpsd,ngvc,nspecies))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(ylmgq(lmmaxo,ngvc),sfacgq(ngvc,natmtot))
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
allocate(wfmt1(npcmtmax,natmtot,nspinor,nstsv),wfir1(ngtc,nspinor,nstsv))
allocate(wfmt2(npcmtmax,natmtot,nspinor,nstsv),wfir2(ngtc,nspinor,nstsv))
allocate(zrhomt(npcmtmax,natmtot,nstsv),zrhoir(ngtc,nstsv))
allocate(zrho(nstsv,ngrf,nstsv),epsi(ngrf,ngrf,nwrf))
allocate(v(nstsv,nstsv),stau(nstsv,nstsv,nwgw),gs(nwgw,nstsv))
! initialise the OpenMP locks
allocate(lock(nwgw))
do it=1,nwgw
  call omp_init_lock(lock(it))
end do
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
call genwfsv(.false.,.true.,nstsv,idx,ngdc,igfc,ngk(1,ikp),igkig(:,1,ikp), &
 apwalm,evecfv,evecsv,wfmt1,ngtc,wfir1)
! local -V_xc and -B_xc matrix elements
if (spinpol) then
  call genvbmatk(vmt,vir,bmt,bir,ngk(1,ikp),igkig(:,1,ikp),wfmt1,ngtc,wfir1,v)
else
  call genvmatk(vmt,vir,ngk(1,ikp),igkig(:,1,ikp),wfmt1,ngtc,wfir1,v)
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
  do ig=1,ngvc
! determine G+q-vectors
    vgqc(:,ig)=vgc(:,ig)+vc(:)
! G+q-vector length and (theta, phi) coordinates
    call sphcrd(vgqc(:,ig),gqc(ig),tp)
! spherical harmonics for G+q-vector
    call genylm(lmaxo,tp,ylmgq(:,ig))
  end do
! structure factors for G+q
  call gensfacgp(ngvc,vgqc,ngvc,sfacgq)
! generate the regularised Coulomb Green's function in G+q-space
  call gengclgq(.true.,iq,ngvc,gqc,gclgq)
! compute the required spherical Bessel functions
  call genjlgprmt(lnpsd,ngvc,gqc,ngvc,jlgqrmt)
  call genjlgqr(gqc,jlgqr)
! find the matching coefficients
  call match(ngk(1,ik),gkc(:,1,ik),tpgkc(:,:,1,ik),sfacgk(:,:,1,ik),apwalm)
! get the eigenvectors from file for non-reduced k-point
  call getevecfv(filext,0,vkl(:,ik),vgkl(:,:,1,ik),evecfv)
  call getevecsv(filext,0,vkl(:,ik),evecsv)
! calculate the wavefunctions for all states
  call genwfsv(.false.,.false.,nstsv,idx,ngdc,igfc,ngk(1,ik),igkig(:,1,ik), &
   apwalm,evecfv,evecsv,wfmt2,ngtc,wfir2)
! determine the complex densities and Fourier transform to G+q-space
  do ist3=1,nstsv
    call omp_hold(nstsv,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zfgq) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
    do ist1=1,nstsv
      allocate(zfgq(ngrf))
      call genzrho(.true.,.true.,ngtc,wfmt2(:,:,:,ist3),wfir2(:,:,ist3), &
       wfmt1(:,:,:,ist1),wfir1(:,:,ist1),zrhomt(:,:,ist1),zrhoir(:,ist1))
      call zftzf(ngrf,jlgqr,ylmgq,ngvc,sfacgq,zrhomt(:,:,ist1),zrhoir(:,ist1), &
       zfgq)
      zrho(ist1,:,ist3)=conjg(zfgq(:))
      deallocate(zfgq)
    end do
!$OMP END DO
!$OMP END PARALLEL
    call omp_free(nthd)
!--------------------------------------!
!     valence Fock matrix elements     !
!--------------------------------------!
    t1=wqptnr*occsv(ist3,jk)/occmax
    if (abs(t1).lt.epsocc) cycle
    call omp_hold(nstsv,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zvclmt,zvclir) &
!$OMP PRIVATE(ist1,z1) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
    do ist2=1,nstsv
      allocate(zvclmt(npcmtmax,natmtot),zvclir(ngtc))
! calculate the Coulomb potential
      call genzvclmt(nrcmt,nrcmti,nrcmtmax,rcmt,npcmtmax,zrhomt(:,:,ist2), &
       zvclmt)
      call zpotcoul(nrcmt,nrcmti,npcmt,npcmti,nrcmtmax,rcmt,ngdc,igfc,ngvc, &
       gqc,gclgq,ngvc,jlgqrmt,ylmgq,sfacgq,zrhoir(:,ist2),npcmtmax,zvclmt, &
       zvclir)
      do ist1=1,ist2
        z1=zfinp(zrhomt(:,:,ist1),zrhoir(:,ist1),zvclmt,zvclir)
        v(ist1,ist2)=v(ist1,ist2)-t1*z1
      end do
      deallocate(zvclmt,zvclir)
    end do
!$OMP END DO
!$OMP END PARALLEL
    call omp_free(nthd)
  end do
!-------------------------------------!
!     correlation matrix elements     !
!-------------------------------------!
! symmetrise the Coulomb Green's function
  gclgq(1:ngrf)=sqrt(gclgq(1:ngrf))
! generate G_s in tau-space
  call omp_hold(nstsv,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(t1,iw) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
  do ist1=1,nstsv
    t1=efermi-evalsv(ist1,jk)
    gs(:,ist1)=0.d0
    do iw=-nwfm,nwfm,2
      gs(iwfft(iw),ist1)=1.d0/cmplx(t1,wgw(iw),8)
    end do
    call zfftifc(1,nwgw,1,gs(:,ist1))
  end do
!$OMP END DO
!$OMP END PARALLEL
  call omp_free(nthd)
! get RPA inverse epsilon from file
  call getcfgq('EPSINV.OUT',vl,ngrf,nwrf,epsi)
  call omp_hold(ngrf,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wc,zv,t1,t2,ig,iw,jw) &
!$OMP PRIVATE(it,ist2,ist3,z1,z2) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
  do jg=1,ngrf
! if epsi is exactly zero then there is no entry for this particular G'+q-vector
! so we cycle to the next
    if (abs(epsi(jg,jg,1)).eq.0.d0) cycle
    allocate(wc(nwgw,ngrf),zv(nstsv))
! subtract one from inverse epsilon to leave just the correlation part
    epsi(jg,jg,:)=epsi(jg,jg,:)-1.d0
! compute the correlation part of the screened interaction W_c
    t1=gclgq(jg)
    do ig=1,ngrf
      t2=t1*gclgq(ig)
      wc(:,ig)=0.d0
      do iw=-nwbs,nwbs,2
        jw=(iw+nwbs)/2+1
        wc(iwfft(iw),ig)=t2*epsi(ig,jg,jw)
      end do
! Fourier transform W_c to tau-space
      call zfftifc(1,nwgw,1,wc(:,ig))
    end do
    do it=1,nwgw
      do ist3=1,nstsv
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
        call omp_set_lock(lock(it))
        if (gwdiag.eq.0) then
! compute the full self-energy matrix
          do ist2=1,nstsv
            z2=conjg(zrho(ist2,jg,ist3))
            call zaxpy(nstsv,z2,zv,1,stau(:,ist2,it),1)
          end do
        else
! compute only the diagonal elements of the self-energy
          do ist2=1,nstsv
            z2=conjg(zrho(ist2,jg,ist3))
            stau(ist2,ist2,it)=stau(ist2,ist2,it)+z2*zv(ist2)
          end do
        end if
        call omp_unset_lock(lock(it))
      end do
    end do
    deallocate(wc,zv)
  end do
!$OMP END DO
!$OMP END PARALLEL
  call omp_free(nthd)
10 continue
! end loop over k-points
end do
! destroy the OpenMP locks
do it=1,nwgw
  call omp_destroy_lock(lock(it))
end do
deallocate(lock)
! Fourier transform the self-energy to frequency space, multiply by GW diagram
! prefactor and store in output array
t1=-wqptnr*omega*kboltz*tempk
call omp_hold(nstsv,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zv,ist1,iw,jw) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ist2=1,nstsv
  allocate(zv(nwgw))
  do ist1=1,nstsv
    zv(1:nwgw)=stau(ist1,ist2,1:nwgw)
    call zfftifc(1,nwgw,-1,zv)
    do iw=-nwfm,nwfm,2
      jw=(iw+nwfm)/2
      swfm(ist1,ist2,jw)=t1*zv(iwfft(iw))
    end do
  end do
  deallocate(zv)
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
! add the local potential and Fock matrix elements to the self-energy for each
! Matsubara frequency
call omp_hold(nwfm+1,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ist1,ist2) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
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
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
deallocate(vgqc,gqc,gclgq,jlgqr,jlgqrmt)
deallocate(ylmgq,sfacgq,apwalm,evecfv,evecsv)
deallocate(wfmt1,wfir1,wfmt2,wfir2)
deallocate(zrhomt,zrhoir,zrho)
deallocate(epsi,v,gs)
return
end subroutine

