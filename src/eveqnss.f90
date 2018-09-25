
! Copyright (C) 2006 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine eveqnss(ngp,igpig,apwalm,evalfv,evecfv,evalsvp,evecsv)
use modmain
use moddftu
use modomp
implicit none
! arguments
integer, intent(in) :: ngp(nspnfv),igpig(ngkmax,nspnfv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
real(8), intent(in) :: evalfv(nstfv,nspnfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv)
real(8), intent(out) :: evalsvp(nstsv)
complex(8), intent(out) :: evecsv(nstsv,nstsv)
! local variables
integer ist,jst,ispn,jspn
integer is,ia,ias,i,j,k
integer nrc,nrci,nrco
integer l,lm,nm,npc,npci
integer igp,ld,nthd
real(8) t1
real(8) ts0,ts1
complex(8) zq
! allocatable arrays
complex(8), allocatable :: wfmt1(:,:,:),wfmt2(:,:),wfmt3(:),wfmt4(:,:)
complex(8), allocatable :: wfir1(:,:),wfir2(:),wfgp(:,:)
! external functions
complex(8) zdotc,zfmtinp
external zdotc,zfmtinp
if (.not.spinpol) then
  write(*,*)
  write(*,'("Error(eveqnss): spin-unpolarised calculation")')
  write(*,*)
  stop
end if
call timesec(ts0)
ld=lmmaxdm*nspinor
! zero the second-variational Hamiltonian (stored in the eigenvector array)
evecsv(:,:)=0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
allocate(wfmt1(npcmtmax,nstfv,nspnfv))
allocate(wfmt2(npcmtmax,nspnfv))
allocate(wfmt3(npcmtmax),wfmt4(npcmtmax,3))
do ias=1,natmtot
  is=idxis(ias)
  ia=idxia(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  nrco=nrc-nrci
  npc=npcmt(is)
  npci=npcmti(is)
! de-phasing factor (FC, FB & LN)
  t1=-0.5d0*dot_product(vqcss(:),atposc(:,ia,is))
  zq=cmplx(cos(t1),sin(t1),8)
! compute the first-variational wavefunctions
  do ispn=1,nspnfv
    if (ispn.eq.2) zq=conjg(zq)
    do ist=1,nstfv
      call wavefmt(lradstp,ias,ngp(ispn),apwalm(:,:,:,ias,ispn), &
       evecfv(:,ist,ispn),wfmt1(:,ist,ispn))
! de-phase if required
      if (ssdph) wfmt1(1:npc,ist,ispn)=zq*wfmt1(1:npc,ist,ispn)
    end do
  end do
  do jst=1,nstfv
! convert wavefunction to spherical coordinates
    do ispn=1,nspnfv
      call zbsht(nrc,nrci,wfmt1(:,jst,ispn),wfmt2(:,ispn))
    end do
! apply effective magnetic field and convert to spherical harmonics
    wfmt3(1:npc)=bsmt(1:npc,ias,3)*wfmt2(1:npc,1)
    call zfsht(nrc,nrci,wfmt3,wfmt4(:,1))
    wfmt3(1:npc)=-bsmt(1:npc,ias,3)*wfmt2(1:npc,2)
    call zfsht(nrc,nrci,wfmt3,wfmt4(:,2))
    wfmt3(1:npc)=cmplx(bsmt(1:npc,ias,1),-bsmt(1:npc,ias,2),8)*wfmt2(1:npc,2)
    call zfsht(nrc,nrci,wfmt3,wfmt4(:,3))
! apply muffin-tin potential matrix if required
    if (tvmatmt) then
      do l=0,lmaxdm
        if (tvmmt(l,ias)) then
          nm=2*l+1
          lm=idxlm(l,-l)
          do k=1,3
            if (k.eq.1) then
              ispn=1
              jspn=1
            else if (k.eq.2) then
              ispn=2
              jspn=2
            else
              ispn=1
              jspn=2
            end if
            if (l.le.lmaxi) then
              call zgemm('N','N',nm,nrci,nm,zone,vmatmt(lm,ispn,lm,jspn,ias), &
               ld,wfmt1(lm,jst,jspn),lmmaxi,zone,wfmt4(lm,k),lmmaxi)
            end if
            i=npci+lm
            call zgemm('N','N',nm,nrco,nm,zone,vmatmt(lm,ispn,lm,jspn,ias),ld, &
             wfmt1(i,jst,jspn),lmmaxo,zone,wfmt4(i,k),lmmaxo)
          end do
        end if
      end do
    end if
! second-variational Hamiltonian matrix
    call omp_hold(nstfv,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(i,j,k,ispn) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
    do ist=1,nstfv
      do k=1,3
        if (k.eq.1) then
          ispn=1
          i=ist
          j=jst
        else if (k.eq.2) then
          ispn=2
          i=ist+nstfv
          j=jst+nstfv
        else
          ispn=1
          i=ist
          j=jst+nstfv
        end if
        if (i.le.j) then
          evecsv(i,j)=evecsv(i,j)+zfmtinp(nrc,nrci,rcmt(:,is),r2cmt(:,is), &
           wfmt1(:,ist,ispn),wfmt4(:,k))
        end if
      end do
    end do
!$OMP END DO
!$OMP END PARALLEL
    call omp_free(nthd)
  end do
end do
deallocate(wfmt1,wfmt2,wfmt3,wfmt4)
!---------------------------!
!     interstitial part     !
!---------------------------!
allocate(wfir1(ngtot,nspnfv),wfir2(ngtot),wfgp(ngkmax,3))
do jst=1,nstfv
  do ispn=1,nspnfv
    wfir1(:,ispn)=0.d0
    do igp=1,ngp(ispn)
      wfir1(igfft(igpig(igp,ispn)),ispn)=evecfv(igp,jst,ispn)
    end do
! Fourier transform wavefunction to real-space
    call zfftifc(3,ngridg,1,wfir1(:,ispn))
  end do
! multiply with magnetic field and transform to G-space
  wfir2(:)=bsir(:,3)*wfir1(:,1)
  call zfftifc(3,ngridg,-1,wfir2)
  do igp=1,ngp(1)
    wfgp(igp,1)=wfir2(igfft(igpig(igp,1)))
  end do
  wfir2(:)=-bsir(:,3)*wfir1(:,2)
  call zfftifc(3,ngridg,-1,wfir2)
  do igp=1,ngp(2)
    wfgp(igp,2)=wfir2(igfft(igpig(igp,2)))
  end do
  wfir2(:)=cmplx(bsir(:,1),-bsir(:,2),8)*wfir1(:,2)
  call zfftifc(3,ngridg,-1,wfir2)
  do igp=1,ngp(1)
    wfgp(igp,3)=wfir2(igfft(igpig(igp,1)))
  end do
! add to Hamiltonian matrix
  call omp_hold(nstfv,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(i,j,k,ispn) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
  do ist=1,nstfv
    do k=1,3
      if (k.eq.1) then
        ispn=1
        i=ist
        j=jst
      else if (k.eq.2) then
        ispn=2
        i=ist+nstfv
        j=jst+nstfv
      else
        ispn=1
        i=ist
        j=jst+nstfv
      end if
      if (i.le.j) then
        evecsv(i,j)=evecsv(i,j)+zdotc(ngp(ispn),evecfv(:,ist,ispn),1, &
         wfgp(:,k),1)
      end if
    end do
  end do
!$OMP END DO
!$OMP END PARALLEL
  call omp_free(nthd)
end do
deallocate(wfir1,wfir2,wfgp)
! add the diagonal first-variational part
i=0
do ispn=1,nspinor
  do ist=1,nstfv
    i=i+1
    evecsv(i,i)=evecsv(i,i)+evalfv(ist,ispn)
  end do
end do
! diagonalise the second-variational Hamiltonian
call eveqnz(nstsv,nstsv,evecsv,evalsvp)
call timesec(ts1)
!$OMP CRITICAL(eveqnss_)
timesv=timesv+ts1-ts0
!$OMP END CRITICAL(eveqnss_)
return
end subroutine

