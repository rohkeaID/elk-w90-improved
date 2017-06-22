
! Copyright (C) 2006 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine eveqnss(ngp,igpig,apwalm,evalfv,evecfv,evalsvp,evecsv)
use modmain
use moddftu
implicit none
! arguments
integer, intent(in) :: ngp(nspnfv),igpig(ngkmax,nspnfv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
real(8), intent(in) :: evalfv(nstfv,nspnfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv)
real(8), intent(out) :: evalsvp(nstsv)
complex(8), intent(out) :: evecsv(nstsv,nstsv)
! local variables
integer ld,ist,jst,ispn,jspn
integer is,ia,ias,nrc,nrci,iro
integer l,lm,lmi,nm,igp,ifg
integer i,j,k,lwork,info
real(8) t1
real(8) ts0,ts1
complex(8) zq
! allocatable arrays
real(8), allocatable :: rwork(:)
complex(8), allocatable :: wfmt1(:,:,:,:),wfmt2(:,:,:),wfmt3(:,:),wfmt4(:,:,:)
complex(8), allocatable :: wfir1(:,:),wfir2(:),z(:,:),work(:)
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
lmi=lmmaxinr
allocate(wfmt1(lmmaxvr,nrcmtmax,nstfv,nspnfv))
allocate(wfmt2(lmmaxvr,nrcmtmax,nspnfv))
allocate(wfmt3(lmmaxvr,nrcmtmax),wfmt4(lmmaxvr,nrcmtmax,3))
do ias=1,natmtot
  is=idxis(ias)
  ia=idxia(ias)
  nrc=nrcmt(is)
  nrci=nrcmtinr(is)
  iro=nrci+1
! de-phasing factor (FC, FB & LN)
  t1=-0.5d0*dot_product(vqcss(:),atposc(:,ia,is))
  zq=cmplx(cos(t1),sin(t1),8)
! compute the first-variational wavefunctions
  do ispn=1,nspnfv
    if (ispn.eq.2) zq=conjg(zq)
    do ist=1,nstfv
      call wavefmt(lradstp,ias,ngp(ispn),apwalm(:,:,:,:,ispn), &
       evecfv(:,ist,ispn),wfmt1(:,:,ist,ispn))
! de-phase if required
      if (ssdph) wfmt1(:,1:nrc,ist,ispn)=zq*wfmt1(:,1:nrc,ist,ispn)
    end do
  end do
  do jst=1,nstfv
! convert wavefunction to spherical coordinates
    do ispn=1,nspnfv
      call zbsht(nrc,nrci,wfmt1(:,:,jst,ispn),wfmt2(:,:,ispn))
    end do
! apply effective magnetic field and convert to spherical harmonics
    wfmt3(1:lmi,1:nrci)=bsmt(1:lmi,1:nrci,ias,3)*wfmt2(1:lmi,1:nrci,1)
    wfmt3(:,iro:nrc)=bsmt(:,iro:nrc,ias,3)*wfmt2(:,iro:nrc,1)
    call zfsht(nrc,nrci,wfmt3,wfmt4(:,:,1))
    wfmt3(1:lmi,1:nrci)=-bsmt(1:lmi,1:nrci,ias,3)*wfmt2(1:lmi,1:nrci,2)
    wfmt3(:,iro:nrc)=-bsmt(:,iro:nrc,ias,3)*wfmt2(:,iro:nrc,2)
    call zfsht(nrc,nrci,wfmt3,wfmt4(:,:,2))
    wfmt3(1:lmi,1:nrci)=cmplx(bsmt(1:lmi,1:nrci,ias,1), &
     -bsmt(1:lmi,1:nrci,ias,2),8)*wfmt2(1:lmi,1:nrci,2)
    wfmt3(:,iro:nrc)=cmplx(bsmt(:,iro:nrc,ias,1), &
     -bsmt(:,iro:nrc,ias,2),8)*wfmt2(:,iro:nrc,2)
    call zfsht(nrc,nrci,wfmt3,wfmt4(:,:,3))
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
            call zgemm('N','N',nm,nrc,nm,zone,vmatmt(lm,ispn,lm,jspn,ias),ld, &
             wfmt1(lm,1,jst,jspn),lmmaxvr,zone,wfmt4(lm,1,k),lmmaxvr)
          end do
        end if
      end do
    end if
! second-variational Hamiltonian matrix
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,ispn)
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
           wfmt1(:,:,ist,ispn),wfmt4(:,:,k))
        end if
      end do
    end do
!$OMP END DO
!$OMP END PARALLEL
  end do
end do
deallocate(wfmt1,wfmt2,wfmt3,wfmt4)
!---------------------------!
!     interstitial part     !
!---------------------------!
allocate(wfir1(ngtot,nspnfv),wfir2(ngtot),z(ngkmax,3))
do jst=1,nstfv
  do ispn=1,nspnfv
    wfir1(:,ispn)=0.d0
    do igp=1,ngp(ispn)
      ifg=igfft(igpig(igp,ispn))
      wfir1(ifg,ispn)=evecfv(igp,jst,ispn)
    end do
! Fourier transform wavefunction to real-space
    call zfftifc(3,ngridg,1,wfir1(:,ispn))
  end do
! multiply with magnetic field and transform to G-space
  wfir2(:)=bsir(:,3)*wfir1(:,1)
  call zfftifc(3,ngridg,-1,wfir2)
  do igp=1,ngp(1)
    ifg=igfft(igpig(igp,1))
    z(igp,1)=wfir2(ifg)
  end do
  wfir2(:)=-bsir(:,3)*wfir1(:,2)
  call zfftifc(3,ngridg,-1,wfir2)
  do igp=1,ngp(2)
    ifg=igfft(igpig(igp,2))
    z(igp,2)=wfir2(ifg)
  end do
  wfir2(:)=cmplx(bsir(:,1),-bsir(:,2),8)*wfir1(:,2)
  call zfftifc(3,ngridg,-1,wfir2)
  do igp=1,ngp(1)
    ifg=igfft(igpig(igp,1))
    z(igp,3)=wfir2(ifg)
  end do
! add to Hamiltonian matrix
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,ispn)
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
        evecsv(i,j)=evecsv(i,j)+zdotc(ngp(ispn),evecfv(:,ist,ispn),1,z(:,k),1)
      end if
    end do
  end do
!$OMP END DO
!$OMP END PARALLEL
end do
deallocate(wfir1,wfir2,z)
! add the diagonal first-variational part
i=0
do ispn=1,nspinor
  do ist=1,nstfv
    i=i+1
    evecsv(i,i)=evecsv(i,i)+evalfv(ist,ispn)
  end do
end do
! diagonalise the second-variational Hamiltonian
allocate(rwork(3*nstsv))
lwork=2*nstsv
allocate(work(lwork))
call zheev('V','U',nstsv,evecsv,nstsv,evalsvp,work,lwork,rwork,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(eveqnss): diagonalisation of the second-variational &
   &Hamiltonian failed")')
  write(*,'(" ZHEEV returned INFO = ",I8)') info
  write(*,*)
  stop
end if
deallocate(rwork,work)
call timesec(ts1)
!$OMP CRITICAL
timesv=timesv+ts1-ts0
!$OMP END CRITICAL
return
end subroutine

