
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Ambrosch-Draxl,
! F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine eveqnsv(ngp,igpig,vgpc,apwalm,evalfv,evecfv,evalsvp,evecsv)
use modmain
use moddftu
use modomp
implicit none
! arguments
integer, intent(in) :: ngp,igpig(ngkmax)
real(8), intent(in) :: vgpc(3,ngkmax)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
real(8), intent(in) :: evalfv(nstfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
real(8), intent(out) :: evalsvp(nstsv)
complex(8), intent(out) :: evecsv(nstsv,nstsv)
! local variables
logical socz
integer nsc,nsd,ld,ist,jst
integer ispn,jspn,is,ias
integer nrc,nrci,nrco,irco,irc
integer l,lm,nm,npc,npci
integer igp,i,j,k,nthd
real(8) ca,t1
real(8) ts0,ts1
complex(8) zsum,z1
! automatic arrays
complex(8) zlflm(lmmaxo,3)
! allocatable arrays
complex(8), allocatable :: wfmt1(:,:),wfmt2(:),wfmt3(:),wfmt4(:,:)
complex(8), allocatable :: gwfmt(:,:,:),gzfmt(:,:)
complex(8), allocatable :: wfir1(:),wfir2(:),wfgp(:,:),gwfgp(:,:)
! external functions
complex(8) zdotc,zfmtinp
external zdotc,zfmtinp
! no calculation of second-variational eigenvectors
if (.not.tevecsv) then
  do i=1,nstsv
    evalsvp(i)=evalfv(i)
  end do
  evecsv(:,:)=0.d0
  do i=1,nstsv
    evecsv(i,i)=1.d0
  end do
  return
end if
call timesec(ts0)
! coupling constant of the external A-field (1/c)
ca=1.d0/solsc
! number of spin combinations after application of Hamiltonian
if (spinpol) then
  if (ncmag.or.spinorb) then
    nsc=3
  else
    nsc=2
  end if
  nsd=2
else
  nsc=1
  nsd=1
end if
! special case of spin-orbit coupling and collinear magnetism
if (spinorb.and.cmagz) then
  socz=.true.
else
  socz=.false.
end if
ld=lmmaxdm*nspinor
! zero the second-variational Hamiltonian (stored in the eigenvector array)
evecsv(:,:)=0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
allocate(wfmt1(npcmtmax,nstfv),wfmt2(npcmtmax))
allocate(wfmt3(npcmtmax),wfmt4(npcmtmax,nsc))
if (tafield.or.(xcgrad.eq.4)) allocate(gzfmt(npcmtmax,3))
if (xcgrad.eq.4) allocate(gwfmt(npcmtmax,3,nstfv))
! begin loop over atoms
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  nrco=nrc-nrci
  irco=nrci+1
  npc=npcmt(is)
  npci=npcmti(is)
! compute the first-variational wavefunctions
  do ist=1,nstfv
    call wavefmt(lradstp,ias,ngp,apwalm(:,:,:,ias),evecfv(:,ist),wfmt1(:,ist))
  end do
! begin loop over states
  do jst=1,nstfv
    if (spinpol) then
! convert wavefunction to spherical coordinates
      call zbsht(nrc,nrci,wfmt1(:,jst),wfmt2)
! apply Kohn-Sham effective magnetic field
      wfmt3(1:npc)=bsmt(1:npc,ias,ndmag)*wfmt2(1:npc)
! convert to spherical harmonics and store in wfmt4
      call zfsht(nrc,nrci,wfmt3,wfmt4(:,1))
      wfmt4(1:npc,2)=-wfmt4(1:npc,1)
! non-collinear magnetic field
      if (ncmag) then
        wfmt3(1:npc)=cmplx(bsmt(1:npc,ias,1),-bsmt(1:npc,ias,2),8)*wfmt2(1:npc)
        call zfsht(nrc,nrci,wfmt3,wfmt4(:,3))
      end if
      if (socz) wfmt4(1:npc,3)=0.d0
! apply spin-orbit coupling if required
      if (spinorb) then
! inner part of muffin-tin
        i=1
        do irc=1,nrci
          call lopzflm(lmaxi,wfmt1(i,jst),lmmaxo,zlflm)
          t1=socfr(irc,ias)
          do j=1,lmmaxi
            wfmt4(i,1)=wfmt4(i,1)+t1*zlflm(j,3)
            wfmt4(i,2)=wfmt4(i,2)-t1*zlflm(j,3)
            wfmt4(i,3)=wfmt4(i,3)+t1*(zlflm(j,1) &
             +cmplx(aimag(zlflm(j,2)),-dble(zlflm(j,2)),8))
            i=i+1
          end do
        end do
! outer part of muffin-tin
        do irc=irco,nrc
          call lopzflm(lmaxo,wfmt1(i,jst),lmmaxo,zlflm)
          t1=socfr(irc,ias)
          do j=1,lmmaxo
            wfmt4(i,1)=wfmt4(i,1)+t1*zlflm(j,3)
            wfmt4(i,2)=wfmt4(i,2)-t1*zlflm(j,3)
            wfmt4(i,3)=wfmt4(i,3)+t1*(zlflm(j,1) &
             +cmplx(aimag(zlflm(j,2)),-dble(zlflm(j,2)),8))
            i=i+1
          end do
        end do
      end if
    else
      do k=1,nsc
        wfmt4(1:npc,k)=0.d0
      end do
    end if
! apply muffin-tin potential matrix if required
    if (tvmatmt) then
      do l=0,lmaxdm
        if (tvmmt(l,ias)) then
          nm=2*l+1
          lm=idxlm(l,-l)
          do k=1,nsc
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
               ld,wfmt1(lm,jst),lmmaxi,zone,wfmt4(lm,k),lmmaxi)
            end if
            i=npci+lm
            call zgemm('N','N',nm,nrco,nm,zone,vmatmt(lm,ispn,lm,jspn,ias),ld, &
             wfmt1(i,jst),lmmaxo,zone,wfmt4(i,k),lmmaxo)
          end do
        end if
      end do
    end if
! apply vector potential if required
    if (tafield) then
      call gradzfmt(nrc,nrci,rcmt(:,is),wfmt1(:,jst),npcmtmax,gzfmt)
      do i=1,npc
        z1=afieldc(1)*gzfmt(i,1)+afieldc(2)*gzfmt(i,2)+afieldc(3)*gzfmt(i,3)
        z1=ca*cmplx(-aimag(z1),dble(z1),8)
        wfmt4(i,1:nsd)=wfmt4(i,1:nsd)+z1
      end do
    end if
! second-variational Hamiltonian matrix
    call omp_hold(nstfv,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(i,j,k) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
    do ist=1,nstfv
      do k=1,nsc
        if (k.eq.1) then
          i=ist
          j=jst
        else if (k.eq.2) then
          i=ist+nstfv
          j=jst+nstfv
        else
          i=ist
          j=jst+nstfv
        end if
        if (i.le.j) then
          evecsv(i,j)=evecsv(i,j)+zfmtinp(nrc,nrci,rcmt(:,is),r2cmt(:,is), &
           wfmt1(:,ist),wfmt4(:,k))
        end if
      end do
    end do
!$OMP END DO
!$OMP END PARALLEL
    call omp_free(nthd)
! end loop over states
  end do
! apply tau-DFT potential if required
  if (xcgrad.eq.4) then
    do ist=1,nstfv
      call gradzfmt(nrc,nrci,rcmt(:,is),wfmt1(:,ist),npcmtmax,gwfmt(:,:,ist))
    end do
    do jst=1,nstfv
      do k=1,3
        call zbsht(nrc,nrci,gwfmt(:,k,jst),wfmt2)
        wfmt2(1:npc)=wsmt(1:npc,ias)*wfmt2(1:npc)
        call zfsht(nrc,nrci,wfmt2,gzfmt(:,k))
      end do
      do ist=1,nstfv
        zsum=0.d0
        do k=1,3
          z1=zfmtinp(nrc,nrci,rcmt(:,is),r2cmt(:,is),gwfmt(:,k,ist),gzfmt(:,k))
          zsum=zsum+z1
        end do
        zsum=0.5d0*zsum
        evecsv(ist,jst)=evecsv(ist,jst)+zsum
        if (spinpol) then
          i=ist+nstfv
          j=jst+nstfv
          evecsv(i,j)=evecsv(i,j)+zsum
        end if
      end do
    end do
  end if
! end loop over atoms
end do
deallocate(wfmt1,wfmt2,wfmt3,wfmt4)
if (tafield.or.(xcgrad.eq.4)) deallocate(gzfmt)
if (xcgrad.eq.4) deallocate(gwfmt)
!---------------------------!
!     interstitial part     !
!---------------------------!
if (spinpol.or.tafield.or.(xcgrad.eq.4)) then
  if (socz) nsc=2
  allocate(wfir1(ngtot),wfir2(ngtot),wfgp(ngkmax,nsc))
  if (xcgrad.eq.4) allocate(gwfgp(ngkmax,nstfv))
! begin loop over states
  do jst=1,nstfv
    wfir1(:)=0.d0
    do igp=1,ngp
      wfir1(igfft(igpig(igp)))=evecfv(igp,jst)
    end do
! Fourier transform wavefunction to real-space
    call zfftifc(3,ngridg,1,wfir1)
! multiply with magnetic field and transform to G-space
    if (spinpol) then
      wfir2(:)=bsir(:,ndmag)*wfir1(:)
      call zfftifc(3,ngridg,-1,wfir2)
      do igp=1,ngp
        wfgp(igp,1)=wfir2(igfft(igpig(igp)))
      end do
      wfgp(1:ngp,2)=-wfgp(1:ngp,1)
      if (ncmag) then
        wfir2(:)=cmplx(bsir(:,1),-bsir(:,2),8)*wfir1(:)
        call zfftifc(3,ngridg,-1,wfir2)
        do igp=1,ngp
          wfgp(igp,3)=wfir2(igfft(igpig(igp)))
        end do
      end if
    else
      wfgp(1:ngp,1:nsd)=0.d0
    end if
! apply vector potential if required
    if (tafield) then
      wfir1(:)=0.d0
      do igp=1,ngp
        t1=-ca*dot_product(afieldc(:),vgpc(:,igp))
        wfir1(igfft(igpig(igp)))=t1*evecfv(igp,jst)
      end do
      call zfftifc(3,ngridg,1,wfir1)
      wfir1(:)=wfir1(:)*cfunir(:)
      call zfftifc(3,ngridg,-1,wfir1)
      do igp=1,ngp
        z1=wfir1(igfft(igpig(igp)))
        wfgp(igp,1:nsd)=wfgp(igp,1:nsd)+z1
      end do
    end if
! add to Hamiltonian matrix
    call omp_hold(nstfv,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(i,j,k) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
    do ist=1,nstfv
      do k=1,nsc
        if (k.eq.1) then
          i=ist
          j=jst
        else if (k.eq.2) then
          i=ist+nstfv
          j=jst+nstfv
        else
          i=ist
          j=jst+nstfv
        end if
        if (i.le.j) then
          evecsv(i,j)=evecsv(i,j)+zdotc(ngp,evecfv(:,ist),1,wfgp(:,k),1)
        end if
      end do
    end do
!$OMP END DO
!$OMP END PARALLEL
    call omp_free(nthd)
! end loop over states
  end do
! apply tau-DFT potential if required
  if (xcgrad.eq.4) then
    do k=1,3
      do ist=1,nstfv
        do igp=1,ngp
          z1=evecfv(igp,ist)
          gwfgp(igp,ist)=vgpc(k,igp)*cmplx(-aimag(z1),dble(z1),8)
        end do
      end do
      do jst=1,nstfv
        wfir1(:)=0.d0
        do igp=1,ngp
          wfir1(igfft(igpig(igp)))=gwfgp(igp,jst)
        end do
        call zfftifc(3,ngridg,1,wfir1)
        wfir1(:)=wsir(:)*wfir1(:)
        call zfftifc(3,ngridg,-1,wfir1)
        do igp=1,ngp
          wfgp(igp,1)=wfir1(igfft(igpig(igp)))
        end do
        do ist=1,nstfv
          z1=0.5d0*zdotc(ngp,gwfgp(:,ist),1,wfgp,1)
          evecsv(ist,jst)=evecsv(ist,jst)+z1
          if (spinpol) then
            i=ist+nstfv
            j=jst+nstfv
            evecsv(i,j)=evecsv(i,j)+z1
          end if
        end do
      end do
    end do
  end if
  deallocate(wfir1,wfir2,wfgp)
  if (xcgrad.eq.4) deallocate(gwfgp)
end if
! add the diagonal first-variational part
i=0
do ispn=1,nspinor
  do ist=1,nstfv
    i=i+1
    evecsv(i,i)=evecsv(i,i)+evalfv(ist)
  end do
end do
if (spcpl.or.(.not.spinpol)) then
! spins are coupled; or spin-unpolarised: full diagonalisation
  call eveqnz(nstsv,nstsv,evecsv,evalsvp)
else
! spins not coupled: block diagonalise H
  call eveqnz(nstfv,nstsv,evecsv,evalsvp)
  i=nstfv+1
  call eveqnz(nstfv,nstsv,evecsv(i,i),evalsvp(i))
  do i=1,nstfv
    do j=1,nstfv
      evecsv(i,j+nstfv)=0.d0
      evecsv(i+nstfv,j)=0.d0
    end do
  end do
end if
call timesec(ts1)
!$OMP CRITICAL(eveqnsv_)
timesv=timesv+ts1-ts0
!$OMP END CRITICAL(eveqnsv_)
return
end subroutine

