
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rhomagk
! !INTERFACE:
subroutine rhomagk(ngp,igpig,wppt,occsvp,apwalm,evecfv,evecsv)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ngp    : number of G+p-vectors (in,integer(nspnfv))
!   igpig  : index from G+p-vectors to G-vectors (in,integer(ngkmax,nspnfv))
!   wppt   : weight of input p-point (in,real)
!   occsvp : occupation number for each state (in,real(nstsv))
!   apwalm : APW matching coefficients
!            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
!   evecfv : first-variational eigenvectors (in,complex(nmatmax,nstfv,nspnfv))
!   evecsv : second-variational eigenvectors (in,complex(nstsv,nstsv))
! !DESCRIPTION:
!   Generates the partial valence charge density and magnetisation from the
!   eigenvectors at a particular $k$-point. In the muffin-tin region, the
!   wavefunction is obtained in terms of its $(l,m)$-components from both the
!   APW and local-orbital functions. Using a backward spherical harmonic
!   transform (SHT), the wavefunction is converted to real-space and the density
!   obtained from its modulus squared. A similar proccess is used for the
!   intersitial density in which the wavefunction in real-space is obtained from
!   a Fourier transform of the APW functions. See routines {\tt wavefmt},
!   {\tt genshtmat} and {\tt eveqn}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!   Removed conversion to spherical harmonics, January 2009 (JKD)
!   Partially de-phased the muffin-tin magnetisation for spin-spirals,
!    February 2009 (FC, FB & LN)
!   Optimisations, July 2010 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ngp(nspnfv),igpig(ngkmax,nspnfv)
real(8), intent(in) :: wppt,occsvp(nstsv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv)
! local variables
integer ispn,jspn,ist,is,ia,ias
integer nr,nri,nrc,nrci,ir,irc
integer igp,ifg,i,j
real(8) wo,wo2,t1
real(8) ts0,ts1
complex(8) zq(2),z1
! automatic arrays
logical done(nstfv,nspnfv)
! allocatable arrays
complex(8), allocatable :: wfmt1(:,:,:,:),wfmt2(:,:)
complex(8), allocatable :: wfmt3(:,:,:),wfir(:,:)
call timesec(ts0)
!----------------------------------------------!
!     muffin-tin density and magnetisation     !
!----------------------------------------------!
if (tevecsv) allocate(wfmt1(lmmaxvr,nrcmtmax,nstfv,nspnfv))
allocate(wfmt2(lmmaxvr,nrcmtmax),wfmt3(lmmaxvr,nrcmtmax,nspinor))
! loop over atoms
do ias=1,natmtot
  is=idxis(ias)
  ia=idxia(ias)
  nr=nrmt(is)
  nri=nrmtinr(is)
  nrc=nrcmt(is)
  nrci=nrcmtinr(is)
! de-phasing factor for spin-spirals
  if (spinsprl.and.ssdph) then
    t1=-0.5d0*dot_product(vqcss(:),atposc(:,ia,is))
    zq(1)=cmplx(cos(t1),sin(t1),8)
    zq(2)=conjg(zq(1))
  end if
  done(:,:)=.false.
  do j=1,nstsv
    if (abs(occsvp(j)).lt.epsocc) cycle
    wo=wppt*occsvp(j)
    wo2=2.d0*wo
    if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
      i=0
      do ispn=1,nspinor
        jspn=jspnfv(ispn)
        wfmt2(:,:)=0.d0
        do ist=1,nstfv
          i=i+1
          z1=evecsv(i,j)
          if (spinsprl.and.ssdph) z1=z1*zq(ispn)
          if (abs(dble(z1))+abs(aimag(z1)).gt.epsocc) then
            if (.not.done(ist,jspn)) then
              call wavefmt(lradstp,ias,ngp(jspn),apwalm(:,:,:,:,jspn), &
               evecfv(:,ist,jspn),wfmt1(:,:,ist,jspn))
              done(ist,jspn)=.true.
            end if
! add to spinor wavefunction
            call zfmtadd(nrc,nrci,z1,wfmt1(:,:,ist,jspn),wfmt2)
          end if
        end do
! convert to spherical coordinates
        call zbsht(nrc,nrci,wfmt2,wfmt3(:,:,ispn))
      end do
    else
! spin-unpolarised wavefunction
      call wavefmt(lradstp,ias,ngp,apwalm,evecfv(:,j,1),wfmt2)
! convert to spherical coordinates
      call zbsht(nrc,nrci,wfmt2,wfmt3)
    end if
! add to density and magnetisation
!$OMP CRITICAL
    if (spinpol) then
! spin-polarised
      if (ncmag) then
! non-collinear
        irc=0
        do ir=1,nri,lradstp
          irc=irc+1
          call rmk1(lmmaxinr,wo,wo2,wfmt3(:,irc,1),wfmt3(:,irc,2), &
           rhomt(:,ir,ias),magmt(:,ir,ias,1),magmt(:,ir,ias,2), &
           magmt(:,ir,ias,3))
        end do
        do ir=nri+lradstp,nr,lradstp
          irc=irc+1
          call rmk1(lmmaxvr,wo,wo2,wfmt3(:,irc,1),wfmt3(:,irc,2), &
           rhomt(:,ir,ias),magmt(:,ir,ias,1),magmt(:,ir,ias,2), &
           magmt(:,ir,ias,3))
        end do
      else
! collinear
        irc=0
        do ir=1,nri,lradstp
          irc=irc+1
          call rmk2(lmmaxinr,wo,wfmt3(:,irc,1),wfmt3(:,irc,2),rhomt(:,ir,ias), &
           magmt(:,ir,ias,1))
        end do
        do ir=nri+lradstp,nr,lradstp
          irc=irc+1
          call rmk2(lmmaxvr,wo,wfmt3(:,irc,1),wfmt3(:,irc,2),rhomt(:,ir,ias), &
           magmt(:,ir,ias,1))
        end do
      end if
    else
! spin-unpolarised
      irc=0
      do ir=1,nri,lradstp
        irc=irc+1
        call rmk3(lmmaxinr,wo,wfmt3(:,irc,1),rhomt(:,ir,ias))
      end do
      do ir=nri+lradstp,nr,lradstp
        irc=irc+1
        call rmk3(lmmaxvr,wo,wfmt3(:,irc,1),rhomt(:,ir,ias))
      end do
    end if
!$OMP END CRITICAL
  end do
! end loop over atoms
end do
if (tevecsv) deallocate(wfmt1)
deallocate(wfmt2,wfmt3)
!------------------------------------------------!
!     interstitial density and magnetisation     !
!------------------------------------------------!
allocate(wfir(ngtot,nspinor))
do j=1,nstsv
  if (abs(occsvp(j)).lt.epsocc) cycle
  wo=wppt*occsvp(j)/omega
  wo2=2.d0*wo
  wfir(:,:)=0.d0
  if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
    i=0
    do ispn=1,nspinor
      jspn=jspnfv(ispn)
      do ist=1,nstfv
        i=i+1
        z1=evecsv(i,j)
        if (abs(dble(z1))+abs(aimag(z1)).gt.epsocc) then
          do igp=1,ngp(jspn)
            ifg=igfft(igpig(igp,jspn))
            wfir(ifg,ispn)=wfir(ifg,ispn)+z1*evecfv(igp,ist,jspn)
          end do
        end if
      end do
    end do
  else
! spin-unpolarised wavefunction
    do igp=1,ngp(1)
      ifg=igfft(igpig(igp,1))
      wfir(ifg,1)=evecfv(igp,j,1)
    end do
  end if
! Fourier transform wavefunction to real-space
  do ispn=1,nspinor
    call zfftifc(3,ngridg,1,wfir(:,ispn))
  end do
! add to density and magnetisation
!$OMP CRITICAL
  if (spinpol) then
! spin-polarised
    if (ncmag) then
! non-collinear
      call rmk1(ngtot,wo,wo2,wfir(:,1),wfir(:,2),rhoir,magir(:,1),magir(:,2), &
       magir(:,3))
    else
! collinear
      call rmk2(ngtot,wo,wfir(:,1),wfir(:,2),rhoir,magir)
    end if
  else
! spin-unpolarised
    call rmk3(ngtot,wo,wfir,rhoir)
  end if
!$OMP END CRITICAL
end do
deallocate(wfir)
call timesec(ts1)
!$OMP CRITICAL
timerho=timerho+ts1-ts0
!$OMP END CRITICAL
return

contains

subroutine rmk1(n,wo,wo2,wf1,wf2,rho,mag1,mag2,mag3)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: wo,wo2
complex(8), intent(in) :: wf1(n),wf2(n)
real(8), intent(inout) :: rho(n),mag1(n),mag2(n),mag3(n)
! local variables
integer i
real(8) t1,t2
complex(8) z1,z2
do i=1,n
  z1=wf1(i)
  z2=wf2(i)
  t1=dble(z1)**2+aimag(z1)**2
  t2=dble(z2)**2+aimag(z2)**2
  z1=conjg(z1)*z2
  rho(i)=rho(i)+wo*(t1+t2)
  mag1(i)=mag1(i)+wo2*dble(z1)
  mag2(i)=mag2(i)+wo2*aimag(z1)
  mag3(i)=mag3(i)+wo*(t1-t2)
end do
return
end subroutine

subroutine rmk2(n,wo,wf1,wf2,rho,mag)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: wo
complex(8), intent(in) :: wf1(n),wf2(n)
real(8), intent(inout) :: rho(n),mag(n)
! local variables
integer i
real(8) t1,t2
do i=1,n
  t1=dble(wf1(i))**2+aimag(wf1(i))**2
  t2=dble(wf2(i))**2+aimag(wf2(i))**2
  rho(i)=rho(i)+wo*(t1+t2)
  mag(i)=mag(i)+wo*(t1-t2)
end do
return
end subroutine

subroutine rmk3(n,wo,wf,rho)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: wo
complex(8), intent(in) :: wf(n)
real(8), intent(inout) :: rho(n)
rho(:)=rho(:)+wo*(dble(wf(:))**2+aimag(wf(:))**2)
return
end subroutine

end subroutine
!EOC
