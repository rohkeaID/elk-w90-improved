
! Copyright (C) 2010 Alexey I. Baranov.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: zftrf
! !INTERFACE:
subroutine zftrf(npv,ivp,vpc,rfmt,rfir,zfp)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   npv  : number of P-vectors (in,integer)
!   ivp  : integer coordinates of the P-vectors (in,integer(3,npv))
!   vpc  : P-vectors in Cartesian coordinates (in,real(3,npv))
!   rfmt : real muffin-tin function (in,real(lmmaxvr,nrmtmax,natmtot))
!   rfir : real interstitial function (in,real(ngtot))
!   zfp  : Fourier expansion coefficients of the real-space function
!          (out,complex(npv))
! !DESCRIPTION:
!   Given a real function periodic in the unit cell, $f({\bf r})$, this routine
!   calculates its complex Fourier expansion coefficients:
!   $$ f({\bf P})=\frac{1}{\Omega}\int d^3r\,f({\bf r})\tilde{\Theta}({\bf r})
!    e^{-i{\bf P}\cdot{\bf r}}
!    +\frac{4\pi}{\Omega}\sum_{\alpha}e^{-i{\bf P}\cdot{\bf R}_{\alpha}}
!    \sum_{lm}(-i)^l Y_{lm}(\hat{\bf P})
!    \int_{0}^{R_{\alpha}}dr\,r^2 j_{l}(|{\bf P}|r)f_{lm}^{\alpha}(r), $$
!   where $\tilde{\Theta}$ is the smooth characteristic function of the
!   interstitial region, $\Omega$ is the unit cell volume and $R_{\alpha}$ is
!   the muffin-tin radius of atom $\alpha$.
!
! !REVISION HISTORY:
!   Created July 2010 (Alexey I. Baranov)
!   Modified, November 2010 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: npv,ivp(3,npv)
real(8), intent(in) :: vpc(3,npv)
real(8), intent(in) :: rfmt(lmmaxvr,nrmtmax,natmtot),rfir(ngtot)
complex(8), intent(out) :: zfp(npv)
! local variables
integer is,ia,ias
integer nrc,nrci,iro,irc
integer lmax,l,m,lm
integer ip,ig,ifg
real(8) p,tp(2),t0,t1,t2
complex(8) zsum1,zsum2
complex(8) z1,z2,z3
! automatic arrays
real(8) jl(0:lmaxvr,nrcmtmax)
real(8) fr1(nrcmtmax),fr2(nrcmtmax),gr(nrcmtmax)
complex(8) ylm(lmmaxvr)
! allocatable arrays
complex(8), allocatable :: zfft(:),zfmt(:,:,:)
allocate(zfft(ngtot))
allocate(zfmt(lmmaxvr,nrcmtmax,natmtot))
! zero the coefficients
zfp(:)=0.d0
!---------------------------!
!     interstitial part     !
!---------------------------!
! Fourier transform to G-space
zfft(:)=rfir(:)
call zfftifc(3,ngridg,-1,zfft)
! find coefficients for all required input vectors
do ip=1,npv
  if ((ivp(1,ip).ge.intgv(1,1)).and.(ivp(1,ip).le.intgv(2,1)).and. &
      (ivp(2,ip).ge.intgv(1,2)).and.(ivp(2,ip).le.intgv(2,2)).and. &
      (ivp(3,ip).ge.intgv(1,3)).and.(ivp(3,ip).le.intgv(2,3))) then
    ig=ivgig(ivp(1,ip),ivp(2,ip),ivp(3,ip))
    zfp(ip)=zfft(igfft(ig))
  end if
end do
!-------------------------!
!     muffin-tin part     !
!-------------------------!
! convert function from real to complex spherical harmonic expansion
do ias=1,natmtot
  is=idxis(ias)
  call rtozfmt(nrcmt(is),nrcmtinr(is),lradstp,rfmt(:,:,ias),1,zfmt(:,:,ias))
end do
! remove continuation of interstitial function into muffin-tin
do ig=1,ngtot
  ifg=igfft(ig)
! spherical coordinates of G
  call sphcrd(vgc(:,ig),t1,tp)
! conjugate spherical harmonics Y_lm*(G)
  call genylm(lmaxvr,tp,ylm)
  ylm(:)=conjg(ylm(:))
  do is=1,nspecies
    nrc=nrcmt(is)
    nrci=nrcmtinr(is)
! generate spherical Bessel functions
    do irc=1,nrc
      t1=gc(ig)*rcmt(irc,is)
      call sbessel(lmaxvr,t1,jl(:,irc))
    end do
    do ia=1,natoms(is)
      ias=idxas(ia,is)
! structure factor
      t1=dot_product(vgc(:,ig),atposc(:,ia,is))
      z1=fourpi*zfft(ifg)*cmplx(cos(t1),sin(t1),8)
      lm=0
      do l=0,lmaxvr
        if (l.le.lmaxinr) then
          iro=1
        else
          iro=nrci+1
        end if
        z2=z1*zil(l)
        do m=-l,l
          lm=lm+1
          z3=z2*ylm(lm)
          do irc=iro,nrc
            zfmt(lm,irc,ias)=zfmt(lm,irc,ias)-z3*jl(l,irc)
          end do
        end do
      end do
    end do
  end do
end do
t0=fourpi/omega
! loop over input P-vectors
do ip=1,npv
! spherical coordinates of P
  call sphcrd(vpc(:,ip),p,tp)
! spherical harmonics Y_lm(P)
  call genylm(lmaxvr,tp,ylm)
  do is=1,nspecies
    nrc=nrcmt(is)
    nrci=nrcmtinr(is)
! generate spherical Bessel functions
    do irc=1,nrc
      t1=p*rcmt(irc,is)
      call sbessel(lmaxvr,t1,jl(:,irc))
    end do
    do ia=1,natoms(is)
      ias=idxas(ia,is)
! conjugate structure factor
      t1=-dot_product(vpc(:,ip),atposc(:,ia,is))
      z1=t0*cmplx(cos(t1),sin(t1),8)
      do irc=1,nrc
        if (irc.le.nrci) then
          lmax=lmaxinr
        else
          lmax=lmaxvr
        end if
        zsum1=jl(0,irc)*zfmt(1,irc,ias)*ylm(1)
        lm=1
        do l=1,lmax
          lm=lm+1
          zsum2=zfmt(lm,irc,ias)*ylm(lm)
          do m=1-l,l
            lm=lm+1
            zsum2=zsum2+zfmt(lm,irc,ias)*ylm(lm)
          end do
          zsum1=zsum1+jl(l,irc)*zilc(l)*zsum2
        end do
        zsum1=zsum1*r2cmt(irc,is)
        fr1(irc)=dble(zsum1)
        fr2(irc)=aimag(zsum1)
      end do
      call fderiv(-1,nrc,rcmt(:,is),fr1,gr)
      t1=gr(nrc)
      call fderiv(-1,nrc,rcmt(:,is),fr2,gr)
      t2=gr(nrc)
      zfp(ip)=zfp(ip)+z1*cmplx(t1,t2,8)
    end do
  end do
end do
deallocate(zfft,zfmt)
return
end subroutine
! EOC

