
! Copyright (C) 2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zftgq0(jlgq0r,ylmgq0,sfacgq0,zrhomt,zrhoir,zgq0)
use modmain
implicit none
! arguments
real(8), intent(in) :: jlgq0r(0:lmaxo,nrcmtmax,nspecies)
complex(8), intent(in) :: ylmgq0(lmmaxo),sfacgq0(natmtot)
complex(8), intent(in) :: zrhomt(npcmtmax,natmtot),zrhoir(ngtot)
complex(8), intent(out) :: zgq0
! local variables
integer is,ias
integer nrc,nrci,ir,irc
integer lmax,l,m,lm,i
real(8) t0,t1,t2
complex(8) zsum1,zsum2
! automatic arrays
real(8) fr1(nrcmtmax),fr2(nrcmtmax)
! external functions
real(8) fintgt
external fintgt
!-----------------------------------!
!     interstitial contribution     !
!-----------------------------------!
! (note that the phase exp(i(G+q).r) is implicit)
zgq0=cfunir(1)*zrhoir(1)
do ir=2,ngtot
  zgq0=zgq0+cfunir(ir)*zrhoir(ir)
end do
zgq0=zgq0/dble(ngtot)
!---------------------------------!
!     muffin-tin contribution     !
!---------------------------------!
! (note that the phase exp(i(G+q).r) is explicit)
t0=fourpi/omega
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  lmax=lmaxi
  i=0
  do irc=1,nrc
    i=i+1
    zsum1=jlgq0r(0,irc,is)*zrhomt(i,ias)*ylmgq0(1)
    lm=1
    do l=1,lmax
      lm=lm+1
      i=i+1
      zsum2=zrhomt(i,ias)*ylmgq0(lm)
      do m=1-l,l
        lm=lm+1
        i=i+1
        zsum2=zsum2+zrhomt(i,ias)*ylmgq0(lm)
      end do
      zsum1=zsum1+jlgq0r(l,irc,is)*zilc(l)*zsum2
    end do
    zsum1=zsum1*r2cmt(irc,is)
    fr1(irc)=dble(zsum1)
    fr2(irc)=aimag(zsum1)
    if (irc.eq.nrci) lmax=lmaxo
  end do
  t1=fintgt(-1,nrc,rcmt(:,is),fr1)
  t2=fintgt(-1,nrc,rcmt(:,is),fr2)
  zgq0=zgq0+t0*conjg(sfacgq0(ias))*cmplx(t1,t2,8)
end do
return
end subroutine

