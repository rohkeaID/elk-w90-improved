
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zftzf(ngp,jlgpr,ylmgp,ld,sfacgp,zfmt,zfir,zfgp)
use modmain
implicit none
! arguments
integer, intent(in) :: ngp
real(8), intent(in) :: jlgpr(njcmax,nspecies,ngp)
complex(8), intent(in) :: ylmgp(lmmaxo,ngp)
integer, intent(in) :: ld
complex(8), intent(in) :: sfacgp(ld,natmtot)
complex(8), intent(in) :: zfmt(npcmtmax,natmtot),zfir(ngtot)
complex(8), intent(out) :: zfgp(ngp)
! local variables
integer is,ia,ias,ig
integer nrc,nrci,irc
integer lmax,l,m,lm,i,j
real(8) t0,t1,t2
complex(8) zsum1,zsum2
! automatic arrays
real(8) fr1(nrcmtmax),fr2(nrcmtmax)
complex(8) ylm(lmmaxo)
! allocatable arrays
complex(8), allocatable :: zfft(:)
! external functions
real(8) fintgt
external fintgt
if ((ngp.lt.1).or.(ngp.gt.ngvec)) then
  write(*,*)
  write(*,'("Error(zftzf): ngp out of range : ",I8)') ngp
  write(*,*)
  stop
end if
! zero the coefficients
zfgp(:)=0.d0
!---------------------------------!
!     muffin-tin contribution     !
!---------------------------------!
t0=fourpi/omega
do ig=1,ngp
  ylm(:)=ylmgp(:,ig)
  do is=1,nspecies
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      lmax=lmaxi
      i=0
      j=0
      do irc=1,nrc
        i=i+1
        j=j+1
        zsum1=jlgpr(j,is,ig)*zfmt(i,ias)*ylm(1)
        lm=1
        do l=1,lmax
          lm=lm+1
          i=i+1
          j=j+1
          zsum2=zfmt(i,ias)*ylm(lm)
          do m=1-l,l
            lm=lm+1
            i=i+1
            zsum2=zsum2+zfmt(i,ias)*ylm(lm)
          end do
          zsum1=zsum1+jlgpr(j,is,ig)*zilc(l)*zsum2
        end do
        zsum1=zsum1*r2cmt(irc,is)
        fr1(irc)=dble(zsum1)
        fr2(irc)=aimag(zsum1)
        if (irc.eq.nrci) lmax=lmaxo
      end do
      t1=fintgt(-1,nrc,rcmt(:,is),fr1)
      t2=fintgt(-1,nrc,rcmt(:,is),fr2)
      zfgp(ig)=zfgp(ig)+t0*conjg(sfacgp(ig,ias))*cmplx(t1,t2,8)
    end do
  end do
end do
!-----------------------------------!
!     interstitial contribution     !
!-----------------------------------!
allocate(zfft(ngtot))
zfft(:)=zfir(:)*cfunir(:)
call zfftifc(3,ngridg,-1,zfft)
do ig=1,ngp
  zfgp(ig)=zfgp(ig)+zfft(igfft(ig))
end do
deallocate(zfft)
return
end subroutine

