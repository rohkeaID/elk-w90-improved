
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zftzf(ngp,gpc,ylmgp,ld,sfacgp,zfmt,zfir,zfgp)
use modmain
implicit none
! arguments
integer, intent(in) :: ngp
real(8), intent(in) :: gpc(ngp)
complex(8), intent(in) :: ylmgp(lmmaxvr,ngp)
integer, intent(in) :: ld
complex(8), intent(in) :: sfacgp(ld,natmtot)
complex(8), intent(in) :: zfmt(lmmaxvr,nrcmtmax,natmtot),zfir(ngtot)
complex(8), intent(out) :: zfgp(ngp)
! local variables
integer is,ia,ias,ig
integer nrc,nrci,irc
integer lmax,l,m,lm
real(8) t0,t1,t2
complex(8) zsum1,zsum2,z1
! automatic arrays
real(8) jl(0:lmaxvr,nrcmtmax)
real(8) fr1(nrcmtmax),fr2(nrcmtmax),gr(nrcmtmax)
complex(8) ylm(lmmaxvr)
! allocatable arrays
complex(8), allocatable :: zfft(:)
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
    nrci=nrcmtinr(is)
! generate spherical Bessel functions
    do irc=1,nrc
      if (irc.le.nrci) then
        lmax=lmaxinr
      else
        lmax=lmaxvr
      end if
      t1=gpc(ig)*rcmt(irc,is)
      call sbessel(lmax,t1,jl(:,irc))
    end do
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      z1=t0*conjg(sfacgp(ig,ias))
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
      zfgp(ig)=zfgp(ig)+z1*cmplx(t1,t2,8)
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

