
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
complex(8), intent(in) :: zfmt(npcmtmax,natmtot),zfir(ngtc)
complex(8), intent(out) :: zfgp(ngp)
! local variables
integer is,ia,ias,ig
integer nrc,nrci,irco,irc
integer l,m,lm,i,j
real(8) t0,t1,t2
complex(8) zsum1,zsum2,z1
! automatic arrays
real(8) fr1(nrcmtmax),fr2(nrcmtmax)
complex(8) ylm(lmmaxo)
! allocatable arrays
complex(8), allocatable :: zfft(:)
! external functions
real(8) fintgt
external fintgt
if ((ngp.lt.1).or.(ngp.gt.ngvc)) then
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
  lm=1
  do l=1,lmaxo
    z1=zilc(l)
    do m=-l,l
      lm=lm+1
      ylm(lm)=z1*ylmgp(lm,ig)
    end do
  end do
  do is=1,nspecies
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    irco=nrci+1
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      i=1
      j=1
! inner part of muffin-tin
      if (lmaxi.eq.1) then
        do irc=1,nrci
          z1=jlgpr(j,is,ig)*zfmt(i,ias)*y00+jlgpr(j+1,is,ig) &
           *(zfmt(i+1,ias)*ylm(2)+zfmt(i+2,ias)*ylm(3)+zfmt(i+3,ias)*ylm(4))
          z1=z1*r2cmt(irc,is)
          fr1(irc)=dble(z1)
          fr2(irc)=aimag(z1)
          i=i+4
          j=j+2
        end do
      else
        do irc=1,nrci
          zsum1=jlgpr(j,is,ig)*zfmt(i,ias)*y00
          i=i+1
          j=j+1
          lm=1
          do l=1,lmaxi
            lm=lm+1
            zsum2=zfmt(i,ias)*ylm(lm)
            i=i+1
            do m=1-l,l
              lm=lm+1
              zsum2=zsum2+zfmt(i,ias)*ylm(lm)
              i=i+1
            end do
            zsum1=zsum1+jlgpr(j,is,ig)*zsum2
            j=j+1
          end do
          zsum1=zsum1*r2cmt(irc,is)
          fr1(irc)=dble(zsum1)
          fr2(irc)=aimag(zsum1)
        end do
      end if
! outer part of muffin-tin
      do irc=irco,nrc
        zsum1=jlgpr(j,is,ig)*zfmt(i,ias)*y00
        i=i+1
        j=j+1
        lm=1
        do l=1,lmaxo
          lm=lm+1
          zsum2=zfmt(i,ias)*ylm(lm)
          i=i+1
          do m=1-l,l
            lm=lm+1
            zsum2=zsum2+zfmt(i,ias)*ylm(lm)
            i=i+1
          end do
          zsum1=zsum1+jlgpr(j,is,ig)*zsum2
          j=j+1
        end do
        zsum1=zsum1*r2cmt(irc,is)
        fr1(irc)=dble(zsum1)
        fr2(irc)=aimag(zsum1)
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
allocate(zfft(ngtc))
! multiply by coarse characteristic function
zfft(:)=zfir(:)*cfrc(:)
! Fourier transform to coarse G-grid
call zfftifc(3,ngdc,-1,zfft)
do ig=1,ngp
  zfgp(ig)=zfgp(ig)+zfft(igfc(ig))
end do
deallocate(zfft)
return
end subroutine

