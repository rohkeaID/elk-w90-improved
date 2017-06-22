
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genexpmt(vpc,expmt)
use modmain
implicit none
! arguments
real(8), intent(in) :: vpc(3)
complex(8), intent(out) :: expmt(lmmaxvr,nrcmtmax,natmtot)
! local variables
integer is,ia,ias
integer nrc,nrci,irc
integer lmax,l,m,lm
real(8) pc,tp(2),t1
complex(8) z1
! automatic arrays
real(8) jl(0:lmaxvr)
complex(8) ylm(lmmaxvr)
! allocatable arrays
complex(8), allocatable :: zfmt1(:,:),zfmt2(:,:)
allocate(zfmt1(lmmaxvr,nrcmtmax),zfmt2(lmmaxvr,nrcmtmax))
! p-vector length and (theta, phi) coordinates
call sphcrd(vpc,pc,tp)
! p-vector spherical harmonics
call genylm(lmaxvr,tp,ylm)
do is=1,nspecies
  nrc=nrcmt(is)
  nrci=nrcmtinr(is)
  do irc=1,nrc
    if (irc.le.nrci) then
      lmax=lmaxinr
    else
      lmax=lmaxvr
    end if
    t1=pc*rcmt(irc,is)
    call sbessel(lmax,t1,jl)
    lm=0
    do l=0,lmax
      z1=fourpi*jl(l)*zil(l)
      do m=-l,l
        lm=lm+1
        zfmt1(lm,irc)=z1*conjg(ylm(lm))
      end do
    end do
  end do
! convert to spherical coordinates
  call zbsht(nrc,nrci,zfmt1,zfmt2)
! mutiply by phase factors and store for all atoms
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    t1=dot_product(vpc(:),atposc(:,ia,is))
    z1=cmplx(cos(t1),sin(t1),8)
    do irc=1,nrci
      expmt(1:lmmaxinr,irc,ias)=z1*zfmt2(1:lmmaxinr,irc)
    end do
    do irc=nrci+1,nrc
      expmt(:,irc,ias)=z1*zfmt2(:,irc)
    end do
  end do
end do
deallocate(zfmt1,zfmt2)
return
end subroutine

