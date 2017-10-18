
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genexpmt(vpc,expmt)
use modmain
implicit none
! arguments
real(8), intent(in) :: vpc(3)
complex(8), intent(out) :: expmt(npcmtmax,natmtot)
! local variables
integer is,ia,ias
integer nrc,nrci,irc,npc
integer lmax,l,m,lm,i
real(8) pc,tp(2),t1
complex(8) z1
! automatic arrays
real(8) jl(0:lmaxo)
complex(8) ylm(lmmaxo)
! allocatable arrays
complex(8), allocatable :: zfmt1(:),zfmt2(:)
allocate(zfmt1(npcmtmax),zfmt2(npcmtmax))
! p-vector length and (theta, phi) coordinates
call sphcrd(vpc,pc,tp)
! p-vector spherical harmonics
call genylm(lmaxo,tp,ylm)
do is=1,nspecies
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
  lmax=lmaxi
  i=0
  do irc=1,nrc
    t1=pc*rcmt(irc,is)
    call sbessel(lmax,t1,jl)
    lm=0
    do l=0,lmax
      z1=fourpi*jl(l)*zil(l)
      do m=-l,l
        lm=lm+1
        i=i+1
        zfmt1(i)=z1*conjg(ylm(lm))
      end do
    end do
    if (irc.eq.nrci) lmax=lmaxo
  end do
! convert to spherical coordinates
  call zbsht(nrc,nrci,zfmt1,zfmt2)
! mutiply by phase factors and store for all atoms
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    t1=dot_product(vpc(:),atposc(:,ia,is))
    z1=cmplx(cos(t1),sin(t1),8)
    expmt(1:npc,ias)=z1*zfmt2(1:npc)
  end do
end do
deallocate(zfmt1,zfmt2)
return
end subroutine

