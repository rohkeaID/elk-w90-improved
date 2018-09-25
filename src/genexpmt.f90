
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genexpmt(ngp,jlgpr,ylmgp,ld,sfacgp,expmt)
use modmain
implicit none
! arguments
integer, intent(in) :: ngp
real(8), intent(in) :: jlgpr(njcmax,nspecies,ngp)
complex(8), intent(in) :: ylmgp(lmmaxo,ngp)
integer, intent(in) :: ld
complex(8), intent(in) :: sfacgp(ld,natmtot)
complex(8), intent(out) :: expmt(npcmtmax,natmtot,ngp)
! local variables
integer ig,is,ia,ias
integer nrc,nrci,irc,npc
integer lmax,l,m,lm,i,j
complex(8) z1
! automatic arrays
complex(8) ylm(lmmaxo)
complex(8) zfmt1(npcmtmax),zfmt2(npcmtmax)
do ig=1,ngp
  ylm(:)=conjg(ylmgp(:,ig))
  do is=1,nspecies
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    npc=npcmt(is)
    lmax=lmaxi
    i=0
    j=0
    do irc=1,nrc
      lm=0
      do l=0,lmax
        j=j+1
        z1=jlgpr(j,is,ig)*zil(l)
        do m=-l,l
          lm=lm+1
          i=i+1
          zfmt1(i)=z1*ylm(lm)
        end do
      end do
      if (irc.eq.nrci) lmax=lmaxo
    end do
! convert to spherical coordinates
    call zbsht(nrc,nrci,zfmt1,zfmt2)
! mutiply by phase factors and store for all atoms
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      z1=fourpi*sfacgp(ig,ias)
      expmt(1:npc,ias,ig)=z1*zfmt2(1:npc)
    end do
  end do
end do
return
end subroutine

