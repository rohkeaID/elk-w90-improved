
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dwavefmt(lrstp,ias,ngp,ngpq,apwalmq,dapwalm,evecfv,devecfv,dwfmt)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: lrstp,ias,ngp,ngpq
complex(8), intent(in) :: apwalmq(ngkmax,apwordmax,lmmaxapw,natmtot)
!**** remove natmtot

complex(8), intent(in) :: dapwalm(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: evecfv(nmatmax),devecfv(nmatmax)
real(8), intent(out) :: dwfmt(2,*)
! local variables
integer is,ldi,ldo,io,ilo
integer nrc,nrci,nrco,iro
integer l,m,lm,npc,npci,i
complex(8) z1
! external functions
complex(8) zdotu
external zdotu
is=idxis(ias)
ldi=2*lmmaxi
ldo=2*lmmaxo
iro=nrmti(is)+lrstp
if (lrstp.eq.1) then
  nrc=nrmt(is)
  nrci=nrmti(is)
  npc=npmt(is)
  npci=npmti(is)
else if (lrstp.eq.lradstp) then
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
  npci=npcmti(is)
else
  write(*,*)
  write(*,'("Error(dwavefmt): invalid lrstp : ",I8)') lrstp
  write(*,*)
  stop
end if
nrco=nrc-nrci
! zero the wavefunction derivative
dwfmt(:,1:npc)=0.d0
!-----------------------!
!     APW functions     !
!-----------------------!
lm=0
do l=0,lmaxo
  do m=-l,l
    lm=lm+1
    i=npci+lm
    do io=1,apword(l,is)
      z1=zdotu(ngpq,devecfv,1,apwalmq(:,io,lm,ias),1)
      if (ias.eq.iasph) then
        z1=z1+zdotu(ngp,evecfv,1,dapwalm(:,io,lm),1)
      end if
      if (abs(dble(z1)).gt.1.d-14) then
        if (l.le.lmaxi) then
          call daxpy(nrci,dble(z1),apwfr(1,1,io,l,ias),lrstp,dwfmt(1,lm),ldi)
        end if
        call daxpy(nrco,dble(z1),apwfr(iro,1,io,l,ias),lrstp,dwfmt(1,i),ldo)
      end if
      if (abs(aimag(z1)).gt.1.d-14) then
        if (l.le.lmaxi) then
          call daxpy(nrci,aimag(z1),apwfr(1,1,io,l,ias),lrstp,dwfmt(2,lm),ldi)
        end if
        call daxpy(nrco,aimag(z1),apwfr(iro,1,io,l,ias),lrstp,dwfmt(2,i),ldo)
      end if
    end do
  end do
end do
!---------------------------------!
!     local-orbital functions     !
!---------------------------------!
do ilo=1,nlorb(is)
  l=lorbl(ilo,is)
  do m=-l,l
    lm=idxlm(l,m)
    i=npci+lm
    z1=devecfv(ngpq+idxlo(lm,ilo,ias))
    if (abs(dble(z1)).gt.1.d-14) then
      if (l.le.lmaxi) then
        call daxpy(nrci,dble(z1),lofr(1,1,ilo,ias),lrstp,dwfmt(1,lm),ldi)
      end if
      call daxpy(nrco,dble(z1),lofr(iro,1,ilo,ias),lrstp,dwfmt(1,i),ldo)
    end if
    if (abs(aimag(z1)).gt.1.d-14) then
      if (l.le.lmaxi) then
        call daxpy(nrci,aimag(z1),lofr(1,1,ilo,ias),lrstp,dwfmt(2,lm),ldi)
      end if
      call daxpy(nrco,aimag(z1),lofr(iro,1,ilo,ias),lrstp,dwfmt(2,i),ldo)
    end if
  end do
end do
return
end subroutine

