
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dwavefmt(lrstp,lmax,ias,ngp,ngpq,apwalmq,dapwalm,evecfv,devecfv,ld, &
 dwfmt)
!******* change me!!!
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: lrstp,lmax,ias,ngp,ngpq
complex(8), intent(in) :: apwalmq(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(in) :: dapwalm(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: evecfv(nmatmax),devecfv(nmatmax)
integer, intent(in) :: ld
real(8), intent(out) :: dwfmt(2,ld,*)
! local variables
integer is,nrc,ir,ld2
integer l,m,lm,io,ilo
complex(8) z1
! external functions
complex(8) zdotu
external zdotu
if (lmax.gt.lmaxapw) then
  write(*,*)
  write(*,'("Error(dwavefmt): lmax > lmaxapw : ",I8)') lmax
  write(*,*)
  stop
end if
ld2=ld*2
is=idxis(ias)
nrc=0
do ir=1,nrmt(is),lrstp
  nrc=nrc+1
  dwfmt(:,:,nrc)=0.d0
end do
! APW functions
lm=0
do l=0,lmax
  do m=-l,l
    lm=lm+1
    do io=1,apword(l,is)
      z1=zdotu(ngpq,devecfv,1,apwalmq(:,io,lm,ias),1)
      if (ias.eq.iasph) then
        z1=z1+zdotu(ngp,evecfv,1,dapwalm(:,io,lm),1)
      end if
      if (abs(dble(z1)).gt.1.d-14) then
        call daxpy(nrc,dble(z1),apwfr(:,1,io,l,ias),lrstp,dwfmt(1,lm,1),ld2)
      end if
      if (abs(aimag(z1)).gt.1.d-14) then
        call daxpy(nrc,aimag(z1),apwfr(:,1,io,l,ias),lrstp,dwfmt(2,lm,1),ld2)
      end if
    end do
  end do
end do
! local-orbital functions
do ilo=1,nlorb(is)
  l=lorbl(ilo,is)
  if (l.le.lmax) then
    do m=-l,l
      lm=idxlm(l,m)
      z1=devecfv(ngpq+idxlo(lm,ilo,ias))
      if (abs(dble(z1)).gt.1.d-14) then
        call daxpy(nrc,dble(z1),lofr(:,1,ilo,ias),lrstp,dwfmt(1,lm,1),ld2)
      end if
      if (abs(aimag(z1)).gt.1.d-14) then
        call daxpy(nrc,aimag(z1),lofr(:,1,ilo,ias),lrstp,dwfmt(2,lm,1),ld2)
      end if
    end do
  end if
end do
return
end subroutine

