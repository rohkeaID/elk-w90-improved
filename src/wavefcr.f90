
! Copyright (C) 2002-2011 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine wavefcr(tsh,lrstp,is,ia,ist,m,ld,wfcr)
use modmain
implicit none
! arguments
logical, intent(in) :: tsh
integer, intent(in) :: lrstp
integer, intent(in) :: is,ia
integer, intent(in) :: ist
! pass in m-1/2
integer, intent(in) :: m
integer, intent(in) :: ld
complex(8), intent(out) :: wfcr(lmmaxvr,ld,2)
! local variables
integer ias,nri,ir,irc
integer k,l,lm,lm1
real(8) c1,c2,t1,t2,t3
l=lsp(ist,is)
k=ksp(ist,is)
if (((k.ne.l+1).and.(k.ne.l)).or.(m.lt.-k).or.(m.gt.k-1)) then
  write(*,*)
  write(*,'("Error(wavefcr): mismatched l, k or m : ",3I4)') l,k,m
  write(*,'(" for species ",I4)') is
  write(*,'(" atom ",I4)') ia
  write(*,'(" and state ",I6)') ist
  write(*,*)
  stop
end if
if (l.gt.lmaxinr) then
  write(*,*)
  write(*,'("Error(wavefcr): l > lmaxinr : ",2I8)') l,lmaxinr
  write(*,*)
  stop
end if
ias=idxas(ia,is)
! calculate the Clebsch-Gordon coefficients
t1=sqrt(dble(l+m+1)/dble(2*l+1))
t2=sqrt(dble(l-m)/dble(2*l+1))
if (k.eq.l+1) then
  c1=t1
  c2=t2
else
  c1=t2
  c2=-t1
end if
if (abs(m).le.l) then
  lm=idxlm(l,m)
else
  lm=0
end if
if (abs(m+1).le.l) then
  lm1=idxlm(l,m+1)
else
  lm1=0
end if
if ((tsh).or.(lm.eq.0)) wfcr(:,:,1)=0.d0
if ((tsh).or.(lm1.eq.0)) wfcr(:,:,2)=0.d0
nri=nrmtinr(is)
irc=0
do ir=1,nrmt(is),lrstp
  irc=irc+1
! major component of radial wavefunction
  t1=rwfcr(ir,1,ist,ias)/rsp(ir,is)
  if (tsh) then
    if (lm.gt.0) wfcr(lm,irc,1)=t1*c1
    if (lm1.gt.0) wfcr(lm1,irc,2)=t1*c2
  else
    t2=t1*c1
    t3=t1*c2
    if (ir.le.nri) then
! inner part of muffin-tin
      if (lm.gt.0) wfcr(1:lmmaxinr,irc,1)=t2*zbshtinr(1:lmmaxinr,lm)
      if (lm1.gt.0) wfcr(1:lmmaxinr,irc,2)=t3*zbshtinr(1:lmmaxinr,lm1)
    else
! outer part of muffin-tin
      if (lm.gt.0) wfcr(:,irc,1)=t2*zbshtvr(:,lm)
      if (lm1.gt.0) wfcr(:,irc,2)=t3*zbshtvr(:,lm1)
    end if
  end if
end do
return
end subroutine

