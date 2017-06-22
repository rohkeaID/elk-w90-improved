
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: zpotclmt
! !INTERFACE:
subroutine zpotclmt(nr,nri,r,zrhomt,zvclmt)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   nr      : number of radial mesh points (in,integer)
!   nri     : number of points on inner part of muffin-tin (in,integer)
!   r       : radial mesh (in,real(nr))
!   zrhomt  : muffin-tin charge density (in,complex(lmmaxvr,nr))
!   zvclmt  : muffin-tin Coulomb potential (out,complex(lmmaxvr,nr))
! !DESCRIPTION:
!   Solves the Poisson equation for the charge density contained in an isolated
!   muffin-tin using the Green's function approach. In other words, the
!   spherical harmonic expansion of the Coulomb potential, $V_{lm}$, is obtained
!   from the density expansion, $\rho_{lm}$, by
!   $$ V_{lm}(r)=\frac{4\pi}{2l+1}\left(\frac{1}{r^{l+1}}\int_0^r\rho_{lm}(r')
!      {r'}^{l+2}dr'+r^l\int_r^R\frac{\rho_{lm}(r')}{{r'}^{l-1}}dr'\right) $$
!   where $R$ is the muffin-tin radius.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: r(nr)
complex(8), intent(in) :: zrhomt(lmmaxvr,nr)
complex(8), intent(out) :: zvclmt(lmmaxvr,nr)
! local variables
integer nro,iro,ir
integer l,m,lm
real(8) t1,t2,t3,t4
! automatic arrays
real(8) ri(nr),rl1(nr),rl2(nr),rl3(nr),rl4(nr)
real(8) fr1(nr),fr2(nr),fr3(nr),fr4(nr),fr5(nr)
! initialise r^l, r^(-l-1), r^(l+2) and r^(-l+1)
do ir=1,nr
  ri(ir)=1.d0/r(ir)
  rl1(ir)=1.d0
  rl2(ir)=ri(ir)
  t1=fourpi*r(ir)
  rl3(ir)=t1*r(ir)
  rl4(ir)=t1
end do
lm=0
do l=0,lmaxvr
  if (l.le.lmaxinr) then
    nro=nr
    iro=1
  else
    nro=nr-nri
    iro=nri+1
  end if
  do m=-l,l
    lm=lm+1
    do ir=iro,nr
      t1=dble(zrhomt(lm,ir))
      t2=aimag(zrhomt(lm,ir))
      fr1(ir)=t1*rl3(ir)
      fr2(ir)=t2*rl3(ir)
      fr3(ir)=t1*rl4(ir)
      fr4(ir)=t2*rl4(ir)
    end do
    call fderiv(-1,nro,r(iro),fr1(iro),fr5(iro))
    call fderiv(-1,nro,r(iro),fr2(iro),fr1(iro))
    call fderiv(-1,nro,r(iro),fr3(iro),fr2(iro))
    call fderiv(-1,nro,r(iro),fr4(iro),fr3(iro))
    t1=fr2(nr)
    t2=fr3(nr)
    do ir=iro,nr
      t3=rl2(ir)*fr5(ir)+rl1(ir)*(t1-fr2(ir))
      t4=rl2(ir)*fr1(ir)+rl1(ir)*(t2-fr3(ir))
      zvclmt(lm,ir)=cmplx(t3,t4,8)
    end do
  end do
  if (l.lt.lmaxvr) then
    t1=fourpi/dble(2*(l+1)+1)
    do ir=iro,nr
      rl1(ir)=rl1(ir)*r(ir)
      rl2(ir)=rl2(ir)*ri(ir)
      t2=t1*r(ir)**2
      rl3(ir)=rl1(ir)*t2
      rl4(ir)=rl2(ir)*t2
    end do
  end if
end do
return
end subroutine
!EOC

