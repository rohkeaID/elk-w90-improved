
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
!   nr     : number of radial mesh points (in,integer)
!   nri    : number of points on inner part of muffin-tin (in,integer)
!   r      : radial mesh (in,real(nr))
!   zrhomt : muffin-tin charge density (in,complex(*))
!   zvclmt : muffin-tin Coulomb potential (out,complex(*))
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
complex(8), intent(in) :: zrhomt(*)
complex(8), intent(out) :: zvclmt(*)
! local variables
integer nro,iro,ir
integer l,m,lm,npi,i
real(8) t1,t2,t3,t4
! automatic arrays
real(8) ri(nr),r1(nr),r2(nr),r3(nr),r4(nr)
real(8) f1(nr),f2(nr),f3(nr),f4(nr),f5(nr)
! initialise r^l, r^(-l-1), r^(l+2) and r^(-l+1)
do ir=1,nr
  ri(ir)=1.d0/r(ir)
  r1(ir)=1.d0
  r2(ir)=ri(ir)
  t1=fourpi*r(ir)
  r3(ir)=t1*r(ir)
  r4(ir)=t1
end do
nro=nr-nri
iro=nri+1
npi=lmmaxi*nri
lm=0
do l=0,lmaxi
  if (l.gt.0) then
    t1=fourpi/dble(2*l+1)
    do ir=1,nr
      r1(ir)=r1(ir)*r(ir)
      r2(ir)=r2(ir)*ri(ir)
      t2=t1*r(ir)**2
      r3(ir)=r1(ir)*t2
      r4(ir)=r2(ir)*t2
    end do
  end if
  do m=-l,l
    lm=lm+1
    i=lm
    do ir=1,nri
      t1=dble(zrhomt(i)); t2=aimag(zrhomt(i))
      f1(ir)=t1*r3(ir); f2(ir)=t2*r3(ir)
      f3(ir)=t1*r4(ir); f4(ir)=t2*r4(ir)
      i=i+lmmaxi
    end do
    do ir=iro,nr
      t1=dble(zrhomt(i)); t2=aimag(zrhomt(i))
      f1(ir)=t1*r3(ir); f2(ir)=t2*r3(ir)
      f3(ir)=t1*r4(ir); f4(ir)=t2*r4(ir)
      i=i+lmmaxo
    end do
    call fderiv(-1,nr,r,f1,f5)
    call fderiv(-1,nr,r,f2,f1)
    call fderiv(-1,nr,r,f3,f2)
    call fderiv(-1,nr,r,f4,f3)
    t1=f2(nr); t2=f3(nr)
    i=lm
    do ir=1,nri
      t3=r2(ir)*f5(ir)+r1(ir)*(t1-f2(ir))
      t4=r2(ir)*f1(ir)+r1(ir)*(t2-f3(ir))
      zvclmt(i)=cmplx(t3,t4,8)
      i=i+lmmaxi
    end do
    do ir=iro,nr
      t3=r2(ir)*f5(ir)+r1(ir)*(t1-f2(ir))
      t4=r2(ir)*f1(ir)+r1(ir)*(t2-f3(ir))
      zvclmt(i)=cmplx(t3,t4,8)
      i=i+lmmaxo
    end do
  end do
end do
do l=lmaxi+1,lmaxo
  t1=fourpi/dble(2*l+1)
  do ir=iro,nr
    r1(ir)=r1(ir)*r(ir)
    r2(ir)=r2(ir)*ri(ir)
    t2=t1*r(ir)**2
    r3(ir)=r1(ir)*t2
    r4(ir)=r2(ir)*t2
  end do
  do m=-l,l
    lm=lm+1
    i=npi+lm
    do ir=iro,nr
      t1=dble(zrhomt(i)); t2=aimag(zrhomt(i))
      f1(ir)=t1*r3(ir); f2(ir)=t2*r3(ir)
      f3(ir)=t1*r4(ir); f4(ir)=t2*r4(ir)
      i=i+lmmaxo
    end do
    call fderiv(-1,nro,r(iro),f1(iro),f5(iro))
    call fderiv(-1,nro,r(iro),f2(iro),f1(iro))
    call fderiv(-1,nro,r(iro),f3(iro),f2(iro))
    call fderiv(-1,nro,r(iro),f4(iro),f3(iro))
    t1=f2(nr); t2=f3(nr)
    i=npi+lm
    do ir=iro,nr
      t3=r2(ir)*f5(ir)+r1(ir)*(t1-f2(ir))
      t4=r2(ir)*f1(ir)+r1(ir)*(t2-f3(ir))
      zvclmt(i)=cmplx(t3,t4,8)
      i=i+lmmaxo
    end do
  end do
end do
return
end subroutine
!EOC

