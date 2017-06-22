
! Copyright (C) 2002-2009 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: gradzfmt
! !INTERFACE:
subroutine gradzfmt(nr,nri,r,zfmt,ld,gzfmt)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   nr    : number of radial mesh points (in,integer)
!   nri   : number of points on inner part of muffin-tin (in,integer)
!   r     : radial mesh (in,real(nr))
!   zfmt  : complex muffin-tin function (in,complex(lmmaxvr,nr))
!   ld    : leading dimension (in,integer)
!   gzfmt : gradient of zfmt (out,complex(lmmaxvr,ld,3))
! !DESCRIPTION:
!   Calculates the gradient of a complex muffin-tin function. In other words,
!   given the spherical harmonic expansion coefficients, $f_{lm}(r)$, of a
!   function $f({\bf r})$, the routine returns ${\bf F}_{lm}$ where
!   $$ \sum_{lm}{\bf F}_{lm}(r)Y_{lm}(\hat{\bf r})=\nabla f({\bf r}). $$
!   This is done using the formula (see, for example, V. Devanathan,
!   {\em Angular Momentum Techniques In Quantum Mechanics})
!   \begin{align*}
!    \nabla_{\mu}^s f_{lm}(r)Y_{lm}(\hat{\bf r})&=\sqrt{\frac{l+1}{2l+3}}
!    C(l,1,l+1|m,\mu,m+\mu)Y_{l+1m+\mu}(\hat{\bf r})\left(\frac{d}{dr}
!    -\frac{l}{r}\right)f_{lm}(r)\\
!    &-\sqrt{\frac{l}{2l-1}}C(l,1,l-1|m,\mu,m+\mu)Y_{l-1,m+\mu}(\hat{\bf r})
!    \left(\frac{d}{dr}+\frac{l+1}{r}\right)f_{lm}(r),
!   \end{align*}
!   where $C$ are Clebsch-Gordan coefficients and the gradient $\nabla_{\mu}^s$
!   is in terms of the spherical unit vectors $\hat{\bf e}_{\mu}$:
!   $$ \hat{\bf e}_{+1}=-\frac{\hat{\bf x}+i\hat{\bf y}}{\sqrt{2}},
!    \qquad\hat{\bf e}_0=\hat{\bf z},\qquad
!    \hat{\bf e}_{-1}=\frac{\hat{\bf x}-i\hat{\bf y}}{\sqrt{2}}. $$
!   Note that the gradient returned is in terms of the global
!   $(\hat{\bf x},\hat{\bf y},\hat{\bf z})$ coordinate system.
!
! !REVISION HISTORY:
!   Rewritten May 2009 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: r(nr)
complex(8), intent(in) :: zfmt(lmmaxvr,nr)
integer, intent(in) :: ld
complex(8), intent(out) :: gzfmt(lmmaxvr,ld,3)
! local variables
integer nro,iro,ir,i,j
integer lmmax,l,m,lm,lm1
! real constant 1/sqrt(2)
real(8), parameter :: c1=0.7071067811865475244d0
real(8) t1,t2,t3
complex(8) z1
! automatic arrays
real(8) ri(nr),f(nr),g1(nr),g2(nr)
! external functions
real(8) clebgor
external clebgor
do ir=1,nr
  ri(ir)=1.d0/r(ir)
end do
do i=1,3
  call zfmtzero(nr,nri,gzfmt(:,:,i))
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
  t1=sqrt(dble(l+1)/dble(2*l+3))
  if (l.gt.0) then
    t2=sqrt(dble(l)/dble(2*l-1))
  else
    t2=0.d0
  end if
  do m=-l,l
    lm=lm+1
! compute the radial derivatives
    f(iro:nr)=dble(zfmt(lm,iro:nr))
    call fderiv(1,nro,r(iro),f(iro),g1(iro))
    f(iro:nr)=aimag(zfmt(lm,iro:nr))
    call fderiv(1,nro,r(iro),f(iro),g2(iro))
    j=1
    do i=-1,1
      if (i.eq.0) j=3
      if (i.eq.1) j=2
      if ((l+1.le.lmaxvr).and.(abs(m+i).le.l+1)) then
! index to (l,m) is l*(l+1)+m+1, therefore index to (l+1,m+i) is
        lm1=(l+1)*(l+2)+(m+i)+1
        t3=t1*clebgor(l,1,l+1,m,i,m+i)
        do ir=iro,nr
          gzfmt(lm1,ir,j)=gzfmt(lm1,ir,j) &
           +t3*(cmplx(g1(ir),g2(ir),8)-dble(l)*ri(ir)*zfmt(lm,ir))
        end do
      end if
      if (abs(m+i).le.l-1) then
! index to (l-1,m+i)
        lm1=(l-1)*l+(m+i)+1
        t3=t2*clebgor(l,1,l-1,m,i,m+i)
        do ir=iro,nr
          gzfmt(lm1,ir,j)=gzfmt(lm1,ir,j) &
           -t3*(cmplx(g1(ir),g2(ir),8)+dble(l+1)*ri(ir)*zfmt(lm,ir))
        end do
      end if
    end do
  end do
end do
! convert from spherical components to Cartesian
lmmax=lmmaxinr
do ir=1,nr
  do lm=1,lmmax
    z1=gzfmt(lm,ir,1)
    gzfmt(lm,ir,1)=c1*(z1-gzfmt(lm,ir,2))
    z1=c1*(z1+gzfmt(lm,ir,2))
    gzfmt(lm,ir,2)=cmplx(-aimag(z1),dble(z1),8)
  end do
  if (ir.eq.nri) lmmax=lmmaxvr
end do
return
end subroutine
!EOC

