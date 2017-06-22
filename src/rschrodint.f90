
! Copyright (C) 2002-2015 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rschrodint
! !INTERFACE:
subroutine rschrodint(sol,l,e,nr,r,vr,nn,p0,p1,q0,q1)
! !INPUT/OUTPUT PARAMETERS:
!   sol : speed of light in atomic units (in,real)
!   l   : angular momentum quantum number (in,integer)
!   e   : energy (in,real)
!   nr  : number of radial mesh points (in,integer)
!   r   : radial mesh (in,real(nr))
!   vr  : potential on radial mesh (in,real(nr))
!   nn  : number of nodes (out,integer)
!   p0  : m th energy derivative of P (out,real(nr))
!   p1  : radial derivative of p0 (out,real(nr))
!   q0  : m th energy derivative of Q (out,real(nr))
!   q1  : radial derivative of q0 (out,real(nr))
! !DESCRIPTION:
!   Integrates the scalar relativistic radial Schr\"{o}dinger equation from
!   $r=0$ outwards. This involves using the predictor-corrector method to solve
!   the coupled first-order equations (in atomic units)
!   \begin{align*}
!    \frac{d}{dr}P_l&=2MQ_l+\frac{1}{r}P_l\\
!    \frac{d}{dr}Q_l&=-\frac{1}{r}Q_l+\left[\frac{l(l+1)}{2Mr^2}
!    +(V-E)\right]P_l,
!   \end{align*}
!   where $V$ is the external potential, $E$ is the eigen energy and
!   $M=1+(E-V)/2c^2$. Following the convention of Koelling and Harmon,
!   {\it J. Phys. C: Solid State Phys.} {\bf 10}, 3107 (1977), the functions
!   $P_l$ and $Q_l$ are defined by
!   \begin{align*}
!    P_l&=rg_l\\
!    Q_l&=\frac{r}{2M}\frac{dg_l}{dr},
!   \end{align*}
!   where $g_l$ is the major component of the Dirac equation (see the routine
!   {\tt rdiracint}).
!
! !REVISION HISTORY:
!   Created October 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: sol
integer, intent(in) :: l
real(8), intent(in) :: e
integer, intent(in) :: nr
real(8), intent(in) :: r(nr),vr(nr)
integer, intent(out) :: nn
real(8), intent(out) :: p0(nr),p1(nr)
real(8), intent(out) :: q0(nr),q1(nr)
! local variables
integer ir,ir0
real(8) ri,t1,t2,t3,t4
t1=1.d0/sol**2
t2=dble(l*(l+1))
! determine the r -> 0 boundary values of P and Q
ri=1.d0/r(1)
t3=2.d0+t1*(e-vr(1))
t4=t2/(t3*r(1)**2)+vr(1)-e
q0(1)=1.d0
q1(1)=0.d0
p0(1)=(q1(1)+q0(1)*ri)/t4
p1(1)=t3*q0(1)+p0(1)*ri
! extrapolate to the first four points
p1(2:4)=p1(1)
q1(2:4)=q1(1)
nn=0
do ir=2,nr
  ri=1.d0/r(ir)
  t3=2.d0+t1*(e-vr(ir))
  t4=t2/(t3*r(ir)**2)+vr(ir)-e
  ir0=ir-3
  if (ir0.lt.1) ir0=1
  p1(ir)=poly3(r(ir0),p1(ir0),r(ir))
  q1(ir)=poly3(r(ir0),q1(ir0),r(ir))
! integrate to find wavefunction
  p0(ir)=poly4i(r(ir0),p1(ir0),r(ir))+p0(ir0)
  q0(ir)=poly4i(r(ir0),q1(ir0),r(ir))+q0(ir0)
! compute the derivatives
  p1(ir)=t3*q0(ir)+p0(ir)*ri
  q1(ir)=t4*p0(ir)-q0(ir)*ri
! integrate for correction
  p0(ir)=poly4i(r(ir0),p1(ir0),r(ir))+p0(ir0)
  q0(ir)=poly4i(r(ir0),q1(ir0),r(ir))+q0(ir0)
! compute the derivatives again
  p1(ir)=t3*q0(ir)+p0(ir)*ri
  q1(ir)=t4*p0(ir)-q0(ir)*ri
! check for overflow
  if ((abs(p0(ir)).gt.1.d100).or.(abs(p1(ir)).gt.1.d100).or. &
      (abs(q0(ir)).gt.1.d100).or.(abs(q1(ir)).gt.1.d100)) then
    p0(ir:nr)=p0(ir)
    p1(ir:nr)=p1(ir)
    q0(ir:nr)=q0(ir)
    q1(ir:nr)=q1(ir)
    return
  end if
! check for node
  if (p0(ir-1)*p0(ir).lt.0.d0) nn=nn+1
end do
return

contains

real(8) function poly3(xa,ya,x)
implicit none
! arguments
real(8) xa(3),ya(3),x
! local variables
real(8) x0,x1,x2,y0,y1,y2
real(8) c1,c2,t0,t1,t2
! evaluate the polynomial coefficients
x0=xa(1)
x1=xa(2)-x0
x2=xa(3)-x0
y0=ya(1)
y1=ya(2)-y0
y2=ya(3)-y0
t0=1.d0/(x1*x2*(x2-x1))
t1=x1*y2
t2=x2*y1
c1=x2*t2-x1*t1
c2=t1-t2
t1=x-x0
! evaluate the polynomial
poly3=y0+t0*t1*(c1+c2*t1)
return
end function

real(8) function poly4i(xa,ya,x)
implicit none
! arguments
real(8), intent(in) :: xa(4),ya(4),x
! local variables
real(8) x0,x1,x2,x3,y0,y1,y2,y3
real(8) c1,c2,c3,t0,t1,t2,t3,t4,t5,t6
! evaluate the polynomial coefficients
x0=xa(1)
x1=xa(2)-x0
x2=xa(3)-x0
x3=xa(4)-x0
y0=ya(1)
y1=ya(2)-y0
y2=ya(3)-y0
y3=ya(4)-y0
t0=1.d0/(x1*x2*x3*(x1-x2)*(x1-x3)*(x2-x3))
t1=x1*x2*y3
t2=x2*x3*y1
t3=x3*x1*y2
c3=t1*(x1-x2)+t2*(x2-x3)+t3*(x3-x1)
t6=x3**2
t5=x2**2
t4=x1**2
c2=t1*(t5-t4)+t2*(t6-t5)+t3*(t4-t6)
c1=t1*(x2*t4-x1*t5)+t2*(x3*t5-x2*t6)+t3*(x1*t6-x3*t4)
t1=x-x0
! integrate the polynomial
poly4i=t1*(y0+t0*t1*(0.5d0*c1+t1*(0.3333333333333333333d0*c2+0.25d0*c3*t1)))
return
end function

end subroutine
!EOC
