
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: gengclq
! !INTERFACE:
subroutine gengclq
! !USES:
use modmain
use modtest
! !DESCRIPTION:
!   The Fock matrix elements
!   $$ V_{ij{\bf k}}\equiv\sum_{l{\bf k'}}\int
!    \frac{\Psi^{\dag}_{i{\bf k}}({\bf r})\cdot\Psi_{l{\bf k}'}({\bf r})
!    \Psi^{\dag}_{l{\bf k}'}({\bf r}')\cdot\Psi_{j{\bf k}}({\bf r}')}
!    {|{\bf r}-{\bf r'}|}\,d^3r\,d^3r' $$
!   contain a divergent term in the sum over ${\bf k}'$ which behaves as
!   $1/q^2$, where ${\bf q}\equiv{\bf k}-{\bf k}'$ is in the first Brillouin
!   zone. The resulting convergence with respect to the number of discrete
!   $q$-points, $N_q$, is very slow. This routine computes the regularised
!   Coulomb Green's function
!   \begin{align}
!    g({\bf q}_i)=\frac{4\pi}{V}\int_{V_i}\frac{1}{q^2}\,d^3q,
!   \end{align}
!   where the integral is over the small parallelepiped with volume
!   $V=\Omega_{\rm BZ}/N_q$ and centered on the discrete point ${\bf q}_i$.
!   This dramatically increases the rate of convergence of methods which involve
!   a summation over the $1/q^2$ part of the Coulomb interaction. The above
!   integral is evaluated numerically on increasingly finer grids and then
!   extrapolated to the continuum.
!
! !REVISION HISTORY:
!   Created August 2004 (JKD,SS)
!   Changed from genwiq2, July 2017 (JKD)
!EOP
!BOC
implicit none
! local variables
integer, parameter :: np=5
integer, parameter :: ns0=10,nss=20
integer ns,iq,i1,i2,i3,ip
real(8) d(3),sum,t1,t2
real(8) v1(3),v2(3),v3(3)
real(8) xa(np),ya(np),c(np)
! external functions
real(8) polynom
external polynom
! allocate global gclq array
if (allocated(gclq)) deallocate(gclq)
allocate(gclq(nqpt))
! begin loop over q-points, note that the vectors vqc are assumed to be in the
! first Brillouin zone
do iq=1,nqpt
! loop over different subdivisions
  ns=ns0
  do ip=1,np
! subdivision vectors in lattice coordinates
    d(:)=1.d0/dble(ngridq(:)*2*ns)
! compute the integral of 1/q^2
    sum=0.d0
    do i1=-ns,ns-1
      t1=dble(i1)*d(1)
      v1(:)=vqc(:,iq)+t1*bvec(:,1)
      do i2=-ns,ns-1
        t1=dble(i2)*d(2)
        v2(:)=v1(:)+t1*bvec(:,2)
        do i3=-ns,ns-1
          t1=dble(i3)*d(3)
          v3(:)=v2(:)+t1*bvec(:,3)
          t2=v3(1)**2+v3(2)**2+v3(3)**2
          if (t2.gt.1.d-14) sum=sum+1.d0/t2
        end do
      end do
    end do
    t1=1.d0/dble(2*ns)
    xa(ip)=t1
    ya(ip)=fourpi*sum*t1**3
! increment number of subdivisions
    ns=ns+nss
  end do
! extrapolate the volume element to zero with a polynomial
  gclq(iq)=polynom(0,np,xa,ya,c,0.d0)
end do
! write gclq to test file
call writetest(800,"regularised Coulomb Green''s function (gclq)",nv=nqpt, &
 tol=1.d-8,rva=gclq)
return
end subroutine
!EOC

