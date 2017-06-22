
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: mixadapt
! !INTERFACE:
subroutine mixadapt(iscl,beta0,betamax,n,nu,mu,beta,f,d)
! !INPUT/OUTPUT PARAMETERS:
!   iscl    : self-consistent loop number (in,integer)
!   beta0   : mixing parameter (in,real)
!   betamax : maximum mixing parameter (in,real)
!   n       : vector length (in,integer)
!   nu      : current output vector as well as the next input vector of the
!             self-consistent loop (inout,real(n))
!   mu      : used internally (inout,real(n))
!   beta    : used internally (inout,real(n))
!   f       : used internally (inout,real(n))
!   d       : RMS difference between old and new output vectors (out,real)
! !DESCRIPTION:
!   Given the input vector $\mu^i$ and output vector $\nu^i$ of the $i$th
!   self-consistent loop, this routine generates the next input vector to the
!   loop using an adaptive mixing scheme. The $j$th component of the output
!   vector is mixed with a fraction of the same component of the input vector:
!   $$ \mu^{i+1}_j=\beta^i_j\nu^i_j+(1-\beta^i_j)\mu^i_j, $$
!   where $\beta^{i+1}_j=\beta^i_j+\beta_0$ if $f^i_j\equiv\nu^i_j-\mu^i_j$ does
!   not change sign between loops. If $f^i_j$ does change sign, then
!   $\beta^{i+1}_j=(\beta^i_j+\beta_0)/2$. Note that the array {\tt nu} serves
!   for both input and output, and the arrays {\tt mu}, {\tt beta} and {\tt f}
!   are used internally and should not be changed between calls. The routine is
!   thread-safe so long as each thread has its own independent work arrays.
!   Complex arrays may be passed as real arrays with $n$ doubled.
!
! !REVISION HISTORY:
!   Created March 2003 (JKD)
!   Modified, September 2008 (JKD)
!   Modified, August 2011 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: iscl
real(8), intent(in) :: beta0,betamax
integer, intent(in) :: n
real(8), intent(inout) :: nu(n),mu(n)
real(8), intent(inout) :: beta(n),f(n)
real(8), intent(out) :: d
! local variables
integer i
real(8) t1
if (n.le.0) return
! initialise mixer
if (iscl.le.0) then
  call dcopy(n,nu,1,mu,1)
  f(:)=0.d0
  beta(:)=beta0
  d=1.d0
  return
end if
d=0.d0
do i=1,n
  t1=nu(i)-mu(i)
  if (t1*f(i).ge.0.d0) then
    beta(i)=beta(i)+beta0
    if (beta(i).gt.betamax) beta(i)=betamax
  else
    beta(i)=(beta(i)+beta0)*0.5d0
  end if
  f(i)=t1
  nu(i)=beta(i)*nu(i)+(1.d0-beta(i))*mu(i)
  d=d+t1**2
  mu(i)=nu(i)
end do
d=sqrt(d/dble(n))
return
end subroutine
!EOC

