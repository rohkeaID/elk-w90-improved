
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: ylmroty
! !INTERFACE:
subroutine ylmroty(beta,lmax,ld,dy)
! !INPUT/OUTPUT PARAMETERS:
!   beta : rotation angle about y-axis (in,real)
!   lmax : maximum angular momentum (in,integer)
!   ld   : leading dimension (in,integer)
!   dy   : rotation matrix for complex spherical harmonics (out,real(ld,*))
! !DESCRIPTION:
!   Returns the rotation matrix in the basis of complex spherical harmonics for
!   a rotation of angle $\beta$ about the $y$-axis. This matrix is real and is
!   given by the formula
!   \begin{align*}
!    d^l_{m_1m_2}(\beta)=&[(l+m_1)!(l-m_1)!(l+m_2)!(l-m_2)!]^{1/2}\\
!    &\times\sum_k(-1)^k\frac{\left(\cos\frac{\beta}{2}\right)^{2(l-k)-m_2+m_1}
!    \left(\sin\frac{\beta}{2}\right)^{2k+m_2-m_1}}
!    {k!(l+m_1-k)!(l-m_2-k)!(m_2-m_1+k)!},
!   \end{align*}
!   where $k$ runs through all integer values for which the factorials exist.
!
! !REVISION HISTORY:
!   Created December 2008 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: beta
integer, intent(in) :: lmax
integer, intent(in) :: ld
real(8), intent(out) :: dy(ld,*)
! local variables
integer j,k,l,m1,m2,lm1,lm2
real(8) cb,sb,sum,t1,t2
! external functions
real(8) factnm
external factnm
cb=cos(beta/2.d0)
sb=sin(beta/2.d0)
lm1=0
do l=0,lmax
! generate rotation operator for m-components of current l
  do m1=-l,l
    lm1=lm1+1
    lm2=l**2
    do m2=-l,l
      lm2=lm2+1
      sum=0.d0
      do k=0,min(l+m1,l-m2)
        if (((l+m1-k).ge.0).and.((l-m2-k).ge.0).and.((m2-m1+k).ge.0)) then
          j=2*(l-k)+m1-m2
          t1=1.d0
          if (j.ne.0) t1=t1*cb**j
          j=2*k+m2-m1
          if (j.ne.0) t1=t1*sb**j
          t2=t1/(factnm(k,1)*factnm(l+m1-k,1)*factnm(l-m2-k,1) &
           *factnm(m2-m1+k,1))
          if (mod(k,2).ne.0) t2=-t2
          sum=sum+t2
        end if
      end do
      t1=sqrt(factnm(l+m1,1)*factnm(l-m1,1)*factnm(l+m2,1)*factnm(l-m2,1))
      dy(lm1,lm2)=t1*sum
    end do
  end do
end do
return
end subroutine
!EOC

