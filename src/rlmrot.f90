
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rlmrot
! !INTERFACE:
subroutine rlmrot(p,ang,lmax,ld,d)
! !INPUT/OUTPUT PARAMETERS:
!   p    : if p=-1 then the rotation matrix is improper (in,integer)
!   ang  : Euler angles; alpha, beta, gamma (in,real(3))
!   lmax : maximum angular momentum (in,integer)
!   ld   : leading dimension (in,integer)
!   d    : real spherical harmonic rotation matrix (out,real(ld,*))
! !DESCRIPTION:
!   Returns the rotation matrix in the basis of real spherical harmonics given
!   the three Euler angles, $(\alpha,\beta,\gamma)$, and the parity, $p$, of the
!   rotation. The matrix is determined using the formula of V. V. Nechaev,
!   [{\it J. Struct. Chem.} {\bf 35}, 115 (1994)], suitably modified for our
!   definition of the real spherical harmonics ($m_1>0$, $m_2>0$):
!   \begin{align*}
!    &\Delta^l_{00}=d^l_{00}, \\
!    &\Delta^l_{m_10}=\sqrt{2}\,(-1)^{m_1}d^l_{0m_1}\cos(m_1\alpha), \\
!    &\Delta^l_{0m_2}=\sqrt{2}\,(-1)^{m_2}d^l_{m_20}\cos(m_2\gamma), \\
!    &\Delta^l_{-m_10}=-\sqrt{2}\,d^l_{0m_1}\sin(m_1\alpha), \\
!    &\Delta^l_{0-m_2}=\sqrt{2}\,d^l_{m_20}\sin(m_2\gamma), \\
!    &\Delta^l_{m_1m_2}=(-1)^{m_1}(-1)^{m_2}\{\cos(m_1\alpha)\cos(m_2\gamma)
!     [d_A+d_B]-\sin(m_1\alpha)\sin(m_2\gamma)[d_A-d_B]\}, \\
!    &\Delta^l_{m_1-m_2}=(-1)^{m_1}\{\sin(m_1\alpha)\cos(m_2\gamma)
!     [d_A-d_B]+\cos(m_1\alpha)\sin(m_2\gamma)[d_A+d_B]\}, \\
!    &\Delta^l_{-m_1m_2}=-(-1)^{m_2}\{\sin(m_1\alpha)\cos(m_2\gamma)
!     [d_A+d_B]+\cos(m_1\alpha)\sin(m_2\gamma)[d_A-d_B]\}, \\
!    &\Delta^l_{-m_1-m_2}=\cos(m_1\alpha)\cos(m_2\gamma)
!     [d_A-d_B]-\sin(m_1\alpha)\sin(m_2\gamma)[d_A+d_B],
!   \end{align*}
!   where $d_A\equiv d^l_{-m_1-m_2}$, $d_B\equiv(-1)^{m_1}d^l_{m_1-m_2}$ and
!   $d$ is the rotation matrix about the $y$-axis for complex spherical
!   harmonics. See the routines {\tt genrlm}, {\tt roteuler} and {\tt ylmroty}.
!
! !REVISION HISTORY:
!   Created December 2008 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: p
real(8), intent(in) :: ang(3)
integer, intent(in) :: lmax
integer, intent(in) :: ld
real(8), intent(out) :: d(ld,*)
! local variables
integer lmmax,l,m1,m2,lm,lm0
real(8), parameter :: sqtwo=1.4142135623730950488d0
real(8) s1,s2,t1,t2,t3,t4,t5,t6,t7,t8
! automatic arrays
integer lmi(-lmax:lmax)
real(8) ca(lmax),sa(lmax),cg(lmax),sg(lmax)
! allocatable arrays
real(8), allocatable :: dy(:,:)
if (lmax.lt.0) then
  write(*,*)
  write(*,'("Error(rlmrot): lmax < 0 : ",I8)') lmax
  write(*,*)
  stop
end if
lmmax=(lmax+1)**2
allocate(dy(lmmax,lmmax))
! generate the complex spherical harmonic rotation matrix about the y-axis
call ylmroty(ang(2),lmax,lmmax,dy)
do m1=1,lmax
  ca(m1)=cos(dble(m1)*ang(1))
  sa(m1)=sin(dble(m1)*ang(1))
  cg(m1)=cos(dble(m1)*ang(3))
  sg(m1)=sin(dble(m1)*ang(3))
end do
lm=0
do l=0,lmax
  do m1=-l,l
    lm=lm+1
    lmi(m1)=lm
  end do
  lm0=lmi(0)
  d(lm0,lm0)=dy(lm0,lm0)
  do m1=1,l
    s1=1.d0
    if (mod(m1,2).ne.0) s1=-1.d0
    t1=sqtwo*dy(lm0,lmi(m1))
    t2=sqtwo*dy(lmi(m1),lm0)
    d(lmi(m1),lm0)=s1*t1*ca(m1)
    d(lm0,lmi(m1))=s1*t2*cg(m1)
    d(lmi(-m1),lm0)=-t1*sa(m1)
    d(lm0,lmi(-m1))=t2*sg(m1)
    do m2=1,l
      s2=1.d0
      if (mod(m2,2).ne.0) s2=-1.d0
      t1=ca(m1)*cg(m2)
      t2=sa(m1)*sg(m2)
      t3=sa(m1)*cg(m2)
      t4=ca(m1)*sg(m2)
      t5=dy(lmi(-m1),lmi(-m2))
      t6=s1*dy(lmi(m1),lmi(-m2))
      t7=t5+t6
      t8=t5-t6
      d(lmi(m1),lmi(m2))=s1*s2*(t1*t7-t2*t8)
      d(lmi(m1),lmi(-m2))=s1*(t3*t8+t4*t7)
      d(lmi(-m1),lmi(m2))=-s2*(t3*t7+t4*t8)
      d(lmi(-m1),lmi(-m2))=t1*t8-t2*t7
    end do
  end do
  if ((p.eq.-1).and.(mod(l,2).ne.0)) then
    do m1=-l,l
      do m2=-l,l
        d(lmi(m1),lmi(m2))=-d(lmi(m1),lmi(m2))
      end do
    end do
  end if
end do
deallocate(dy)
return
end subroutine
!EOC

