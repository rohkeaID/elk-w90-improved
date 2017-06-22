
! Copyright (C) 2010 T. Bjorkman and O. Granas.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: findswidth
! !INTERFACE:
subroutine findswidth
! !USES:
use modmain
! !DESCRIPTION:
!   Calculates the smearing width from the $k$-point density,
!   $V_{\rm BZ}/n_k$; the valence band width, $W$; and an effective mass
!   parameter, $m^{*}$; according to
!   $$ \sigma=\frac{\sqrt{2W}}{m^{*}}\left(\frac{3}{4\pi}
!    \frac{V_{\rm BZ}}{n_k}\right)^{1/3}. $$
!   The valence bandwidth is determined by stepping down in energy from the
!   Fermi level until a gap larger than a given tolerance is found. This method
!   was presented in T. Bj\"{o}rkman and O. Gr\aa n\"{a}s, {\it Int. J. Quant.
!   Chem.} DOI: 10.1002/qua.22476.
!
! !REVISION HISTORY:
!   Created April 2010 (Torbjorn Bjorkman and JKD)
!EOP
!BOC
implicit none
! local variables
integer i,j,m,n,ik,ist
real(8), parameter :: eps=0.05d0
real(8) e,vbw
! allocatable arrays
integer, allocatable :: idx(:)
n=nstsv*nkpt
allocate(idx(n))
! find the index which sorts the eigenvalues in ascending order
call sortidx(n,evalsv,idx)
! find the highest eigenvalue < efermi
m=n
e=efermi
do i=n,1,-1
  j=idx(i)
  ik=(j-1)/nstsv+1
  ist=mod(j-1,nstsv)+1
  if (evalsv(ist,ik).lt.efermi) then
    m=i
    e=evalsv(ist,ik)
    goto 10
  end if
end do
10 continue
! search downwards until a gap larger than eps is found
do i=m,1,-1
  j=idx(i)
  ik=(j-1)/nstsv+1
  ist=mod(j-1,nstsv)+1
  if ((e-evalsv(ist,ik)).gt.eps) goto 20
  e=evalsv(ist,ik)
end do
20 continue
! valence band width
vbw=efermi-e
vbw=max(vbw,1.d-2)
! determine swidth
swidth=(sqrt(2.d0*vbw)/mstar)*(6.d0*pi**2/(omega*dble(nkptnr)))**(1.d0/3.d0)
deallocate(idx)
return
end subroutine
!EOC

