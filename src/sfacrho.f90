
! Copyright (C) 2010 A. I. Baranov and F. Wagner.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: sfacrho
! !INTERFACE:
subroutine sfacrho
! !USES:
use modmain
use modpw
! !DESCRIPTION:
!   Outputs X-ray structure factors, i.e. the Fourier transform coefficients of
!   the total electron density
!   $$ F({\bf H})=\int_{\Omega}d^3r\,\rho({\bf r})e^{i{\bf H}\cdot{\bf r}}, $$
!   to the file {\tt SFACRHO.OUT}. The lattice coordinates $(h,k,l)$ of the
!   $\bf H$-vectors in this file are transformed by the matrix {\tt vhmat}. If
!   and energy window is set using the variable {\tt wsfac}, then only those
!   states within the window are used to compute the density. See also routines
!   {\tt zftrf} and {\tt genhvec}.
!
! !REVISION HISTORY:
!   Created July 2010 (Alexey I. Baranov)
!   Added multiplicity of the H-vectors, Oct. 2010 (Alexey I. Baranov)
!EOP
!BOC
implicit none
! local variables
integer ih,iv(3)
real(8) v(3),a,b,r
! allocatable arrays
complex(8), allocatable :: zrhoh(:)
! initialise the structure factor specific variables
call sfacinit
! calculate the density structure factors
allocate(zrhoh(nhvec))
call zftrf(nhvec,ivh,vhc,rhomt,rhoir,zrhoh)
open(50,file='SFACRHO.OUT',action='WRITE',form='FORMATTED')
write(50,*)
write(50,'("h k l indices transformed by vhmat matrix:")')
write(50,'(3G18.10)') vhmat(:,1)
write(50,'(3G18.10)') vhmat(:,2)
write(50,'(3G18.10)') vhmat(:,3)
write(50,*)
write(50,'("      h      k      l  multipl.   |H|            Re(F)&
 &            Im(F)           |F|")')
write(50,*)
do ih=1,nhvec
! apply transformation matrix
  v(:)=vhmat(:,1)*dble(ivh(1,ih)) &
      +vhmat(:,2)*dble(ivh(2,ih)) &
      +vhmat(:,3)*dble(ivh(3,ih))
! in crystallography the forward Fourier transform of real-space density is
! usually done with positive phase and without 1/omega prefactor
  a=dble(zrhoh(ih))*omega
  b=-aimag(zrhoh(ih))*omega
  r=abs(zrhoh(ih))*omega
  iv(:)=nint(v(:))
  if ((abs(v(1)-iv(1)).le.epslat).and. &
      (abs(v(2)-iv(2)).le.epslat).and. &
      (abs(v(3)-iv(3)).le.epslat)) then
! integer hkl
    write(50,'(4I7,4G16.8)') iv(:),mulh(ih),hc(ih),a,b,r
  else
! non-integer hkl
    write(50,'(3F7.2,I7,4G16.8)') v(:),mulh(ih),hc(ih),a,b,r
  end if
end do
deallocate(zrhoh)
close(50)
write(*,*)
write(*,'("Info(sfacrho): density structure factors written to SFACRHO.OUT")')
write(*,*)
write(*,'(" Energy window : ",2G18.10)') wsfac(:)
return
end subroutine
!EOC

