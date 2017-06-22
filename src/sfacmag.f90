
! Copyright (C) 2010 A. I. Baranov and F. Wagner.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: sfacmag
! !INTERFACE:
subroutine sfacmag
! !USES:
use modmain
use modpw
! !DESCRIPTION:
!   Outputs magnetic structure factors, i.e. the Fourier transform coefficients
!   of each component $j$ of magnetization ${\bf m}({\bf r})$,
!   $$ F_j({\bf H})=\int_{\Omega}d^3r\,m_j({\bf r})e^{i{\bf H}\cdot{\bf r}}, $$
!   to the files {\tt SFACMAG\_j.OUT}. The lattice coordinates $(h,k,l)$ of
!   $\bf H$-vectors in this file are transformed by the matrix {\tt vhmat}. See
!   also routines {\tt zftrf} and {\tt genhvec}.
!
! !REVISION HISTORY:
!   Created July 2010 (Alexey I. Baranov)
!   Added multiplicity of the H-vectors, Oct. 2010 (Alexey I. Baranov)
!EOP
!BOC
implicit none
! local variables
integer idm,ih,iv(3)
real(8) v(3),a,b,r
character(256) fname
! allocatable arrays
complex(8), allocatable :: zmagh(:)
if (.not.spinpol) return
! initialise the structure factor specific variables
call sfacinit
! generate the magnetic structure factors
allocate(zmagh(nhvec))
do idm=1,ndmag
  call zftrf(nhvec,ivh,vhc,magmt(:,:,:,idm),magir(:,idm),zmagh)
  write(fname,'("SFACMAG_",I1.1,".OUT")') idm
  open(50,file=trim(fname),action='WRITE',form='FORMATTED')
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
    a=dble(zmagh(ih))*omega
    b=-aimag(zmagh(ih))*omega
    r=abs(zmagh(ih))*omega
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
  close(50)
end do
deallocate(zmagh)
write(*,*)
write(*,'("Info(sfacmag): magnetic structure factors written to &
 &SFACMAG_j.OUT")')
write(*,'(" for magnetic components j = ",3I2)') (idm,idm=1,ndmag)
if (ndmag.eq.1) then
  write(*,'(" (this corresponds to the z-component of the magnetisation)")')
end if
write(*,*)
write(*,'(" Energy window : ",2G18.10)') wsfac(:)
return
end subroutine
!EOC

