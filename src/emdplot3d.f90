
! Copyright (C) 2015 D. Ernsting, S. Dugdale and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine emdplot3d(emds)
use modmain
use modpw
implicit none
! arguments
real(4), intent(in) :: emds(nhkmax,nkpt)
! local variables
integer np,ip
real(8) v1(3),t1
! allocatable arrays
real(8), allocatable :: vpl(:,:)
! external functions
real(8) rfhkintp
external rfhkintp
! total number of plot points
np=np3d(1)*np3d(2)*np3d(3)
! generate the 3D plotting points
allocate(vpl(3,np))
call plotpt3d(vpl)
open(50,file='EMD3D.OUT',action='WRITE',form='FORMATTED')
write(50,'(3I6," : grid size")') np3d(:)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(t1,v1)
!$OMP DO ORDERED
do ip=1,np
  t1=rfhkintp(vpl(:,ip),emds)
  call r3mv(bvec,vpl(:,ip),v1)
!$OMP ORDERED
  write(50,'(4G18.10)') v1(:),t1
!$OMP END ORDERED
end do
!$OMP END DO
!$OMP END PARALLEL
close(50)
deallocate(vpl)
return
end subroutine

