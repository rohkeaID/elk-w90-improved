
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genstrain
use modmain
implicit none
! local variables
integer i,j,k
real(8) a(3,3),b(3,3),t1
! external functions
real(8) ddot,dnrm2
external ddot,dnrm2
nstrain=0
do i=1,3
  do j=1,3
! set strain tensor in lattice coordinates to delta_ij
    a(:,:)=0.d0
    a(i,j)=1.d0
! symmetrise strain tensor
    call symmat(a)
! convert to mixed Cartesian-lattice coordinates
    call r3mtm(ainv,a,b)
! orthogonalise strain tensor to previous tensors
    do k=1,nstrain
      t1=-ddot(9,strain(:,:,k),1,b,1)
      call daxpy(9,t1,strain(:,:,k),1,b,1)
    end do
! compute the norm
    t1=dnrm2(9,b,1)
    if (t1.lt.epslat) cycle
! normalise tensor and store in global array
    nstrain=nstrain+1
    strain(:,:,nstrain)=b(:,:)/t1
  end do
end do
! zero small components
do k=1,nstrain
  do i=1,3
    do j=1,3
      if (abs(strain(i,j,k)).lt.epslat) strain(i,j,k)=0.d0
    end do
  end do
end do
return
end subroutine

