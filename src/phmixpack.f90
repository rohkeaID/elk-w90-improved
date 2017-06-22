
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phmixpack(tpack,n,v)
use modmain
use modphonon
implicit none
! arguments
logical, intent(in) :: tpack
integer, intent(out) :: n
real(8), intent(inout) :: v(*)
! local variables
integer idm
n=0
call zfpack(tpack,n,nrmt,nrmtinr,nrmtmax,dvsmt,dvsir,v)
do idm=1,ndmag
  call zfpack(tpack,n,nrcmt,nrcmtinr,nrcmtmax,dbsmt(:,:,:,idm),dbsir(:,idm),v)
end do
return
end subroutine

