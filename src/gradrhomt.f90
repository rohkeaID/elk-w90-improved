
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gradrhomt
use modmain
use modphonon
implicit none
! local variables
integer nr,nri
! allocatable arrays
real(8), allocatable :: rfmt(:,:),grfmt(:,:,:)
! add gradient contribution from rigid shift of muffin-tin
allocate(rfmt(lmmaxvr,nrmtmax),grfmt(lmmaxvr,nrmtmax,3))
nr=nrmt(isph)
nri=nrmtinr(isph)
call gradrfmt(nr,nri,rsp(:,isph),rhomt(:,:,iasph),nrmtmax,grfmt)
! convert to spherical coordinates
call rbsht(nr,nri,1,grfmt(:,:,ipph),1,rfmt)
! subtract from density derivative
drhomt(:,1:nr,iasph)=drhomt(:,1:nr,iasph)-rfmt(:,1:nr)
deallocate(rfmt,grfmt)
return
end subroutine

