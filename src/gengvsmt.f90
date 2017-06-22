
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gengvsmt
use modmain
use modphonon
implicit none
! local variables
integer nr,nri
! allocatable arrays
complex(8), allocatable :: zfmt(:,:),gzfmt(:,:,:)
allocate(zfmt(lmmaxvr,nrmtmax),gzfmt(lmmaxvr,nrmtmax,3))
nr=nrmt(isph)
nri=nrmtinr(isph)
! convert potential to complex spherical harmonics
call rtozfmt(nr,nri,1,vsmt(:,:,iasph),1,zfmt)
! calculate the gradient
call gradzfmt(nr,nri,rsp(:,isph),zfmt,nrmtmax,gzfmt)
! copy current polarisation component to global array
gvsmt(:,1:nr)=gzfmt(:,1:nr,ipph)
deallocate(zfmt,gzfmt)
return
end subroutine

