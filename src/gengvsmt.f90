
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gengvsmt
use modmain
use modphonon
implicit none
! local variables
integer nr,nri,np
! allocatable arrays
complex(8), allocatable :: zfmt(:),gzfmt(:,:)
allocate(zfmt(npmtmax),gzfmt(npmtmax,3))
nr=nrmt(isph)
nri=nrmti(isph)
np=npmt(isph)
! convert potential to complex spherical harmonics
call rtozfmt(nr,nri,vsmt(:,iasph),zfmt)
! calculate the gradient
call gradzfmt(nr,nri,rsp(:,isph),zfmt,npmtmax,gzfmt)
! copy current polarisation component to global array
gvsmt(1:np)=gzfmt(1:np,ipph)
deallocate(zfmt,gzfmt)
return
end subroutine

