
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gendvsig
use modmain
use modphonon
implicit none
! local variables
integer ig,ifg
! allocatable arrays
complex(8), allocatable :: zfft(:)
allocate(zfft(ngtot))
zfft(:)=dvsir(:)*cfunir(:)+vsir(:)*dcfunir(:)
! Fourier transform to G+q-space
call zfftifc(3,ngridg,-1,zfft)
! store in global array
do ig=1,ngvec
  ifg=igfft(ig)
  dvsig(ig)=zfft(ifg)
end do
deallocate(zfft)
return
end subroutine

