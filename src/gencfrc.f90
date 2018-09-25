
! Copyright (C) 2017 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gencfrc
use modmain
implicit none
! allocatable arrays
complex(8), allocatable :: zfft(:)
if (allocated(cfrc)) deallocate(cfrc)
allocate(cfrc(ngtc))
allocate(zfft(ngtc))
zfft(:)=0.d0
zfft(igfc(1:ngvc))=cfunig(1:ngvc)
! Fourier transform to real-space
call zfftifc(3,ngdc,1,zfft)
cfrc(:)=dble(zfft(:))
deallocate(zfft)
return
end subroutine

