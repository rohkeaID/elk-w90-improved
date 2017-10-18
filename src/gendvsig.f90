
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gendvsig
use modmain
use modphonon
implicit none
! local variables
integer ig
! allocatable arrays
complex(8), allocatable :: zfft1(:),zfft2(:)
allocate(zfft1(ngtot))
if (trimvg) then
! trim the Fourier components |G| > 3*gkmax
  zfft1(:)=dvsir(:)
  call zfftifc(3,ngridg,-1,zfft1)
  do ig=ng3gk+1,ngtot
    zfft1(igfft(ig))=0.d0
  end do
  call zfftifc(3,ngridg,1,zfft1)
  allocate(zfft2(ngtot))
  zfft2(:)=vsir(:)
  call zfftifc(3,ngridg,-1,zfft2)
  do ig=ng3gk+1,ngtot
    zfft2(igfft(ig))=0.d0
  end do
  call zfftifc(3,ngridg,1,zfft2)
  zfft1(:)=zfft1(:)*cfunir(:)+zfft2(:)*dcfunir(:)
  deallocate(zfft2)
else
  zfft1(:)=dvsir(:)*cfunir(:)+vsir(:)*dcfunir(:)
end if
! Fourier transform to G+q-space
call zfftifc(3,ngridg,-1,zfft1)
! store in global array
do ig=1,ngvec
  dvsig(ig)=zfft1(igfft(ig))
end do
deallocate(zfft1)
return
end subroutine

