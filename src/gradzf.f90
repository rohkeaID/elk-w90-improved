
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gradzf(zfmt,zfir,gzfmt,gzfir)
use modmain
implicit none
! arguments
complex(8), intent(in) :: zfmt(lmmaxvr,nrmtmax,natmtot),zfir(ngtot)
complex(8), intent(out) :: gzfmt(lmmaxvr,nrmtmax,natmtot,3),gzfir(ngtot,3)
! local variables
integer is,ias,nr
integer ig,ifg,i
complex(8) z1
! allocatable arrays
complex(8), allocatable :: gzfmt1(:,:,:),zfft(:)
! muffin-tin gradient
allocate(gzfmt1(lmmaxvr,nrmtmax,3))
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  call gradzfmt(nr,nrmtinr(is),rsp(:,is),zfmt(:,:,ias),nrmtmax,gzfmt1)
  do i=1,3
    gzfmt(:,1:nr,ias,i)=gzfmt1(:,1:nr,i)
  end do
end do
deallocate(gzfmt1)
! interstitial gradient
allocate(zfft(ngtot))
call zcopy(ngtot,zfir,1,zfft,1)
call zfftifc(3,ngridg,-1,zfft)
do i=1,3
  gzfir(:,i)=0.d0
  do ig=1,ngvec
    ifg=igfft(ig)
    z1=zfft(ifg)
    gzfir(ifg,i)=vgc(i,ig)*cmplx(-aimag(z1),dble(z1),8)
  end do
  call zfftifc(3,ngridg,1,gzfir(:,i))
end do
deallocate(zfft)
return
end subroutine

