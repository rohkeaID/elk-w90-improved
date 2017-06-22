
! Copyright (C) 2012 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvfxcg(gqc,vfxc)
use modmain
implicit none
! arguments
real(8), intent(in) :: gqc(ngrf)
complex(8), intent(out) :: vfxc(ngrf,ngrf,nwrf)
! local variables
integer ig,jg,kg,iv(3)
real(8) t0
complex(8) z1
! allocatable arrays
real(8), allocatable :: fxcmt(:,:,:),fxcir(:)
complex(8), allocatable :: fxcg(:)
allocate(fxcmt(lmmaxvr,nrmtmax,natmtot),fxcir(ngtot))
allocate(fxcg(ngtot))
! generate the kernel f_xc in real-space
call genfxcr(.true.,fxcmt,fxcir)
! Fourier transform the kernel to G-space
call zftrf(ngtot,ivg,vgc,fxcmt,fxcir,fxcg)
t0=1.d0/fourpi
do ig=1,ngrf
  do jg=1,ngrf
    iv(:)=ivg(:,ig)-ivg(:,jg)
    if ((iv(1).ge.intgv(1,1)).and.(iv(1).le.intgv(2,1)).and. &
        (iv(2).ge.intgv(1,2)).and.(iv(2).le.intgv(2,2)).and. &
        (iv(3).ge.intgv(1,3)).and.(iv(3).le.intgv(2,3))) then
      kg=ivgig(iv(1),iv(2),iv(3))
      z1=t0*fxcg(kg)*(gqc(ig)*gqc(jg))
      vfxc(ig,jg,:)=z1
    end if
  end do
end do
deallocate(fxcmt,fxcir,fxcg)
return
end subroutine

