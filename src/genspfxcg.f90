
! Copyright (C) 2013 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genspfxcg(fxc)
use modmain
implicit none
! arguments
complex(8), intent(out) :: fxc(ngrf,4,ngrf,4)
! local variables
integer ig,jg,kg
integer iv(3),i,j
complex(8) z1
! allocatable arrays
real(8), allocatable :: fxcmt(:,:,:,:,:),fxcir(:,:,:)
complex(8), allocatable :: fxcg(:)
allocate(fxcmt(lmmaxvr,nrmtmax,natmtot,4,4),fxcir(ngtot,4,4))
allocate(fxcg(ngtot))
! generate the kernel f_xc in real-space
call genspfxcr(.true.,fxcmt,fxcir)
! Fourier transform the kernel to G-space
do i=1,4
  do j=i,4
    call zftrf(ngtot,ivg,vgc,fxcmt(:,:,:,i,j),fxcir(:,i,j),fxcg)
    do ig=1,ngrf
      do jg=1,ngrf
        iv(:)=ivg(:,ig)-ivg(:,jg)
        if ((iv(1).ge.intgv(1,1)).and.(iv(1).le.intgv(2,1)).and. &
            (iv(2).ge.intgv(1,2)).and.(iv(2).le.intgv(2,2)).and. &
            (iv(3).ge.intgv(1,3)).and.(iv(3).le.intgv(2,3))) then
          kg=ivgig(iv(1),iv(2),iv(3))
          z1=fxcg(kg)
          fxc(ig,i,jg,j)=z1
          fxc(jg,j,ig,i)=conjg(z1)
        end if
      end do
    end do
  end do
end do
deallocate(fxcmt,fxcir,fxcg)
return
end subroutine

