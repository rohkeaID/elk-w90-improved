
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genachi(evecbdg,achi)
use modmain
use modscdft
implicit none
! arguments
complex(8), intent(in) :: evecbdg(nmbdg,nmbdg)
complex(8), intent(out) :: achi(nbdg,nbdg)
! local variables
integer i,j,k
complex(8) z1
! zero the anomalous density
achi(:,:)=0.d0
! loop over the BdG states
do k=1,nmbdg
  if (evalbdg(k).lt.0.d0) then
    do i=1,nbdg
      z1=conjg(evecbdg(i,k))
      do j=1,nbdg
        achi(i,j)=achi(i,j)+z1*evecbdg(nbdg+j,k)
      end do
    end do
  else
    do i=1,nbdg
      z1=evecbdg(nbdg+i,k)
      do j=1,nbdg
        achi(i,j)=achi(i,j)+z1*conjg(evecbdg(j,k))
      end do
    end do
  end if
end do
return
end subroutine

