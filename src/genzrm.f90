
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genzrm(n,wf11,wf12,wf21,wf22,zrho,ld,zmag)
use modmain
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(in) :: wf11(n),wf12(n),wf21(n),wf22(n)
complex(8), intent(out) :: zrho(n)
integer, intent(in) :: ld
complex(8), intent(out) :: zmag(ld,ndmag)
! local variables
integer i
complex(8) z11,z12,z21,z22,z1,z2
if (ncmag) then
! non-collinear case
  do i=1,n
    z11=wf11(i)
    z12=wf12(i)
    z21=wf21(i)
    z22=wf22(i)
! up-dn spin density
    z1=conjg(z11)*z22
! dn-up spin density
    z2=conjg(z12)*z21
! x-component: up-dn + dn-up
    zmag(i,1)=z1+z2
! y-component: i*(dn-up - up-dn)
    z1=z2-z1
    zmag(i,2)=cmplx(-aimag(z1),dble(z1),8)
    z1=conjg(z11)*z21
    z2=conjg(z12)*z22
! z-component: up-up - dn-dn
    zmag(i,3)=z1-z2
! density: up-up + dn-dn
    zrho(i)=z1+z2
  end do
else
! collinear case
  do i=1,n
    z1=conjg(wf11(i))*wf21(i)
    z2=conjg(wf12(i))*wf22(i)
    zmag(i,1)=z1-z2
    zrho(i)=z1+z2
  end do
end if
return
end subroutine

