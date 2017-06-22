
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modfxcifc

use libxcifc

contains

subroutine fxcifc(fxctype,n,rho,rhoup,rhodn,fxc,fxcuu,fxcud,fxcdd)
implicit none
! mandatory arguments
integer, intent(in) :: fxctype(3),n
! optional arguments
real(8), optional, intent(in) :: rho(n),rhoup(n),rhodn(n)
real(8), optional, intent(out) :: fxc(n),fxcuu(n),fxcud(n),fxcdd(n)
! allocatable arrays
real(8), allocatable :: ra(:,:)
if (n.le.0) then
  write(*,*)
  write(*,'("Error(fxcifc): n <= 0 : ",I8)') n
  write(*,*)
  stop
end if
select case(abs(fxctype(1)))
case(0,1)
! f_xc = 0
  if (present(fxcuu).and.present(fxcud).and.present(fxcdd)) then
    fxcuu(:)=0.d0
    fxcud(:)=0.d0
    fxcdd(:)=0.d0
  else if (present(fxc)) then
    fxc(:)=0.d0
  else
    goto 10
  end if
case(3)
! Perdew-Wang-Ceperley-Alder
  if (present(rhoup).and.present(rhodn).and.present(fxcuu).and.present(fxcud) &
   .and.present(fxcdd)) then
! spin-polarised density
    call fxc_pwca(n,rhoup,rhodn,fxcuu,fxcud,fxcdd)
  else if (present(rho).and.present(fxc)) then
! divide spin-unpolarised density into up and down
    allocate(ra(n,4))
    ra(:,1)=0.5d0*rho(:)
    call fxc_pwca(n,ra(:,1),ra(:,1),ra(:,2),ra(:,3),ra(:,4))
    fxc(:)=0.5d0*(ra(:,2)+ra(:,3))
    deallocate(ra)
  else
    goto 10
  end if
case(100)
! libxc library functionals
  if (present(rhoup).and.present(rhodn).and.present(fxcuu).and.present(fxcud) &
   .and.present(fxcdd)) then
! LSDA
    call fxcifc_libxc(fxctype,n,rhoup=rhoup,rhodn=rhodn,fxcuu=fxcuu, &
     fxcud=fxcud,fxcdd=fxcdd)
  else if (present(rho).and.present(fxc)) then
! LDA
    call fxcifc_libxc(fxctype,n,rho=rho,fxc=fxc)
  else
    goto 10
  end if
case default
  write(*,*)
  write(*,'("Error(fxcifc): response function unavailable for fxctype ",3I8)') &
   fxctype
  write(*,*)
  stop
end select
return
10 continue
write(*,*)
write(*,'("Error(fxcifc): missing arguments for exchange-correlation type ",&
 &3I6)') fxctype(:)
write(*,*)
stop
end subroutine

end module

