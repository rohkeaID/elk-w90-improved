
! Copyright (C) 2009 T. McQueen and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

! Stub routines for libxc. See Elk manual for libxc installation instructions.

module libxcifc

integer libxcv(3)

contains

subroutine xcifc_libxc(xctype,n,c_tb09,rho,rhoup,rhodn,g2rho,g2up,g2dn,grho2, &
 gup2,gdn2,gupdn,tau,tauup,taudn,ex,ec,vx,vc,vxup,vxdn,vcup,vcdn,dxdgr2, &
 dxdgu2,dxdgd2,dxdgud,dcdgr2,dcdgu2,dcdgd2,dcdgud,dxdg2r,dxdg2u,dxdg2d,dcdg2r, &
 dcdg2u,dcdg2d,wx,wxup,wxdn,wc,wcup,wcdn)
implicit none
! mandatory arguments
integer, intent(in) :: xctype(3),n
! optional arguments
real(8), optional :: c_tb09
real(8), optional :: rho(n),rhoup(n),rhodn(n)
real(8), optional :: g2rho(n),g2up(n),g2dn(n)
real(8), optional :: grho2(n),gup2(n),gdn2(n),gupdn(n)
real(8), optional :: tau(n),tauup(n),taudn(n)
real(8), optional :: ex(n),ec(n),vx(n),vc(n)
real(8), optional :: vxup(n),vxdn(n),vcup(n),vcdn(n)
real(8), optional :: dxdgr2(n),dxdgu2(n),dxdgd2(n),dxdgud(n)
real(8), optional :: dcdgr2(n),dcdgu2(n),dcdgd2(n),dcdgud(n)
real(8), optional :: dxdg2r(n),dxdg2u(n),dxdg2d(n)
real(8), optional :: dcdg2r(n),dcdg2u(n),dcdg2d(n)
real(8), optional :: wx(n),wxup(n),wxdn(n)
real(8), optional :: wc(n),wcup(n),wcdn(n)
write(*,*)
write(*,'("Error(libxcifc): libxc not or improperly installed")')
write(*,*)
stop
end subroutine

subroutine fxcifc_libxc(fxctype,n,rho,rhoup,rhodn,fxc,fxcuu,fxcud,fxcdd)
implicit none
! mandatory arguments
integer, intent(in) :: fxctype(3),n
! optional arguments
real(8), optional :: rho(n),rhoup(n),rhodn(n)
real(8), optional :: fxc(n),fxcuu(n),fxcud(n),fxcdd(n)
write(*,*)
write(*,'("Error(libxcifc): libxc not or improperly installed")')
write(*,*)
stop
end subroutine

subroutine xcdata_libxc(xctype,xcdescr,xcspin,xcgrad,hybrid,hybridc)
implicit none
! arguments
integer :: xctype(3)
character(512) :: xcdescr
integer :: xcspin
integer :: xcgrad
logical :: hybrid
real(8) :: hybridc
write(*,*)
write(*,'("Error(libxcifc):  libxc not or improperly installed")')
write(*,*)
stop
end subroutine
!EOC

end module

