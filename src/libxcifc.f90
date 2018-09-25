
! Copyright (C) 2009 T. McQueen and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module libxcifc

use xc_f90_lib_m

! libxc version number
integer libxcv(3)

contains

!BOP
! !ROUTINE: xcifc_libxc
! !INTERFACE:
subroutine xcifc_libxc(xctype,n,c_tb09,rho,rhoup,rhodn,g2rho,g2up,g2dn,grho2, &
 gup2,gdn2,gupdn,tau,tauup,taudn,ex,ec,vx,vc,vxup,vxdn,vcup,vcdn,dxdgr2, &
 dxdgu2,dxdgd2,dxdgud,dcdgr2,dcdgu2,dcdgd2,dcdgud,dxdg2r,dxdg2u,dxdg2d,dcdg2r, &
 dcdg2u,dcdg2d,wx,wxup,wxdn,wc,wcup,wcdn)
! !INPUT/OUTPUT PARAMETERS:
!   xctype : type of exchange-correlation functional (in,integer(3))
!   n      : number of density points (in,integer)
!   c_tb09 : Tran-Blaha '09 constant c (in,real,optional)
!   rho    : spin-unpolarised charge density (in,real(n),optional)
!   rhoup  : spin-up charge density (in,real(n),optional)
!   rhodn  : spin-down charge density (in,real(n),optional)
!   g2rho  : grad^2 rho (in,real(n),optional)
!   g2up   : grad^2 rhoup (in,real(n),optional)
!   g2dn   : grad^2 rhodn (in,real(n),optional)
!   grho2  : |grad rho|^2 (in,real(n),optional)
!   gup2   : |grad rhoup|^2 (in,real(n),optional)
!   gdn2   : |grad rhodn|^2 (in,real(n),optional)
!   gupdn  : (grad rhoup).(grad rhodn) (in,real(n),optional)
!   tau    : kinetic energy density (in,real(n),optional)
!   tauup  : spin-up kinetic energy density (in,real(n),optional)
!   taudn  : spin-down kinetic energy density (in,real(n),optional)
!   ex     : exchange energy density (out,real(n),optional)
!   ec     : correlation energy density (out,real(n),optional)
!   vx     : spin-unpolarised exchange potential (out,real(n),optional)
!   vc     : spin-unpolarised correlation potential (out,real(n),optional)
!   vxup   : spin-up exchange potential (out,real(n),optional)
!   vxdn   : spin-down exchange potential (out,real(n),optional)
!   vcup   : spin-up correlation potential (out,real(n),optional)
!   vcdn   : spin-down correlation potential (out,real(n),optional)
!   dxdgr2 : de_x/d(|grad rho|^2) (out,real(n),optional)
!   dxdgu2 : de_x/d(|grad rhoup|^2) (out,real(n),optional)
!   dxdgd2 : de_x/d(|grad rhodn|^2) (out,real(n),optional)
!   dxdgud : de_x/d((grad rhoup).(grad rhodn)) (out,real(n),optional)
!   dcdgr2 : de_c/d(|grad rho|^2) (out,real(n),optional)
!   dcdgu2 : de_c/d(|grad rhoup|^2) (out,real(n),optional)
!   dcdgd2 : de_c/d(|grad rhodn|^2) (out,real(n),optional)
!   dcdgud : de_c/d((grad rhoup).(grad rhodn)) (out,real(n),optional)
!   dxdg2r : de_x/d(grad^2 rho) (out,real(n),optional)
!   dxdg2u : de_x/d(grad^2 rhoup) (out,real(n),optional)
!   dxdg2d : de_x/d(grad^2 rhodn) (out,real(n),optional)
!   dcdg2r : de_c/d(grad^2 rho) (out,real(n),optional)
!   dcdg2u : de_c/d(grad^2 rhoup) (out,real(n),optional)
!   dcdg2d : de_c/d(grad^2 rhodn) (out,real(n),optional)
!   wx     : de_x/dtau (out,real(n),optional)
!   wxup   : de_x/dtauup (out,real(n),optional)
!   wxdn   : de_x/dtaudn (out,real(n),optional)
!   wc     : de_c/dtau (out,real(n),optional)
!   wcup   : de_c/dtauup (out,real(n),optional)
!   wcdn   : de_c/dtaudn (out,real(n),optional)
! !DESCRIPTION:
!   Interface to the ETSF {\tt libxc} exchange-correlation functional library:
!   \newline{\tt http://www.tddft.org/programs/octopus/wiki/index.php/Libxc}.
!   The second and third integers in {\tt xctype} define the exchange and
!   correlation functionals in {\tt libxc}, respectively.
!
! !REVISION HISTORY:
!   Created April 2009 (Tyrel McQueen)
!   Modified September 2009 (JKD and TMQ)
!   Updated for Libxc 1, July 2010 (JKD)
!   Updated for Libxc 4, March 2018 (JKD)
!EOP
!BOC
implicit none
! mandatory arguments
integer, intent(in) :: xctype(3),n
! optional arguments
real(8), optional, intent(in) :: c_tb09
real(8), optional, intent(in) :: rho(n),rhoup(n),rhodn(n)
real(8), optional, intent(in) :: g2rho(n),g2up(n),g2dn(n)
real(8), optional, intent(in) :: grho2(n),gup2(n),gdn2(n),gupdn(n)
real(8), optional, intent(in) :: tau(n),tauup(n),taudn(n)
real(8), optional, intent(out) :: ex(n),ec(n),vx(n),vc(n)
real(8), optional, intent(out) :: vxup(n),vxdn(n),vcup(n),vcdn(n)
real(8), optional, intent(out) :: dxdgr2(n),dxdgu2(n),dxdgd2(n),dxdgud(n)
real(8), optional, intent(out) :: dxdg2r(n),dxdg2u(n),dxdg2d(n)
real(8), optional, intent(out) :: wx(n),wxup(n),wxdn(n)
real(8), optional, intent(out) :: dcdgr2(n),dcdgu2(n),dcdgd2(n),dcdgud(n)
real(8), optional, intent(out) :: dcdg2r(n),dcdg2u(n),dcdg2d(n)
real(8), optional, intent(out) :: wc(n),wcup(n),wcdn(n)
! local variables
integer nspin,xcf,id,k
type(xc_f90_pointer_t) p,info
! allocatable arrays
real(8), allocatable :: r(:,:),sigma(:,:),vrho(:,:),vsigma(:,:)
real(8), allocatable :: lapl(:,:),t(:,:),vlapl(:,:),vtau(:,:)
if (present(rho)) then
  nspin=XC_UNPOLARIZED
else if (present(rhoup).and.present(rhodn)) then
  nspin=XC_POLARIZED
else
  write(*,*)
  write(*,'("Error(xcifc_libxc): missing arguments")')
  write(*,*)
  stop
end if
if (xctype(2).ne.0) then
  if (xctype(2).eq.xctype(3)) then
    write(*,*)
    write(*,'("Error(xcifc_libxc): Libxc exchange and correlation &
     &functionals")')
    write(*,'(" are the same : ",2I8)') xctype(2:3)
    write(*,*)
    stop
  end if
end if
! loop over functional kinds (exchange or correlation)
do k=2,3
  id=xctype(k)
  if (id.gt.0) then
    xcf=xc_f90_family_from_id(id)
! initialise functional
    call xc_f90_func_init(p,info,id,nspin)
    select case(xcf)
    case(XC_FAMILY_LDA)
!-------------------------!
!     LDA functionals     !
!-------------------------!
      if (k.eq.2) then
! exchange
        if (present(rho)) then
          call xc_f90_lda_exc_vxc(p,n,rho(1),ex(1),vx(1))
        else
          allocate(r(2,n),vrho(2,n))
          r(1,:)=rhoup(:); r(2,:)=rhodn(:)
          call xc_f90_lda_exc_vxc(p,n,r(1,1),ex(1),vrho(1,1))
          vxup(:)=vrho(1,:); vxdn(:)=vrho(2,:)
          deallocate(r,vrho)
        end if
      else
! correlation
        if (present(rho)) then
          call xc_f90_lda_exc_vxc(p,n,rho(1),ec(1),vc(1))
        else
          allocate(r(2,n),vrho(2,n))
          r(1,:)=rhoup(:); r(2,:)=rhodn(:)
          call xc_f90_lda_exc_vxc(p,n,r(1,1),ec(1),vrho(1,1))
          vcup(:)=vrho(1,:); vcdn=vrho(2,:)
          deallocate(r,vrho)
        end if
      end if
    case(XC_FAMILY_GGA,XC_FAMILY_HYB_GGA)
!-------------------------!
!     GGA functionals     !
!-------------------------!
      if (k.eq.2) then
! exchange
        if (present(rho)) then
          call xc_f90_gga_exc_vxc(p,n,rho(1),grho2(1),ex(1),vx(1),dxdgr2(1))
        else
          allocate(r(2,n),sigma(3,n),vrho(2,n),vsigma(3,n))
          r(1,:)=rhoup(:); r(2,:)=rhodn(:)
          sigma(1,:)=gup2(:); sigma(2,:)=gupdn(:); sigma(3,:)=gdn2(:)
          call xc_f90_gga_exc_vxc(p,n,r(1,1),sigma(1,1),ex(1),vrho(1,1), &
           vsigma(1,1))
          vxup(:)=vrho(1,:); vxdn(:)=vrho(2,:)
          dxdgu2(:)=vsigma(1,:); dxdgud(:)=vsigma(2,:); dxdgd2(:)=vsigma(3,:)
          deallocate(r,sigma,vrho,vsigma)
        end if
      else
! correlation
        if (present(rho)) then
          call xc_f90_gga_exc_vxc(p,n,rho(1),grho2(1),ec(1),vc(1),dcdgr2(1))
        else
          allocate(r(2,n),sigma(3,n),vrho(2,n),vsigma(3,n))
          r(1,:)=rhoup(:); r(2,:)=rhodn(:)
          sigma(1,:)=gup2(:); sigma(2,:)=gupdn(:); sigma(3,:)=gdn2(:)
          call xc_f90_gga_exc_vxc(p,n,r(1,1),sigma(1,1),ec(1),vrho(1,1), &
           vsigma(1,1))
          vcup(:)=vrho(1,:); vcdn(:)=vrho(2,:)
          dcdgu2(:)=vsigma(1,:); dcdgud(:)=vsigma(2,:); dcdgd2(:)=vsigma(3,:)
          deallocate(r,sigma,vrho,vsigma)
        end if
      end if
    case(XC_FAMILY_MGGA)
!------------------------------!
!     meta-GGA functionals     !
!------------------------------!
! set Tran-Blaha '09 constant if required
      if (id.eq.XC_MGGA_X_TB09) then
        if (present(c_tb09)) call xc_f90_func_set_ext_params(p,c_tb09)
      end if
      if (k.eq.2) then
! exchange
        if (present(rho)) then
          if (present(ex)) then
! spin-unpolarised energy functional
            call xc_f90_mgga_exc_vxc(p,n,rho(1),grho2(1),g2rho(1),tau(1), &
             ex(1),vx(1),dxdgr2(1),dxdg2r(1),wx(1))
          else
! spin-unpolarised potential-only functional
            allocate(vsigma(1,n),vlapl(1,n),vtau(1,n))
            call xc_f90_mgga_vxc(p,n,rho(1),grho2(1),g2rho(1),tau(1),vx(1), &
             vsigma(1,1),vlapl(1,1),vtau(1,1))
            deallocate(vsigma,vlapl,vtau)
          end if
        else
          allocate(r(2,n),sigma(3,n),lapl(2,n),t(2,n))
          allocate(vrho(2,n),vsigma(3,n),vlapl(2,n),vtau(2,n))
          r(1,:)=rhoup(:); r(2,:)=rhodn(:)
          sigma(1,:)=gup2(:); sigma(2,:)=gupdn(:); sigma(3,:)=gdn2(:)
          lapl(1,:)=g2up(:); lapl(2,:)=g2dn(:)
          t(1,:)=tauup(:); t(2,:)=taudn(:)
          if (present(ex)) then
! spin-polarised energy functional
            call xc_f90_mgga_exc_vxc(p,n,r(1,1),sigma(1,1),lapl(1,1),t(1,1), &
             ex(1),vrho(1,1),vsigma(1,1),vlapl(1,1),vtau(1,1))
            dxdgu2(:)=vsigma(1,:); dxdgud(:)=vsigma(2,:); dxdgd2(:)=vsigma(3,:)
            dxdg2u(:)=vlapl(1,:); dxdg2d(:)=vlapl(2,:)
            wxup(:)=vtau(1,:); wxdn(:)=vtau(2,:)
          else
! spin-polarised potential-only functional
            call xc_f90_mgga_vxc(p,n,r(1,1),sigma(1,1),lapl(1,1),t(1,1), &
             vrho(1,1),vsigma(1,1),vlapl(1,1),vtau(1,1))
          end if
          vxup(:)=vrho(1,:); vxdn(:)=vrho(2,:)
          deallocate(r,sigma,lapl,t)
          deallocate(vrho,vsigma,vlapl,vtau)
        end if
      else
! correlation
        if (present(rho)) then
          if (present(ec)) then
            call xc_f90_mgga_exc_vxc(p,n,rho(1),grho2(1),g2rho(1),tau(1), &
             ec(1),vc(1),dcdgr2(1),dcdg2r(1),wc(1))
          else
            allocate(vsigma(1,n),vlapl(1,n),vtau(1,n))
            call xc_f90_mgga_vxc(p,n,rho(1),grho2(1),g2rho(1),tau(1),vc(1), &
             vsigma(1,1),vlapl(1,1),vtau(1,1))
            deallocate(vsigma,vlapl,vtau)
          end if
        else
          allocate(r(2,n),sigma(3,n),lapl(2,n),t(2,n))
          allocate(vrho(2,n),vsigma(3,n),vlapl(2,n),vtau(2,n))
          r(1,:)=rhoup(:); r(2,:)=rhodn(:)
          sigma(1,:)=gup2(:); sigma(2,:)=gupdn(:); sigma(3,:)=gdn2(:)
          lapl(1,:)=g2up(:); lapl(2,:)=g2dn(:)
          t(1,:)=tauup(:); t(2,:)=taudn(:)
          if (present(ec)) then
            call xc_f90_mgga_exc_vxc(p,n,r(1,1),sigma(1,1),lapl(1,1),t(1,1), &
             ec(1),vrho(1,1),vsigma(1,1),vlapl(1,1),vtau(1,1))
            dcdgu2(:)=vsigma(1,:); dcdgud(:)=vsigma(2,:); dcdgd2(:)=vsigma(3,:)
            dcdg2u(:)=vlapl(1,:); dcdg2d(:)=vlapl(2,:)
            wcup(:)=vtau(1,:); wcdn(:)=vtau(2,:)
          else
            call xc_f90_mgga_vxc(p,n,r(1,1),sigma(1,1),lapl(1,1),t(1,1), &
             vrho(1,1),vsigma(1,1),vlapl(1,1),vtau(1,1))
          end if
          vcup(:)=vrho(1,:); vcdn(:)=vrho(2,:)
          deallocate(r,sigma,lapl,t)
          deallocate(vrho,vsigma,vlapl,vtau)
        end if
      end if
    case default
      write(*,*)
      write(*,'("Error(xcifc_libxc): unsupported Libxc functional family : ",&
       &I8)') xcf
      write(*,*)
      stop
    end select
! destroy functional
    call xc_f90_func_end(p)
  else
! case when id=0
    if (k.eq.2) then
      if (present(ex)) ex(:)=0.d0
      if (present(vx)) vx(:)=0.d0
      if (present(vxup)) vxup(:)=0.d0
      if (present(vxdn)) vxdn(:)=0.d0
      if (present(dxdgr2)) dxdgr2(:)=0.d0
      if (present(dxdgu2)) dxdgu2(:)=0.d0
      if (present(dxdgd2)) dxdgd2(:)=0.d0
      if (present(dxdgud)) dxdgud(:)=0.d0
    else
      if (present(ec)) ec(:)=0.d0
      if (present(vc)) vc(:)=0.d0
      if (present(vcup)) vcup(:)=0.d0
      if (present(vcdn)) vcdn(:)=0.d0
      if (present(dcdgr2)) dcdgr2(:)=0.d0
      if (present(dcdgu2)) dcdgu2(:)=0.d0
      if (present(dcdgd2)) dcdgd2(:)=0.d0
      if (present(dcdgud)) dcdgud(:)=0.d0
    end if
  end if
end do
return
end subroutine

subroutine fxcifc_libxc(fxctype,n,rho,rhoup,rhodn,fxc,fxcuu,fxcud,fxcdd)
implicit none
! mandatory arguments
integer, intent(in) :: fxctype(3),n
! optional arguments
real(8), optional, intent(in) :: rho(n),rhoup(n),rhodn(n)
real(8), optional, intent(out) :: fxc(n),fxcuu(n),fxcud(n),fxcdd(n)
! local variables
integer nspin,xcf,id,k
type(xc_f90_pointer_t) p,info
! allocatable arrays
real(8), allocatable :: r(:,:),f(:,:)
if (present(rho)) then
  nspin=XC_UNPOLARIZED
else if (present(rhoup).and.present(rhodn)) then
  nspin=XC_POLARIZED
else
  write(*,*)
  write(*,'("Error(fxcifc_libxc): missing arguments")')
  write(*,*)
  stop
end if
! zero the kernel
if (present(fxc)) fxc(:)=0.d0
if (present(fxcuu)) fxcuu(:)=0.d0
if (present(fxcud)) fxcud(:)=0.d0
if (present(fxcdd)) fxcdd(:)=0.d0
! loop over functional kinds (exchange or correlation)
do k=2,3
  id=fxctype(k)
  if (id.le.0) cycle
  xcf=xc_f90_family_from_id(id)
! initialise functional
  call xc_f90_func_init(p,info,id,nspin)
  select case(xcf)
  case(XC_FAMILY_LDA)
!-------------------------!
!     LDA functionals     !
!-------------------------!
    if (present(rho)) then
      allocate(f(1,n))
      call xc_f90_lda_fxc(p,n,rho(1),f(1,1))
      fxc(:)=fxc(:)+f(1,:)
      deallocate(f)
    else
      allocate(r(2,n),f(3,n))
      r(1,:)=rhoup(:); r(2,:)=rhodn(:)
      call xc_f90_lda_fxc(p,n,r(1,1),f(1,1))
      fxcuu(:)=fxcuu(:)+f(1,:)
      fxcud(:)=fxcud(:)+f(2,:)
      fxcdd(:)=fxcdd(:)+f(3,:)
      deallocate(r,f)
    end if
  case default
    write(*,*)
    write(*,'("Error(fxcifc_libxc): unsupported Libxc functional family : ",&
     &I8)') xcf
    write(*,*)
    stop
  end select
end do
return
end subroutine

subroutine xcdata_libxc(xctype,xcdescr,xcspin,xcgrad,hybrid,hybridc)
implicit none
! arguments
integer, intent(in) :: xctype(3)
character(512), intent(out) :: xcdescr
integer, intent(out) :: xcspin
integer, intent(out) :: xcgrad
logical, intent(inout) :: hybrid
real(8), intent(inout) :: hybridc
! local variables
integer j,k,id
integer fmly,flg
character(256) name
type(xc_f90_pointer_t) p,info
! check version is compatible
call xc_f90_version(libxcv(1),libxcv(2),libxcv(3))
if (libxcv(1).ne.4) then
  write(*,*)
  write(*,'("Error(xcdata_libxc): incompatible Libxc version : ",I2.2,".",&
   &I3.3,".",I3.3)') libxcv(:)
  write(*,*)
  stop
end if
! unknown spin polarisation
xcspin=-1
! no gradients by default
xcgrad=0
! not hybrid by default
hybrid=.false.
do k=2,3
  id=xctype(k)
  if (id.gt.0) then
    call xc_f90_func_init(p,info,id,XC_UNPOLARIZED)
    call xc_f90_info_name(info,name)
    fmly=xc_f90_family_from_id(id)
    if (fmly.eq.XC_FAMILY_HYB_GGA) then
      call xc_f90_hyb_exx_coef(p,hybridc)
      hybrid=.true.
    end if
    call xc_f90_func_end(p)
! hybrids should have only correlation part to avoid double-scaling with hybridc
    if ((fmly.eq.XC_FAMILY_HYB_GGA).and.(k.eq.2)) then
      write(*,*)
      write(*,'("Error(xcdata_libxc): set only correlation part of xctype for &
       &Libxc hybrids")')
      write(*,*)
      stop
    end if
! functional family
    if ((fmly.ne.XC_FAMILY_LDA).and.(fmly.ne.XC_FAMILY_GGA).and. &
     (fmly.ne.XC_FAMILY_HYB_GGA).and.(fmly.ne.XC_FAMILY_MGGA)) then
      write(*,*)
      write(*,'("Error(xcdata_libxc): unsupported Libxc family : ",I8)') fmly
      write(*,*)
      stop
    end if
! post-processed gradients required for GGA functionals
    if (fmly.eq.XC_FAMILY_GGA.or.fmly.eq.XC_FAMILY_HYB_GGA) xcgrad=2
! kinetic energy density required for meta-GGA functionals
    if (fmly.eq.XC_FAMILY_MGGA) then
      flg=xc_f90_info_flags(info)
      if ((iand(flg,XC_FLAGS_HAVE_VXC).ne.0).and. &
       (iand(flg,XC_FLAGS_HAVE_EXC).eq.0)) then
! potential-only functional
        xcgrad=3
      else if (iand(flg,XC_FLAGS_HAVE_EXC).ne.0) then
! energy functional
        xcgrad=4
      else
        write(*,*)
        write(*,'("Error(xcdata_libxc): unsupported Libxc meta-GGA type")')
        write(*,*)
        stop
      end if
! check if the density laplacian is required
      if ((xcgrad.eq.4).and.(iand(flg,XC_FLAGS_NEEDS_LAPLACIAN).ne.0)) then
        write(*,*)
        write(*,'("Error(xcdata_libxc): energy meta-GGAs requiring the density &
         &laplacian are not supported")')
        write(*,*)
        stop
      end if
    end if
  else
    name='none'
  end if
  if (k.eq.2) then
    xcdescr='exchange: '//trim(name)
  else
    xcdescr=trim(xcdescr)//'; correlation: '//trim(name)
  end if
end do
return
end subroutine
!EOC

end module

