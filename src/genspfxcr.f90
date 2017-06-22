
! Copyright (C) 2012 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genspfxcr(tsh,fxcmt,fxcir)
use modmain
use modtddft
use modfxcifc
implicit none
! arguments
logical, intent(in) :: tsh
real(8), intent(out) :: fxcmt(lmmaxvr,nrmtmax,natmtot,4,4)
real(8), intent(out) :: fxcir(ngtot,4,4)
! local variables
integer idm,is,ia,ias
integer nr,nri,ir,ld,i,j,n
real(8) t1
! allocatable arrays
real(8), allocatable :: rho(:),rhoup(:),rhodn(:)
real(8), allocatable :: mag(:,:),magu(:,:),magm(:)
real(8), allocatable :: bxc(:,:),bxcp(:)
real(8), allocatable :: fxcuu(:),fxcud(:),fxcdd(:)
real(8), allocatable :: fxc(:,:,:)
if (.not.spinpol) then
  write(*,*)
  write(*,'("Error(genspfxcr): spin-unpolarised calculation")')
  write(*,*)
  stop
end if
! allocate local arrays
n=lmmaxvr*nrmtmax
allocate(rho(n),mag(n,ndmag))
allocate(bxc(n,ndmag),fxc(n,4,4))
n=max(n,ngtot)
allocate(rhoup(n),rhodn(n))
allocate(magu(3,n),magm(n),bxcp(n))
allocate(fxcuu(n),fxcud(n),fxcdd(n))
!---------------------------!
!     muffin-tin kernel     !
!---------------------------!
ld=lmmaxvr*nrmtmax
do is=1,nspecies
  nr=nrmt(is)
  nri=nrmtinr(is)
  n=lmmaxvr*nr
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! compute the density in spherical coordinates
    call rbsht(nr,nri,1,rhomt(:,:,ias),1,rho)
    do idm=1,ndmag
! magnetisation in spherical coordinates
      call rbsht(nr,nri,1,magmt(:,:,ias,idm),1,mag(:,idm))
! B_xc in spherical coordinates
      call rbsht(nr,nri,1,bxcmt(:,:,ias,idm),1,bxc(:,idm))
    end do
    if (ncmag) then
! non-collinear (use Kubler's trick)
      do i=1,n
! compute |m|
        magm(i)=sqrt(mag(i,1)**2+mag(i,2)**2+mag(i,3)**2)
! compute rhoup=(rho+|m|)/2 and rhodn=(rho-|m|)/2
        rhoup(i)=0.5d0*(rho(i)+magm(i))
        rhodn(i)=0.5d0*(rho(i)-magm(i))
! unit vector m/|m|
        t1=1.d0/(magm(i)+1.d-8)
        magu(1,i)=t1*mag(i,1)
        magu(2,i)=t1*mag(i,2)
        magu(3,i)=t1*mag(i,3)
! compute B_xc.(m/|m|)
        bxcp(i)=bxc(i,1)*magu(1,i)+bxc(i,2)*magu(2,i)+bxc(i,3)*magu(3,i)
      end do
    else
! collinear
      do i=1,n
! compute |m| = |m_z|
        magm(i)=abs(mag(i,1))
! compute rhoup=(rho+|m|)/2 and rhodn=(rho-|m|)/2
        rhoup(i)=0.5d0*(rho(i)+magm(i))
        rhodn(i)=0.5d0*(rho(i)-magm(i))
! unit vector m/|m|
        magu(1,i)=0.d0
        magu(2,i)=0.d0
        if (mag(i,1).gt.0.d0) then
          magu(3,i)=1.d0
        else
          magu(3,i)=-1.d0
        end if
! compute B_xc.(m/|m|)
        bxcp(i)=bxc(i,1)*magu(3,i)
      end do
    end if
! compute f_xc in U(2) x U(2) basis
    call fxcifc(fxctype,n=n,rhoup=rhoup,rhodn=rhodn,fxcuu=fxcuu,fxcud=fxcud, &
     fxcdd=fxcdd)
! transform f_xc to O(1) x O(3) basis (upper triangular part)
    call tfm2213(n,fxcuu,fxcud,fxcdd,magu,magm,bxcp,ld,fxc)
    do i=1,4
      do j=i,4
        if (tsh) then
! convert to spherical harmonics if required
          call rfsht(nr,nri,1,fxc(:,i,j),1,fxcmt(:,:,ias,i,j))
        else
          call dcopy(n,fxc(:,i,j),1,fxcmt(:,:,ias,i,j),1)
        end if
      end do
    end do
  end do
end do
!-----------------------------!
!     interstitial kernel     !
!-----------------------------!
if (ncmag) then
! non-collinear
  do ir=1,ngtot
    magm(ir)=sqrt(magir(ir,1)**2+magir(ir,2)**2+magir(ir,3)**2)
    rhoup(ir)=0.5d0*(rhoir(ir)+magm(ir))
    rhodn(ir)=0.5d0*(rhoir(ir)-magm(ir))
    t1=1.d0/(magm(ir)+1.d-8)
    magu(1,ir)=t1*magir(ir,1)
    magu(2,ir)=t1*magir(ir,2)
    magu(3,ir)=t1*magir(ir,3)
! compute B_xc.(m/|m|)
    bxcp(ir)=bxcir(ir,1)*magu(1,ir) &
            +bxcir(ir,2)*magu(2,ir) &
            +bxcir(ir,3)*magu(3,ir)
  end do
else
! collinear
  do ir=1,ngtot
    magm(ir)=abs(magir(ir,1))
    rhoup(ir)=0.5d0*(rhoir(ir)+magm(ir))
    rhodn(ir)=0.5d0*(rhoir(ir)-magm(ir))
    magu(1,ir)=0.d0
    magu(2,ir)=0.d0
    if (magir(ir,1).gt.0.d0) then
      magu(3,ir)=1.d0
    else
      magu(3,ir)=-1.d0
    end if
! compute B_xc.(m/|m|)
    bxcp(ir)=bxcir(ir,1)*magu(3,ir)
  end do
end if
! compute f_xc in U(2) x U(2) basis
call fxcifc(fxctype,n=ngtot,rhoup=rhoup,rhodn=rhodn,fxcuu=fxcuu,fxcud=fxcud, &
 fxcdd=fxcdd)
! transform f_xc to O(1) x O(3) basis
call tfm2213(ngtot,fxcuu,fxcud,fxcdd,magu,magm,bxcp,ngtot,fxcir)
deallocate(rho,mag,bxc,fxc)
deallocate(rhoup,rhodn)
deallocate(magu,magm,bxcp)
deallocate(fxcuu,fxcud,fxcdd)
return

contains

subroutine tfm2213(n,fxcuu,fxcud,fxcdd,magu,magm,bxcp,ld,fxc)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: fxcuu(n),fxcud(n),fxcdd(n)
real(8), intent(in) :: magu(3,n),magm(n),bxcp(n)
integer, intent(in) :: ld
real(8), intent(out) :: fxc(ld,4,4)
! local variables
integer i
real(8) t1,t2
do i=1,n
! charge-charge
  fxc(i,1,1)=0.25d0*(fxcuu(i)+2.d0*fxcud(i)+fxcdd(i))
! charge-spin
  t1=0.25d0*(fxcuu(i)-fxcdd(i))
  fxc(i,1,2)=t1*magu(1,i)
  fxc(i,1,3)=t1*magu(2,i)
  fxc(i,1,4)=t1*magu(3,i)
! spin-spin
  if (magm(i).gt.1.d-14) then
    t1=bxcp(i)/magm(i)
  else
    t1=0.d0
  end if
  t2=0.25d0*(fxcuu(i)-2.d0*fxcud(i)+fxcdd(i))-t1
  fxc(i,2,2)=t2*magu(1,i)*magu(1,i)+t1
  fxc(i,2,3)=t2*magu(1,i)*magu(2,i)
  fxc(i,2,4)=t2*magu(1,i)*magu(3,i)
  fxc(i,3,3)=t2*magu(2,i)*magu(2,i)+t1
  fxc(i,3,4)=t2*magu(2,i)*magu(3,i)
  fxc(i,4,4)=t2*magu(3,i)*magu(3,i)+t1
end do
return
end subroutine

end subroutine

