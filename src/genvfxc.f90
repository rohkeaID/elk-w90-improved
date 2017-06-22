
! Copyright (C) 2011 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvfxc(gqc,vchi0,eps0,eps,vfxc)
use modmain
use modtddft
implicit none
! arguments
real(8), intent(in) :: gqc(ngrf)
complex(8), intent(in) :: vchi0(nwrf,ngrf,ngrf)
complex(8), intent(in) :: eps0(ngrf,ngrf,nwrf)
complex(8), intent(in) :: eps(ngrf,ngrf,nwrf)
complex(8), intent(out) :: vfxc(ngrf,ngrf,nwrf)
! local variables
integer ig,jg,iw,info
real(8) t1
complex(8) z1
! allocatable arrays
integer, allocatable :: ipiv(:)
complex(8), allocatable :: a(:,:),b(:,:),work(:)
! compute v^(-1/2) f_xc v^(-1/2)
select case(fxctype(1))
case(0,1)
! RPA
  vfxc(:,:,:)=0.d0
  return
case(3)
! ALDA
  call genvfxcg(gqc,vfxc)
case(200)
! long-range contribution with dynamic correlations
  vfxc(:,:,:)=0.d0
  do ig=1,ngrf
    vfxc(ig,ig,:)=-(fxclrc(1)+fxclrc(2)*dble(wrf(:))**2)/fourpi
  end do
case(210)
! bootstrap
  vfxc(:,:,:)=0.d0
  t1=-1.d0/(dble(eps0(1,1,1))-1.d0)
  do ig=1,ngrf
    do jg=1,ngrf
      vfxc(ig,jg,:)=t1*eps(ig,jg,1)
    end do
  end do
case(211)
! single iteration bootstrap
  vfxc(:,:,:)=0.d0
  allocate(ipiv(ngrf),a(ngrf,ngrf),work(ngrf))
  a(:,:)=eps0(:,:,1)
! invert RPA epsilon
  call zgetrf(ngrf,ngrf,a,ngrf,ipiv,info)
  if (info.eq.0) call zgetri(ngrf,a,ngrf,ipiv,work,ngrf,info)
  if (info.ne.0) then
    write(*,*)
    write(*,'("Error(genvfxc): unable to invert RPA epsilon")')
    write(*,*)
    stop
  end if
  do ig=1,ngrf
    do jg=1,ngrf
      vfxc(ig,jg,:)=a(ig,jg)/vchi0(1,1,1)
    end do
  end do
  deallocate(ipiv,a,work)
case(212)
! Revised bootstrap (RBO)
! See: S. Rigamonti, et al., Phys. Rev. Lett. 114, 146402
  vfxc(:,:,:)=0.d0
  allocate(ipiv(ngrf),a(ngrf,ngrf),b(ngrf,ngrf),work(ngrf))
  a(:,:)=eps0(:,:,1)
! invert RPA epsilon
  call zgetrf(ngrf,ngrf,a,ngrf,ipiv,info)
  if (info.eq.0) call zgetri(ngrf,a,ngrf,ipiv,work,ngrf,info)
  if (info.ne.0) then
    write(*,*)
    write(*,'("Error(genvfxc): unable to invert RPA epsilon")')
    write(*,*)
    stop
  end if
! calculate RPA chi
  b(:,:)=a(:,:)
  do ig=1,ngrf
    b(ig,ig)=b(ig,ig)-1.d0
  end do
! invert RPA chi
  call zgetrf(ngrf,ngrf,a,ngrf,ipiv,info)
  if (info.eq.0) call zgetri(ngrf,b,ngrf,ipiv,work,ngrf,info)
  if (info.ne.0) then
    write(*,*)
    write(*,'("Error(genvfxc): unable to invert RPA chi")')
    write(*,*)
    stop
  end if
  z1=a(1,1)
  do ig=1,ngrf
    do jg=1,ngrf
      vfxc(ig,jg,:)=z1*b(ig,jg)
    end do
  end do
  vfxc(1,1,:)=a(1,1)/vchi0(1,1,1)
  deallocate(ipiv,a,b,work)
case default
  write(*,*)
  write(*,'("Error(genvfxc): fxctype not defined : ",I8)') fxctype
  write(*,*)
  stop
end select
! right multiply by v^1/2 chi0 v^1/2
allocate(a(ngrf,ngrf),b(ngrf,ngrf))
do iw=1,nwrf
  a(:,:)=vfxc(:,:,iw)
  b(:,:)=vchi0(iw,:,:)
  call zgemm('N','N',ngrf,ngrf,ngrf,zone,a,ngrf,b,ngrf,zzero,vfxc(:,:,iw),ngrf)
end do
deallocate(a,b)
return
end subroutine

