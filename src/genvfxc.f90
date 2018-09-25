
! Copyright (C) 2011 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvfxc(tq0,t3hw,gclgq,nm,vchi0,eps0,epsi,vfxc)
use modmain
use modtddft
implicit none
! arguments
logical, intent(in) :: tq0,t3hw
real(8), intent(in) :: gclgq(ngrf)
integer, intent(in) :: nm
complex(8), intent(in) :: vchi0(nm,nm,nwrf)
complex(8), intent(in) :: eps0(nm,nm,nwrf)
complex(8), intent(in) :: epsi(nm,nm,nwrf)
complex(8), intent(out) :: vfxc(nm,nm,nwrf)
! local variables
integer iw,i,j
complex(8) z1
! allocatable arrays
complex(8), allocatable :: a(:,:)
! compute v^(-1/2) f_xc v^(-1/2)
select case(fxctype(1))
case(0,1)
! RPA
  vfxc(:,:,:)=0.d0
  return
case(3)
! ALDA
  if (tq0.and.t3hw) then
    call genvfxcg(gclgq,nm,vfxc(3,3,1))
! the head and wings are zero
    vfxc(1:3,:,:)=0.d0
    vfxc(4:,1:3,:)=0.d0
  else
    call genvfxcg(gclgq,nm,vfxc)
  end if
case(200)
! long-range contribution with dynamic correlations
  vfxc(:,:,:)=0.d0
  do i=1,nm
    vfxc(i,i,:)=-(fxclrc(1)+fxclrc(2)*dble(wrf(:))**2)/fourpi
  end do
case(210,211)
! bootstrap
  vfxc(:,:,:)=0.d0
  if (tq0.and.t3hw) then
    z1=(eps0(1,1,1)+eps0(2,2,1)+eps0(3,3,1))/3.d0
  else
    z1=eps0(1,1,1)
  end if
  z1=-1.d0/(z1-1.d0)
  do i=1,nm
    do j=1,nm
      vfxc(i,j,:)=z1*epsi(i,j,1)
    end do
  end do
case default
  write(*,*)
  write(*,'("Error(genvfxc): fxctype not defined : ",I8)') fxctype
  write(*,*)
  stop
end select
! right multiply by v^1/2 chi0 v^1/2
allocate(a(nm,nm))
do iw=1,nwrf
  a(:,:)=vfxc(:,:,iw)
  call zgemm('N','N',nm,nm,nm,zone,a,nm,vchi0(:,:,iw),nm,zzero,vfxc(:,:,iw),nm)
end do
deallocate(a)
return
end subroutine


