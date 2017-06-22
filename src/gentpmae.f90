
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gentpmae
use modmain
implicit none
! local variables
integer na,n,i1,i2,i3,i
integer isym,lspl
real(8) v1(3),v2(3),t1
! allocatable arrays
real(8), allocatable :: vp(:,:)
if (allocated(tpmae)) deallocate(tpmae)
select case(npmae0)
case(4:)
! distribute points evenly on a sphere
  npmae=npmae0
  allocate(tpmae(2,npmae))
  call sphcover(npmae,tpmae)
case(-4:-1)
! use symmetry reduced cardinal directions
  na=abs(npmae0)
  n=(2*na+1)**3
  allocate(vp(3,n))
  npmae=0
  do i1=-na,na
    v1(1)=dble(i1)
    do i2=-na,na
      v1(2)=dble(i2)
      do i3=-na,na
        v1(3)=dble(i3)
        if ((i1.eq.0).and.(i2.eq.0).and.(i3.eq.0)) cycle
        do isym=1,nsymcrys
          lspl=lsplsymc(isym)
          v2(:)=symlat(:,1,lspl)*v1(1) &
               +symlat(:,2,lspl)*v1(2) &
               +symlat(:,3,lspl)*v1(3)
          do i=1,npmae
            t1=abs(vp(1,i)-v2(1))+abs(vp(2,i)-v2(2))+abs(vp(3,i)-v2(3))
            if (t1.lt.epslat) goto 10
          end do
        end do
        npmae=npmae+1
        vp(:,npmae)=v1(:)
10 continue
      end do
    end do
  end do
! convert vectors to spherical coordinates
  allocate(tpmae(2,npmae))
  do i=1,npmae
    call sphcrd(vp(:,i),t1,tpmae(:,i))
  end do
  deallocate(vp)
case(2)
! use x- and z-directions
  npmae=2
  allocate(tpmae(2,npmae))
  tpmae(1,1)=pi/2.d0
  tpmae(2,1)=0.d0
  tpmae(1,2)=0.d0
  tpmae(2,2)=0.d0
case(3)
! use x-, y- and z-directions
  npmae=3
  allocate(tpmae(2,npmae))
  tpmae(1,1)=pi/2.d0
  tpmae(2,1)=0.d0
  tpmae(1,2)=pi/2.d0
  tpmae(2,2)=pi/2.d0
  tpmae(1,3)=0.d0
  tpmae(2,3)=0.d0
case default
  write(*,*)
  write(*,'("Error(gentpmae): invalid npmae : ",I8)') npmae0
  write(*,*)
  stop
end select
return
end subroutine

