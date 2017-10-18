
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genkpakq
use modmain
use modultra
implicit none
! local variables
integer i1,i2,i3,ikpa
integer j1,j2,j3
integer ik,jk
integer iq
real(8) v(3)

!******

! allocatable arrays
real(8), allocatable :: vkl0(:,:),vkc0(:,:),wkpt0(:)
if (.not.ultracell) return
!-----------------------------------!
!     generate the kappa-points     !
!-----------------------------------!
! number of kappa-points
nkpa=ngridkpa(1)*ngridkpa(2)*ngridkpa(3)
! allocate the global kappa-point arrays
if (allocated(ivkpa)) deallocate(ivkpa)
allocate(ivkpa(3,nkpa))
if (allocated(vkpalu)) deallocate(vkpalu)
allocate(vkpalu(3,nkpa))
if (allocated(vkpal)) deallocate(vkpal)
allocate(vkpal(3,nkpa))
if (allocated(vkpac)) deallocate(vkpac)
allocate(vkpac(3,nkpa))
! integer grid intervals for kappa-points
intkpa(1,:)=ngridkpa(:)/2-ngridkpa(:)+1
intkpa(2,:)=ngridkpa(:)/2
! generate the kappa-point grid
ikpa=0
do i3=intkpa(1,3),intkpa(2,3)
  v(3)=dble(i3)/dble(ngridkpa(3))
  do i2=intkpa(1,2),intkpa(2,2)
    v(2)=dble(i2)/dble(ngridkpa(2))
    do i1=intkpa(1,1),intkpa(2,1)
      v(1)=dble(i1)/dble(ngridkpa(1))
      ikpa=ikpa+1
! find index to kappa = 0 point
      if ((i1.eq.0).and.(i2.eq.0).and.(i3.eq.0)) ikpa0=ikpa
! store the kappa-point locations on the integer grid
      ivkpa(1,ikpa)=i1
      ivkpa(2,ikpa)=i2
      ivkpa(3,ikpa)=i3
! store the kappa-point in ultracell lattice coordinates
      vkpalu(:,ikpa)=v(:)
! kappa-point in Cartesian coordinates
      call r3mv(bvecu,v,vkpac(:,ikpa))
! kappa-point in unit cell lattice coordinates
      call r3mv(binv,vkpac(:,ikpa),vkpal(:,ikpa))
    end do
  end do
end do
! store the existing k-point arrays
allocate(vkl0(3,nkpt),vkc0(3,nkpt),wkpt0(nkpt))
vkl0(:,1:nkpt)=vkl(:,1:nkpt)
vkc0(:,1:nkpt)=vkc(:,1:nkpt)
wkpt0(1:nkpt)=wkpt(1:nkpt)
! number of k+kappa-points
nkpt0=nkpt
nkpt=nkpt0*nkpa
! deallocate and reallocate k-point arrays
deallocate(vkl,vkc,wkpt)
allocate(vkl(3,nkpt),vkc(3,nkpt),wkpt(nkpt))
jk=0
do ik=1,nkpt0
  do ikpa=1,nkpa
    jk=jk+1
    vkl(:,jk)=vkl0(:,ik)+vkpal(:,ikpa)
    vkc(:,jk)=vkc0(:,ik)+vkpac(:,ikpa)
    wkpt(jk)=wkpt0(ik)/dble(nkpa)
  end do
end do
deallocate(vkl0,vkc0,wkpt0)
!-------------------------------!
!     generate the Q-points     !
!-------------------------------!
! Q-point grid should be twice the kappa-point grid
ngridq(:)=2*ngridkpa(:)
! find next largest FFT-compatible grid size
call nfftifc(ngridq(1))
call nfftifc(ngridq(2))
call nfftifc(ngridq(3))
! total number of Q-points
nqpt=ngridq(1)*ngridq(2)*ngridq(3)
! integer grid intervals for the Q-points
intq(1,:)=ngridq(:)/2-ngridq(:)+1
intq(2,:)=ngridq(:)/2
! allocate global Q-point arrays
if (allocated(ivqiq)) deallocate(ivqiq)
allocate(ivqiq(intq(1,1):intq(2,1),intq(1,2):intq(2,2),intq(1,3):intq(2,3)))
if (allocated(vqc)) deallocate(vqc)
allocate(vqc(3,nqpt))
! generate the FFT-compatible Q-point grid
do i1=intq(1,1),intq(2,1)
  if (i1.ge.0) then
    j1=i1
  else
    j1=ngridq(1)+i1
  end if
  v(1)=dble(i1)/dble(ngridkpa(1))
  do i2=intq(1,2),intq(2,2)
    if (i2.ge.0) then
      j2=i2
    else
      j2=ngridq(2)+i2
    end if
    v(2)=dble(i2)/dble(ngridkpa(2))
    do i3=intq(1,3),intq(2,3)
      if (i3.ge.0) then
        j3=i3
      else
        j3=ngridq(3)+i3
      end if
      v(3)=dble(i3)/dble(ngridkpa(3))
      iq=j3*ngridq(2)*ngridq(1)+j2*ngridq(1)+j1+1
! map from (i1,i2,i3) to Q-point index
      ivqiq(i1,i2,i3)=iq
! store the Q-point in Cartesian coordinates
      call r3mv(bvecu,v,vqc(:,iq))
    end do
  end do
end do
return
end subroutine

