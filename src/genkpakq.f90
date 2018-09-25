
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genkpakq
use modmain
use modulr
implicit none
! local variables
integer i1,i2,i3,j1,j2,j3
integer ikpa,ik,jk,iq,jq
real(8) v1(3),v2(3),v3(3)
! allocatable arrays
integer, allocatable :: idx(:),ivq0(:,:)
real(8), allocatable :: vqc0(:,:),qc0(:)
real(8), allocatable :: vkl0(:,:),vkc0(:,:)
! determine the kappa-point grid sizes
ngridkpa(:)=int(kpamax*sqrt(avecu(1,:)**2+avecu(2,:)**2+avecu(3,:)**2)/pi)+1
! Q-point grid should be twice the kappa-point grid
ngridq(:)=2*ngridkpa(:)-1
! find next largest FFT-compatible grid size
call nfftifc(ngridq(1))
call nfftifc(ngridq(2))
call nfftifc(ngridq(3))
! total number of Q-points
nqpt=ngridq(1)*ngridq(2)*ngridq(3)
! integer grid intervals for the Q-points
intq(1,:)=ngridq(:)/2-ngridq(:)+1
intq(2,:)=ngridq(:)/2
! allocate local arrays
allocate(idx(nqpt),ivq0(3,nqpt))
allocate(vqc0(3,nqpt),qc0(nqpt))
! generate Q-point grid
iq=0
do i1=intq(1,1),intq(2,1)
  v1(:)=dble(i1)*bvecu(:,1)
  do i2=intq(1,2),intq(2,2)
    v2(:)=v1(:)+dble(i2)*bvecu(:,2)
    do i3=intq(1,3),intq(2,3)
      v3(:)=v2(:)+dble(i3)*bvecu(:,3)
      iq=iq+1
      ivq0(1,iq)=i1
      ivq0(2,iq)=i2
      ivq0(3,iq)=i3
! store the Q-point in Cartesian coordinates
      vqc0(:,iq)=v3(:)
! length of Q-vector
      qc0(iq)=sqrt(v3(1)**2+v3(2)**2+v3(3)**2)
    end do
  end do
end do
! sort the Q-points according to length
call sortidx(nqpt,qc0,idx)
! allocate global Q-point arrays
if (allocated(ivq)) deallocate(ivq)
allocate(ivq(3,nqpt))
if (allocated(ivqiq)) deallocate(ivqiq)
allocate(ivqiq(intq(1,1):intq(2,1),intq(1,2):intq(2,2),intq(1,3):intq(2,3)))
if (allocated(iqfft)) deallocate(iqfft)
allocate(iqfft(nqpt))
if (allocated(vql)) deallocate(vql)
allocate(vql(3,nqpt))
if (allocated(vqc)) deallocate(vqc)
allocate(vqc(3,nqpt))
if (allocated(qc)) deallocate(qc)
allocate(qc(nqpt))
! reorder and store in global arrays
do iq=1,nqpt
  jq=idx(iq)
  ivq(:,iq)=ivq0(:,jq)
  vqc(:,iq)=vqc0(:,jq)
  qc(iq)=qc0(jq)
end do
! determine the Q-vectors in lattice coordinates
do iq=1,nqpt
  call r3mv(binv,vqc(:,iq),vql(:,iq))
end do
deallocate(idx,ivq0,vqc0,qc0)
do iq=1,nqpt
  i1=ivq(1,iq)
  i2=ivq(2,iq)
  i3=ivq(3,iq)
! map from (i1,i2,i3) to Q-vector index
  ivqiq(i1,i2,i3)=iq
! Fourier transform index
  if (i1.ge.0) then
    j1=i1
  else
    j1=ngridq(1)+i1
  end if
  if (i2.ge.0) then
    j2=i2
  else
    j2=ngridq(2)+i2
  end if
  if (i3.ge.0) then
    j3=i3
  else
    j3=ngridq(3)+i3
  end if
  iqfft(iq)=j3*ngridq(2)*ngridq(1)+j2*ngridq(1)+j1+1
end do
! determine the number of kappa-points
nkpa=0
do iq=1,nqpt
  if (qc(iq).gt.kpamax) exit
  nkpa=nkpa+1
end do
! store the existing k-point arrays
allocate(vkl0(3,nkpt),vkc0(3,nkpt))
vkl0(:,1:nkpt)=vkl(:,1:nkpt)
vkc0(:,1:nkpt)=vkc(:,1:nkpt)
! number of k+kappa-points
nkpt0=nkpt
nkpt=nkpt0*nkpa
! deallocate and reallocate k-point arrays
deallocate(vkl,vkc)
allocate(vkl(3,nkpt),vkc(3,nkpt))
jk=0
do ik=1,nkpt0
  do ikpa=1,nkpa
    jk=jk+1
    vkl(:,jk)=vkl0(:,ik)+vql(:,ikpa)
    vkc(:,jk)=vkc0(:,ik)+vqc(:,ikpa)
  end do
end do
deallocate(vkl0,vkc0)
return
end subroutine

