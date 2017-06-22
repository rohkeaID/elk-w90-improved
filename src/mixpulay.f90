
! Copyright (C) 2008 S. Suehara. This file is distributed under the terms of the
! GNU Lesser General Public License. See the file COPYING for license details.

!BOP
! !ROUTINE: mixpulay
! !INTERFACE:
subroutine mixpulay(iscl,n,msd,nu,mu,f,d)
! !INPUT/OUTPUT PARAMETERS:
!   iscl : self-consistent loop number (in,integer)
!   n    : vector length (in,integer)
!   msd  : maximum subspace dimension (in,integer)
!   nu   : current output vector as well as the next input vector of the
!          self-consistent loop (inout,real(n))
!   mu   : used internally (inout,real(n,msd))
!   f    : used internally (inout,real(n,msd))
!   d    : RMS difference between old and new output vectors (out,real)
! !DESCRIPTION:
!   Pulay's mixing scheme which uses direct inversion in the iterative subspace
!   (DIIS). See {\it Chem. Phys. Lett.} {\bf 73}, 393 (1980).
!
! !REVISION HISTORY:
!   Created October 2008 (S. Suehara, NIMS)
!   Modified, October 2008 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: iscl
integer, intent(in) :: n
integer, intent(in) :: msd
real(8), intent(inout) :: nu(n)
real(8), intent(inout) :: mu(n,msd)
real(8), intent(inout) :: f(n,msd)
real(8), intent(out) :: d
! local variables
integer jc,jn,i,j,m,info
! initial mixing parameter
real(8), parameter :: beta=0.25d0
! allocatable arrays
integer, allocatable :: ipiv(:)
real(8), allocatable :: alpha(:),a(:,:),work(:)
! external functions
real(8) ddot
external ddot
if (n.lt.1) then
  write(*,*)
  write(*,'("Error(mixpulay): n < 1 : ",I8)') n
  write(*,*)
  stop
end if
if (msd.lt.2) then
  write(*,*)
  write(*,'("Error(mixpulay): msd < 2 : ",I8)') msd
  write(*,*)
  stop
end if
if (iscl.le.0) then
  mu(:,1)=nu(:)
  f(:,1)=0.d0
  d=1.d0
  return
end if
if (iscl.le.(msd-1)) then
  nu(:)=beta*nu(:)+(1.d0-beta)*mu(:,iscl)
  f(:,iscl+1)=nu(:)-mu(:,iscl)
  mu(:,iscl+1)=nu(:)
  d=sum(f(:,iscl+1)**2)
  d=sqrt(d/dble(n))
  return
end if
! current index
jc=mod(iscl,msd)+1
! next index
jn=mod(iscl+1,msd)+1
! matrix size
m=min(iscl+1,msd)+1
allocate(ipiv(m),alpha(m),a(m,m),work(m))
! compute f and RMS difference
f(:,jc)=nu(:)-mu(:,jc)
d=sum(f(:,jc)**2)
d=sqrt(d/dble(n))
! solve the linear system
a(:,:)=0.d0
do i=1,m-1
  do j=i,m-1
    a(i,j)=a(i,j)+ddot(n,f(:,i),1,f(:,j),1)
  end do
  a(i,m)=1.d0
end do
alpha(:)=0.d0
alpha(m)=1.d0
call dsysv('U',m,1,a,m,ipiv,alpha,m,work,m,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(mixpulay): could not solve linear system")')
  write(*,'(" DSYSV returned INFO = ",I8)') info
  write(*,*)
  stop
end if
nu(:)=0.d0
do i=1,m-1
  nu(:)=nu(:)+alpha(i)*(mu(:,i)+f(:,i))
end do
mu(:,jn)=nu(:)
deallocate(ipiv,alpha,a,work)
return
end subroutine
!EOC

