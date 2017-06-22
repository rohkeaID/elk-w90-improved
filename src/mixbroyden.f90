
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine mixbroyden(iscl,n,msd,alpha,w0,nu,mu,f,df,u,a,d)
implicit none
! arguments
integer, intent(in) :: iscl,n,msd
real(8), intent(in) :: alpha,w0
real(8), intent(inout) :: nu(n),mu(n,2)
real(8), intent(inout) :: f(n,2),df(n,msd)
real(8), intent(inout) :: u(n,msd)
real(8), intent(inout) :: a(msd,msd)
real(8), intent(out) :: d
! local variables
integer jc,kp,kc
integer k,l,m,info
real(8) t1
! automatic arrays
integer ipiv(msd)
real(8) c(msd),beta(msd,msd),gamma(msd)
real(8) work(msd)
! external functions
real(8) ddot,dnrm2
external ddot,dnrm2
if (n.lt.1) then
  write(*,*)
  write(*,'("Error(mixbroyden): n < 1 : ",I8)') n
  write(*,*)
  stop
end if
if (msd.lt.2) then
  write(*,*)
  write(*,'("Error(mixbroyden): msd < 2 : ",I8)') msd
  write(*,*)
  stop
end if
! initialise mixer
if (iscl.le.0) then
  call dcopy(n,nu,1,mu(:,1),1)
  call dcopy(n,nu,1,mu(:,2),1)
  f(:,1)=0.d0
  df(:,1)=0.d0
  u(:,1)=0.d0
  a(:,:)=0.d0
  d=1.d0
  return
end if
! current subspace dimension
m=min(iscl+1,msd)
! current index modulo m
jc=mod(iscl,m)+1
! previous index modulo 2
kp=mod(iscl-1,2)+1
! current index modulo 2
kc=mod(iscl,2)+1
f(:,kc)=nu(:)-mu(:,kp)
d=sum(f(:,kc)**2)
d=sqrt(d/dble(n))
df(:,jc)=f(:,kc)-f(:,kp)
t1=dnrm2(n,df(:,jc),1)
if (t1.gt.1.d-8) t1=1.d0/t1
call dscal(n,t1,df(:,jc),1)
u(:,jc)=alpha*df(:,jc)+t1*(mu(:,kp)-mu(:,kc))
do k=1,m
  c(k)=ddot(n,df(:,k),1,f(:,kc),1)
end do
do k=1,m
  a(k,jc)=ddot(n,df(:,jc),1,df(:,k),1)
  a(jc,k)=a(k,jc)
end do
beta(:,:)=a(:,:)
do k=1,m
  beta(k,k)=beta(k,k)+w0**2
end do
! invert beta
call dgetrf(m,m,beta,msd,ipiv,info)
if (info.eq.0) call dgetri(m,beta,msd,ipiv,work,m,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(mixbroyden): could not invert matrix")')
  write(*,*)
  stop
end if
do l=1,m
  gamma(l)=0.d0
  do k=1,m
    gamma(l)=gamma(l)+c(k)*beta(k,l)
  end do
end do
nu(:)=mu(:,kp)+alpha*f(:,kc)
do l=1,m
  call daxpy(n,-gamma(l),u(:,l),1,nu,1)
end do
call dcopy(n,nu,1,mu(:,kc),1)
return
end subroutine

