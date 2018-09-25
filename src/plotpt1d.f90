
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: plotpt1d
! !INTERFACE:
subroutine plotpt1d(cvec,nv,np,vvl,vpl,dv,dp)
! !INPUT/OUTPUT PARAMETERS:
!   cvec : matrix of (reciprocal) lattice vectors stored column-wise
!         (in,real(3,3))
!   nv   : number of vertices (in,integer)
!   np   : number of connecting points (in,integer)
!   vvl  : vertex vectors in lattice coordinates (in,real(3,nv))
!   vpl  : connecting point vectors in lattice coordinates (out,real(3,np))
!   dv   : cummulative distance to each vertex (out,real(nv))
!   dp   : cummulative distance to each connecting point (out,real(np))
! !DESCRIPTION:
!   Generates a set of points which interpolate between a given set of vertices.
!   Vertex points are supplied in lattice coordinates in the array {\tt vvl} and
!   converted to Cartesian coordinates with the matrix {\tt cvec}. Interpolating
!   points are stored in the array {\tt vpl}. The cummulative distances to the
!   vertices and points along the path are stored in arrays {\tt dv} and
!   {\tt dp}, respectively.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!   Improved September 2007 (JKD)
!   Improved again, July 2010 (T. McQueen and JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: cvec(3,3)
integer, intent(in) :: nv,np
real(8), intent(in) :: vvl(3,nv)
real(8), intent(out) :: vpl(3,np),dv(nv),dp(np)
! local variables
integer i,j,k,m,n
real(8) vl(3),vc(3)
real(8) dt,f,t1
! alloctable arrays
real(8), allocatable :: seg(:)
if (nv.lt.1) then
  write(*,*)
  write(*,'("Error(plotpt1d): nv < 1 : ",I8)') nv
  write(*,*)
  stop
end if
if (np.lt.nv) then
  write(*,*)
  write(*,'("Error(plotpt1d): np < nv : ",2I8)') np,nv
  write(*,*)
  stop
end if
! special case of 1 vertex
if (nv.eq.1) then
  dv(1)=0.d0
  dp(:)=0.d0
  do i=1,np
    vpl(:,i)=vvl(:,1)
  end do
  return
end if
allocate(seg(nv))
! find the length of each segment and total distance
dt=0.d0
do i=1,nv-1
  dv(i)=dt
  vl(:)=vvl(:,i+1)-vvl(:,i)
  call r3mv(cvec,vl,vc)
  seg(i)=sqrt(vc(1)**2+vc(2)**2+vc(3)**2)
  dt=dt+seg(i)
end do
dv(nv)=dt
! add small amount to total distance to avoid 0/0 condition
dt=dt+1.d-8
! number of points to use between vertices
n=np-nv
! construct the interpolating path
k=0
do i=1,nv-1
  t1=dble(n)*seg(i)/dt
  m=nint(t1)
  if ((m.gt.n).or.(i.eq.(nv-1))) m=n
  do j=1,m+1
    k=k+1
    f=dble(j-1)/dble(m+1)
    dp(k)=dv(i)+f*seg(i)
    vpl(:,k)=vvl(:,i)*(1.d0-f)+vvl(:,i+1)*f
  end do
  dt=dt-seg(i)
  n=n-m
end do
dp(np)=dv(nv)
vpl(:,np)=vvl(:,nv)
deallocate(seg)
return
end subroutine
!EOC

