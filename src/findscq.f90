
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine findscq(iq,avec0,nsc,vsc)
use modmain
implicit none
! arguments
integer, intent(in) :: iq
real(8), intent(in) :: avec0(3,3)
integer, intent(out) :: nsc
real(8), intent(out) :: vsc(3,nqptnr)
! local variables
integer i1,i2,i3
integer scl(3,3),i,n
real(8) dmin,t1
real(8) v1(3),v2(3)
! check for Gamma-point phonon
if (iq.eq.iq0) then
  scl(:,:)=0
  scl(1,1)=1
  scl(2,2)=1
  scl(3,3)=1
  nsc=1
  goto 10
end if
! find the first lattice vector
dmin=1.d8
do i1=-ngridq(1),ngridq(1)
  do i2=-ngridq(2),ngridq(2)
    do i3=-ngridq(3),ngridq(3)
      t1=dble(i1)*vql(1,iq)+dble(i2)*vql(2,iq)+dble(i3)*vql(3,iq)
      if (abs(t1-nint(t1)).lt.epslat) then
        v1(:)=dble(i1)*avec0(:,1)+dble(i2)*avec0(:,2)+dble(i3)*avec0(:,3)
        t1=sqrt(v1(1)**2+v1(2)**2+v1(3)**2)
        if ((t1.lt.dmin).and.(t1.gt.epslat)) then
          scl(1,1)=i1
          scl(2,1)=i2
          scl(3,1)=i3
          dmin=t1
        end if
      end if
    end do
  end do
end do
! find the second lattice vector
dmin=1.d8
do i1=-ngridq(1),ngridq(1)
  do i2=-ngridq(2),ngridq(2)
    do i3=-ngridq(3),ngridq(3)
      t1=dble(i1)*vql(1,iq)+dble(i2)*vql(2,iq)+dble(i3)*vql(3,iq)
      if (abs(t1-nint(t1)).lt.epslat) then
! area defined by first two lattice vectors
        n=(i2*scl(3,1)-i3*scl(2,1))**2 &
         +(i3*scl(1,1)-i1*scl(3,1))**2 &
         +(i1*scl(2,1)-i2*scl(1,1))**2
        if (n.ne.0) then
          v1(:)=dble(i1)*avec0(:,1)+dble(i2)*avec0(:,2)+dble(i3)*avec0(:,3)
          t1=v1(1)**2+v1(2)**2+v1(3)**2
          if (t1.lt.dmin) then
            scl(1,2)=i1
            scl(2,2)=i2
            scl(3,2)=i3
            dmin=t1
          end if
        end if
      end if
    end do
  end do
end do
! find the third lattice vector
nsc=0
dmin=1.d8
do i1=-ngridq(1),ngridq(1)
  do i2=-ngridq(2),ngridq(2)
    do i3=-ngridq(3),ngridq(3)
      t1=dble(i1)*vql(1,iq)+dble(i2)*vql(2,iq)+dble(i3)*vql(3,iq)
      if (abs(t1-nint(t1)).lt.epslat) then
! number of primitive unit cells in supercell
        n=scl(1,2)*(i2*scl(3,1)-i3*scl(2,1)) &
         +scl(2,2)*(i3*scl(1,1)-i1*scl(3,1)) &
         +scl(3,2)*(i1*scl(2,1)-i2*scl(1,1))
        if (n.ne.0) then
          v1(:)=dble(i1)*avec0(:,1)+dble(i2)*avec0(:,2)+dble(i3)*avec0(:,3)
          t1=v1(1)**2+v1(2)**2+v1(3)**2
          if (t1.lt.dmin) then
            nsc=abs(n)
            scl(1,3)=i1
            scl(2,3)=i2
            scl(3,3)=i3
            dmin=t1
          end if
        end if
      end if
    end do
  end do
end do
if (nsc.eq.0) goto 30
10 continue
! new lattice vectors
do i=1,3
  avec(:,i)=dble(scl(1,i))*avec0(:,1) &
           +dble(scl(2,i))*avec0(:,2) &
           +dble(scl(3,i))*avec0(:,3)
end do
! inverse of lattice vector matrix
call r3minv(avec,ainv)
! generate offset vectors for each primitive cell in the supercell
n=1
vsc(:,1)=0.d0
do i1=-ngridq(1),ngridq(1)
  do i2=-ngridq(2),ngridq(2)
    do i3=-ngridq(3),ngridq(3)
      if (n.eq.nsc) return
      v1(:)=dble(i1)*avec0(:,1)+dble(i2)*avec0(:,2)+dble(i3)*avec0(:,3)
      call r3mv(ainv,v1,v2)
      call r3frac(epslat,v2)
      call r3mv(avec,v2,v1)
      do i=1,n
        t1=abs(v1(1)-vsc(1,i))+abs(v1(2)-vsc(2,i))+abs(v1(3)-vsc(3,i))
        if (t1.lt.epslat) goto 20
      end do
      n=n+1
      vsc(:,n)=v1(:)
20 continue
    end do
  end do
end do
30 continue
write(*,*)
write(*,'("Error(findscq): unable to generate supercell")')
write(*,*)
stop
end subroutine

