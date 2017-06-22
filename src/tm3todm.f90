
! Copyright (C) 2007 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine tm3todm(l,nspinor,k,p,r,ld,tm3,dmat)
implicit none
! arguments
integer, intent(in) :: l,nspinor
integer, intent(in) :: k,p,r
integer, intent(in) :: ld
complex(8), intent(in) :: tm3(-ld:ld)
complex(8), intent(out) :: dmat(ld,nspinor,ld,nspinor)
! parameters
complex(8), parameter :: zi=(0.d0,1.d0)
! local variables
integer x,y,g,t
real(8) t1,t2
complex(8) z1
! allocatable arrays
complex(8), allocatable :: tm2(:,:)
! external functions
real(8) wigner3j,factnm,factr
external wigner3j,factnm,factr
if (l.lt.0) then
  write(*,*)
  write(*,'("Error(tm3todm): l < 0 : ",I8)') l
  write(*,*)
  stop
end if
if ((nspinor.lt.1).or.(nspinor.gt.2)) then
  write(*,*)
  write(*,'("Error(tm3todm): nspinor should be 1 or 2 : ",I8)') nspinor
  write(*,*)
  stop
end if
if (k.lt.0) then
  write(*,*)
  write(*,'("Error(tm3todm): k < 0 : ",I8)') k
  write(*,*)
  stop
end if
if (k.gt.2*l) then
  write(*,*)
  write(*,'("Error(tm3todm): k > 2*l : ",2I8)') k,2*l
  write(*,*)
  stop
end if
if ((p.lt.0).or.(p.gt.1)) then
  write(*,*)
  write(*,'("Error(tm3todm): p should be 0 or 1 : ",I8)') p
  write(*,*)
  stop
end if
if (r.lt.abs(k-p)) then
  write(*,*)
  write(*,'("Error(tm3todm): r < |k-p| : ",2I8)') r,abs(k-p)
  write(*,*)
  stop
end if
if (r.gt.(k+p)) then
  write(*,*)
  write(*,'("Error(tm3todm): r > k+p : ",2I8)') r,k+p
  write(*,*)
  stop
end if
! compute 2-index tensor moment from 3-index tensor moment
allocate(tm2(-ld:ld,-1:1))
tm2(:,:)=0.d0
do x=-k,k
  do y=-p,p
    t1=dble((-1)**(x+y-k-p))
    g=k+p+r
    if (mod(g,2).eq.0) then
      z1=t1*wigner3j(k,p,r,0,0,0)
    else
      t1=t1*sqrt(factr(g-2*k,g+1)*factnm(g-2*p,1)*factnm(g-2*r,1))
      t1=t1*factnm(g,2)/(factnm(g-2*k,2)*factnm(g-2*p,2)*factnm(g-2*r,2))
      z1=t1*zi**g
    end if
    do t=-r,r
      t2=wigner3j(k,r,p,-x,t,-y)*(2*r+1)
      tm2(x,y)=tm2(x,y)+t2*z1*tm3(t)
    end do
  end do
end do
! compute the density matrix from the 2-index tensor moment
call tm2todm(l,nspinor,k,p,ld,tm2,dmat)
deallocate(tm2)
return
end subroutine

