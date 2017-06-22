
! Copyright (C) 2007 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine tm2todm(l,nspinor,k,p,ld,tm2,dmat)
implicit none
! arguments
integer, intent(in) :: l,nspinor
integer, intent(in) :: k,p
integer, intent(in) :: ld
complex(8), intent(in) :: tm2(-ld:ld,-1:1)
complex(8), intent(out) :: dmat(ld,nspinor,ld,nspinor)
! local variables
integer ispn,jspn,x,y
integer m1,m2,lm1,lm2
real(8) nlk,nsp,t1,t2,t3
! external functions
real(8) wigner3j,wigner3jf,factnm
external wigner3j,wigner3jf,factnm
if (l.lt.0) then
  write(*,*)
  write(*,'("Error(tm2todm): l < 0 : ",I8)') l
  write(*,*)
  stop
end if
if ((nspinor.lt.1).or.(nspinor.gt.2)) then
  write(*,*)
  write(*,'("Error(tm2todm): nspinor should be 1 or 2 : ",I8)') nspinor
  write(*,*)
  stop
end if
if (k.lt.0) then
  write(*,*)
  write(*,'("Error(tm2todm): k < 0 : ",I8)') k
  write(*,*)
  stop
end if
if (k.gt.2*l) then
  write(*,*)
  write(*,'("Error(tm2todm): k > 2*l : ",2I8)') k,2*l
  write(*,*)
  stop
end if
if ((p.lt.0).or.(p.gt.1)) then
  write(*,*)
  write(*,'("Error(tm2todm): p should be 0 or 1 : ",I8)') p
  write(*,*)
  stop
end if
! factors n_lk and n_sp; Eq. (21) in article
nlk=factnm(2*l,1)/sqrt(factnm(2*l-k,1)*factnm(2*l+k+1,1))
nsp=1.d0/sqrt(factnm(2+p,1))
t1=dble((2*k+1)*(2*p+1))*nlk*nsp
! compute density matrix from 2-index tensor moment
dmat(:,:,:,:)=0.d0
do x=-k,k
  do y=-p,p
    do ispn=1,nspinor
      do jspn=1,nspinor
        t2=t1*wigner3jf(1,2*p,1,2*jspn-3,2*y,3-2*ispn)
        lm1=l**2
        do m1=-l,l
          lm1=lm1+1
          lm2=l**2
          do m2=-l,l
            lm2=lm2+1
            t3=t2*dble((-1)**(m2-l+jspn-1))*wigner3j(l,k,l,-m2,x,m1)
            dmat(lm1,ispn,lm2,jspn)=dmat(lm1,ispn,lm2,jspn)+t3*tm2(x,y)
          end do
        end do
      end do
    end do
  end do
end do
return
end subroutine

