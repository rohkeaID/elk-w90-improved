
! Copyright (C) 2008 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

! !ROUTINE: dmtotm2
! !INTERFACE:
subroutine dmtotm2(l,nspinor,k,p,ld,dmat,tm2)
! !INPUT/OUTPUT PARAMETERS:
!   l       : angular momentum (in,integer)
!   nspinor : number of spinors (in,integer)
!   k       : k-index of tensor moment (in,integer)
!   p       : p-index of tensor moment (in,integer)
!   ld      : leading dimension (in,integer)
!   dmat    : density matrix (in,complex(ld,nspinor,ls,nspinor))
!   tm2     : 2-index tensor moment (out,complex(-ld:ld,-1:1))
! !DESCRIPTION:
!   Transform the density matrix to a 2-index tensor moment representation, see
!   {\it Phys. Rev. B} {\bf 80}, 035121 (2009).
!
! !REVISION HISTORY:
!   Created April 2008 (F. Cricchio and L. Nordstrom)
!   Modified, January 2014 (JKD)
!EOP
!BOC
implicit none
integer, intent(in) :: l,nspinor
integer, intent(in) :: k,p
integer, intent(in) :: ld
complex(8),intent(in) :: dmat(ld,nspinor,ld,nspinor)
complex(8),intent(out) :: tm2(-ld:ld,-1:1)
! local variables
integer ispn,jspn,x,y
integer m1,m2,lm1,lm2
real(8) nlk,nsp,t1,t2,t3
! external functions
real(8) wigner3j,wigner3jf,factnm,factr
external wigner3j,wigner3jf,factnm,factr
if (l.lt.0) then
  write(*,*)
  write(*,'("Error(dmtotm2): l < 0 : ",I8)') l
  write(*,*)
  stop
end if
if ((nspinor.lt.1).or.(nspinor.gt.2)) then
  write(*,*)
  write(*,'("Error(dmtotm2): nspinor should be 1 or 2 : ",I8)') nspinor
  write(*,*)
  stop
end if
if (k.lt.0) then
  write(*,*)
  write(*,'("Error(dmtotm2): k < 0 : ",I8)') k
  write(*,*)
  stop
end if
if (k.gt.2*l) then
  write(*,*)
  write(*,'("Error(dmtotm2): k > 2*l : ",2I8)') k,2*l
  write(*,*)
  stop
end if
if ((p.lt.0).or.(p.gt.1)) then
  write(*,*)
  write(*,'("Error(dmtotm2): p should be 0 or 1 : ",I8)') p
  write(*,*)
  stop
end if
! calculate the 2-index tensor moment; see Eq. (23) in article
nlk=factnm(2*l,1)/sqrt(factnm(2*l-k,1)*factnm(2*l+k+1,1))
nsp=1.d0/sqrt(factnm(2+p,1))
t1=dble((-1)**(l))/(nlk*nsp)
tm2(:,:)=0.d0
do x=-k,k
  do y=-p,p
    do ispn=1,nspinor
      do jspn=1,nspinor
        t2=t1*wigner3jf(1,2*p,1,2*jspn-3,2*y,3-2*ispn)
        if (abs(t2).gt.1.d-10) then
          lm1=l**2
          do m1=-l,l
            lm1=lm1+1
            lm2=l**2
            do m2=-l,l
              lm2=lm2+1
              t3=t2*dble((-1)**(1+jspn-m2))*wigner3j(l,k,l,-m2,x,m1)
              tm2(x,y)=tm2(x,y)+t3*dmat(lm1,ispn,lm2,jspn)
            end do
          end do
        end if
      end do
    end do
  end do
end do
return
end subroutine
!EOC

