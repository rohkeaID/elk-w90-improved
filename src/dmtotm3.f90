
! Copyright (C) 2008 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

! !ROUTINE: dmtotm3
! !INTERFACE:
subroutine dmtotm3(l,nspinor,k,p,r,ld,dmat,tm3)
! !INPUT/OUTPUT PARAMETERS:
!   l       : angular momentum (in,integer)
!   nspinor : number of spinor components (in,integer)
!   k       : k-index of tensor moment (in,integer)
!   p       : p-index of tensor moment (in,integer)
!   r       : r-index of tensor moment (in,integer)
!   ld      : leading dimension (in,integer)
!   dmat    : density matrix (in,complex(ld,nspinor,ld,nspinor))
!   tm3     : 3-index spherical tensor moment (out,complex(-ld:ld))
! !DESCRIPTION:
!   Transform the density matrix to a 3-index spherical tensor moment
!   representation, see {\it Phys. Rev. B} {\bf 80}, 035121 (2009).
!
! !REVISION HISTORY:
!   Created April 2008 (F. Cricchio and L. Nordstrom)
!   Modified, January 2014 (JKD)
!EOP
!BOC
implicit none
integer, intent(in) :: l,nspinor
integer, intent(in) :: k,p,r
integer, intent(in) :: ld
complex(8), intent(in) :: dmat(ld,nspinor,ld,nspinor)
complex(8), intent(out) :: tm3(-ld:ld)
! local variables
integer ispn,jspn,g,t,x,y
integer m1,m2,lm1,lm2
real(8) nlk,nsp,t1,t2,t3
complex(8), parameter :: zi=(0.d0,1.d0)
complex(8) z1
! external functions
real(8) wigner3j,wigner3jf,factnm,factr
external wigner3j,wigner3jf,factnm,factr
if (l.lt.0) then
  write(*,*)
  write(*,'("Error(dmtotm3): l < 0 : ",I8)') l
  write(*,*)
  stop
end if
if ((nspinor.lt.1).or.(nspinor.gt.2)) then
  write(*,*)
  write(*,'("Error(dmtotm3): nspinor should be 1 or 2 : ",I8)') nspinor
  write(*,*)
  stop
end if
if (k.lt.0) then
  write(*,*)
  write(*,'("Error(dmtotm3): k < 0 : ",I8)') k
  write(*,*)
  stop
end if
if (k.gt.2*l) then
  write(*,*)
  write(*,'("Error(dmtotm3): k > 2*l : ",2I8)') k,2*l
  write(*,*)
  stop
end if
if ((p.lt.0).or.(p.gt.1)) then
  write(*,*)
  write(*,'("Error(dmtotm3): p should be 0 or 1 : ",I8)') p
  write(*,*)
  stop
end if
if (r.lt.abs(k-p)) then
  write(*,*)
  write(*,'("Error(dmtotm3): r < |k-p| : ",2I8)') r,abs(k-p)
  write(*,*)
  stop
end if
if (r.gt.(k+p)) then
  write(*,*)
  write(*,'("Error(dmtotm3): r > k+p : ",2I8)') r,k+p
  write(*,*)
  stop
end if
! calculate the 3-index tensor moment; see Eqs. (23), (26), (27) in article
g=k+p+r
if (mod(g,2).eq.0) then
  z1=1.d0/wigner3j(k,p,r,0,0,0)
else
  t1=sqrt(factr(g+1,g-2*k)/(factnm(g-2*p,1)*factnm(g-2*r,1)))
  t1=t1*factnm(g-2*k,2)*factnm(g-2*p,2)*factnm(g-2*r,2)/factnm(g,2)
  z1=t1*zi**(-g)
end if
nlk=factnm(2*l,1)/sqrt(factnm(2*l-k,1)*factnm(2*l+k+1,1))
nsp=1.d0/sqrt(factnm(2+p,1))
z1=z1*dble((-1)**(k+p+l))/(nlk*nsp)
tm3(:)=0.d0
do t=-r,r
  do x=-k,k
    do y=-p,p
      t1=dble((-1)**(x+y))*wigner3j(k,r,p,-x,t,-y)
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
                tm3(t)=tm3(t)+t3*z1*dmat(lm1,ispn,lm2,jspn)
              end do
            end do
          end if
        end do
      end do
    end do
  end do
end do
return
end subroutine
!EOC

