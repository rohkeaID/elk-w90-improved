
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: hmlaa
! !INTERFACE:
subroutine hmlaa(ias,ngp,apwalm,ld,h)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ias    : joint atom and species number (in,integer)
!   ngp    : number of G+p-vectors (in,integer)
!   apwalm : APW matching coefficients
!            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
!   ld     : leading dimension of h (in,integer)
!   h      : Hamiltonian matrix (inout,complex(*))
! !DESCRIPTION:
!   Calculates the APW-APW contribution to the Hamiltonian matrix.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ias,ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
integer, intent(in) :: ld
complex(8), intent(inout) :: h(*)
! local variables
integer is,io,jo
integer l1,l2,l3,m1,m2,m3
integer lm1,lm2,lm3
real(8) t1
complex(8) z1,zsum
! automatic arrays
complex(8) x(ngp),y(ngp)
is=idxis(ias)
lm1=0
do l1=0,lmaxmat
  do m1=-l1,l1
    lm1=lm1+1
    do io=1,apword(l1,is)
      y(:)=0.d0
      lm3=0
      do l3=0,lmaxmat
        do m3=-l3,l3
          lm3=lm3+1
          do jo=1,apword(l3,is)
            zsum=0.d0
            do l2=0,lmaxvr
              if (mod(l1+l2+l3,2).eq.0) then
                do m2=-l2,l2
                  lm2=idxlm(l2,m2)
                  zsum=zsum+gntyry(lm1,lm2,lm3)*haa(lm2,jo,l3,io,l1,ias)
                end do
              end if
            end do
            if (abs(dble(zsum))+abs(aimag(zsum)).gt.1.d-14) then
              call zaxpy(ngp,zsum,apwalm(:,jo,lm3,ias),1,y,1)
            end if
          end do
        end do
      end do
      x(1:ngp)=conjg(apwalm(1:ngp,io,lm1,ias))
      call zher2i(ngp,zone,x,y,ld,h)
    end do
  end do
end do
! kinetic surface contribution
t1=0.5d0*rmt(is)**2
lm1=0
do l1=0,lmaxmat
  do m1=-l1,l1
    lm1=lm1+1
    do io=1,apword(l1,is)
      x(1:ngp)=conjg(apwalm(1:ngp,io,lm1,ias))
      do jo=1,apword(l1,is)
        z1=t1*apwfr(nrmt(is),1,io,l1,ias)*apwdfr(jo,l1,ias)
        call zher2i(ngp,z1,x,apwalm(:,jo,lm1,ias),ld,h)
      end do
    end do
  end do
end do
return

contains

subroutine zher2i(n,alpha,x,y,ld,a)
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(in) :: alpha
complex(8), intent(in) :: x(n),y(n)
integer, intent(in) :: ld
complex(8), intent(inout) :: a(*)
! local variables
integer j,k
! numbers less than eps are considered to be zero
real(8), parameter :: eps=1.d-10
real(8) a1,b1
complex(8) z1
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,z1,a1,b1)
!$OMP DO
do j=1,n
  k=(j-1)*ld
  z1=alpha*y(j)
  if (abs(dble(z1)).gt.eps) then
    if (abs(aimag(z1)).gt.eps) then
! complex prefactor
      call zaxpy(j-1,z1,x,1,a(k+1),1)
      a(k+j)=dble(a(k+j))+dble(z1*x(j))
    else
! real prefactor
      a1=dble(z1)
      call daxpy(2*(j-1),a1,x,1,a(k+1),1)
      a(k+j)=dble(a(k+j))+a1*dble(x(j))
    end if
  else if (abs(aimag(z1)).gt.eps) then
! imaginary prefactor
    b1=aimag(z1)
    a(k+1:k+j-1)=a(k+1:k+j-1)+b1*cmplx(-aimag(x(1:j-1)),dble(x(1:j-1)),8)
    a(k+j)=dble(a(k+j))-b1*aimag(x(j))
  end if
end do
!$OMP END DO
!$OMP END PARALLEL
return
end subroutine

end subroutine
!EOC

