
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hmlalo(ias,ngp,apwalm,ld,h)
use modmain
implicit none
! arguments
integer, intent(in) :: ias
integer, intent(in) :: ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
integer, intent(in) :: ld
complex(8), intent(inout) :: h(*)
! local variables
integer is,io,ilo
integer l1,l2,l3,m1,m2,m3
integer lm1,lm2,lm3,i,j,k
complex(8) zsum
is=idxis(ias)
do ilo=1,nlorb(is)
  l1=lorbl(ilo,is)
  do m1=-l1,l1
    lm1=idxlm(l1,m1)
    j=ngp+idxlo(lm1,ilo,ias)
    lm3=0
    do l3=0,lmaxmat
      do m3=-l3,l3
        lm3=lm3+1
        do io=1,apword(l3,is)
          zsum=0.d0
          do l2=0,lmaxvr
            if (mod(l1+l2+l3,2).eq.0) then
              do m2=-l2,l2
                lm2=idxlm(l2,m2)
                zsum=zsum+gntyry(lm1,lm2,lm3)*hloa(lm2,io,l3,ilo,ias)
              end do
            end if
          end do
! note that what is actually computed is the Hermitian conjugate of <lo|H|APW>
          if (abs(dble(zsum))+abs(aimag(zsum)).gt.1.d-14) then
            k=(j-1)*ld
            do i=1,ngp
              k=k+1
              h(k)=h(k)+conjg(zsum*apwalm(i,io,lm3,ias))
            end do
          end if
        end do
      end do
    end do
  end do
end do
return
end subroutine

