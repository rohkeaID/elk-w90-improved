
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeftm
use modmain
use moddftu
implicit none
! local variables
integer is,ia,ias,i
integer l,n,k,p,r,t,x,y
! allocatable arrays
complex(8), allocatable :: tm2(:,:),tm3(:)
allocate(tm2(-lmmaxdm:lmmaxdm,-1:1),tm3(-lmmaxdm:lmmaxdm))
! open FTM.OUT
open(50,file='FTM.OUT',action='WRITE',form='FORMATTED')
do i=1,ntmfix
  is=itmfix(1,i)
  ia=itmfix(2,i)
  ias=idxas(ia,is)
  l=itmfix(3,i)
  n=itmfix(4,i)
  k=itmfix(5,i)
  p=itmfix(6,i)
  write(50,*)
  write(50,'("Species : ",I4," (",A,"), atom : ",I4)') is,trim(spsymb(is)),ia
  write(50,'(" l = ",I2,", n = ",I2)') l,n
  if (n.eq.2) then
    x=itmfix(7,i)
    y=itmfix(8,i)
    write(50,'(" k = ",I2,", p = ",I2,", x = ",I2,", y = ",I2)') k,p,x,y
! decompose density matrix in 2-index tensor moment components
    call dmtotm2(l,nspinor,k,p,lmmaxdm,dmatmt(:,:,:,:,ias),tm2)
    write(50,'(" tensor moment")')
    write(50,'("  current : ",2G18.10)') tm2(x,y)
    write(50,'("  target  : ",2G18.10)') tmfix(i)
  else
    r=itmfix(7,i)
    t=itmfix(8,i)
    write(50,'(" k = ",I2,", p = ",I2,", r = ",I2,", t = ",I2)') k,p,r,t
! decompose density matrix in 3-index tensor moment components
    call dmtotm3(l,nspinor,k,p,r,lmmaxdm,dmatmt(:,:,:,:,ias),tm3)
    write(50,'(" tensor moment")')
    write(50,'("  current : ",2G18.10)') tm3(t)
    write(50,'("  target  : ",2G18.10)') tmfix(i)
  end if
end do
close(50)
deallocate(tm2,tm3)
return
end subroutine

