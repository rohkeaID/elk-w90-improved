
! Copyright (C) 2014 L. Nordstrom, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine ftmfield
use modmain
use moddftu
implicit none
! local variables
integer is,ia,ias,n2,i
integer l,n,k,p,r,t,x,y
complex(8) z1
! allocatable arrays
complex(8), allocatable :: tm2(:,:),tm3(:)
complex(8), allocatable :: dmat(:,:,:,:)
if (ftmtype.le.0) return
if (mod(iscl,ftmstep).ne.1) return
allocate(tm2(-lmmaxdm:lmmaxdm,-1:1),tm3(-lmmaxdm:lmmaxdm))
allocate(dmat(lmmaxdm,nspinor,lmmaxdm,nspinor))
n2=(lmmaxdm*nspinor)**2
! loop over FTM entries
do i=1,ntmfix
  is=itmfix(1,i)
  ia=itmfix(2,i)
  ias=idxas(ia,is)
  l=itmfix(3,i)
  n=itmfix(4,i)
  k=itmfix(5,i)
  p=itmfix(6,i)
  if (n.eq.2) then
    x=itmfix(7,i)
    y=itmfix(8,i)
! decompose density matrix in 2-index tensor moment components
    call dmtotm2(l,nspinor,k,p,lmmaxdm,dmatmt(:,:,:,:,ias),tm2)
! take difference between current and target moment
    z1=tm2(x,y)-tmfix(i)
    tm2(:,:)=0.d0
    tm2(x,y)=tauftm*z1
! compute new density matrix
    call tm2todm(l,nspinor,k,p,lmmaxdm,tm2,dmat)
  else
    r=itmfix(7,i)
    t=itmfix(8,i)
! decompose density matrix in 3-index tensor moment components
    call dmtotm3(l,nspinor,k,p,r,lmmaxdm,dmatmt(:,:,:,:,ias),tm3)
! take difference between current and target moment
    z1=tm3(t)-tmfix(i)
    tm3(:)=0.d0
    tm3(t)=tauftm*z1
! compute new density matrix
    call tm3todm(l,nspinor,k,p,r,lmmaxdm,tm3,dmat)
  end if
! make density matrix Hermitian
  call hrmdmat(lmmaxdm*nspinor,dmat)
! apply rotation matrices and add to potential matrix global arrays
  call rotdmat(rtmfix(:,:,1,i),rtmfix(:,:,2,i),lmaxdm,nspinor,lmmaxdm,dmat, &
   vmftm(:,:,:,:,ias))
  call zaxpy(n2,zone,vmftm(:,:,:,:,ias),1,vmatmt(:,:,:,:,ias),1)
end do
deallocate(tm2,tm3,dmat)
return
end subroutine

