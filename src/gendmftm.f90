
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma, E. K. U. Gross and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gendmftm
use modmain
use moddftu
implicit none
! local variables
integer is,ia,ias,i
integer l,n,k,p,r,t,x,y
! allocatable arrays
complex(8), allocatable :: tm2(:,:),tm3(:)
complex(8), allocatable :: dmat(:,:,:,:)
! external functions
real(8) dznrm2
external dznrm2
if (ftmtype.eq.0) return
! allocate global array
if (allocated(dmftm)) deallocate(dmftm)
allocate(dmftm(lmmaxdm,nspinor,lmmaxdm,nspinor,natmtot))
! allocate local arrays
allocate(tm2(-lmmaxdm:lmmaxdm,-1:1),tm3(-lmmaxdm:lmmaxdm))
allocate(dmat(lmmaxdm,nspinor,lmmaxdm,nspinor))
! zero the fixed tensor moment density matrices
dmftm(:,:,:,:,:)=0.d0
do i=1,ntmfix
  is=itmfix(1,i)
  if (is.gt.nspecies) then
    write(*,*)
    write(*,'("Error(gendmftm): invalid species number : ",I8)') is
    write(*,*)
    stop
  end if
  ia=itmfix(2,i)
  if (ia.gt.natoms(is)) then
    write(*,*)
    write(*,'("Error(gendmftm): invalid atom number : ",I8)') ia
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  ias=idxas(ia,is)
  l=itmfix(3,i)
  if (l.gt.lmaxdm) then
    write(*,*)
    write(*,'("Error(gendmftm): l > lmaxdm ",2I8)') l,lmaxdm
    write(*,'(" for species ",I4," and atom ",I4)') is,ia
    write(*,*)
    stop
  end if
  n=itmfix(4,i)
  k=itmfix(5,i)
  p=itmfix(6,i)
  if (n.eq.2) then
! generate the 2-index density matrix
    x=itmfix(7,i)
    y=itmfix(8,i)
    if ((abs(x).gt.lmmaxdm).or.(abs(y).gt.1)) then
      write(*,*)
      write(*,'("Error(gendmftm): invalid x or y : ",2I8)') x,y
      write(*,'(" for tensor moment entry ",I3)') i
      write(*,*)
      stop
    end if
    tm2(:,:)=0.d0
    tm2(x,y)=tmfix(i)
    call tm2todm(l,nspinor,k,p,lmmaxdm,tm2,dmat)
  else
! generate the 3-index density matrix
    r=itmfix(7,i)
    t=itmfix(8,i)
    if (abs(t).gt.lmmaxdm) then
      write(*,*)
      write(*,'("Error(gendmftm): invalid t : ",I8)') t
      write(*,'(" for tensor moment entry ",I3)') i
      write(*,*)
      stop
    end if
    tm3(:)=0.d0
    tm3(t)=tmfix(i)
    call tm3todm(l,nspinor,k,p,r,lmmaxdm,tm3,dmat)
  end if
! make density matrix Hermitian
  call hrmdmat(lmmaxdm*nspinor,dmat)
! apply rotation matrices and add to density matrix global array
  call rotdmat(rtmfix(:,:,1,i),rtmfix(:,:,2,i),lmaxdm,nspinor,lmmaxdm,dmat, &
   dmftm(:,:,:,:,ias))
end do
deallocate(tm2,tm3,dmat)
return
end subroutine

