
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeevbse
use modmain
implicit none
! local variables
integer ik,a
integer iostat,nmbse_
integer lwork,info
! allocatable arrays
real(8), allocatable :: rwork(:)
complex(8), allocatable :: w(:),vl(:,:),vr(:,:)
complex(8), allocatable :: work(:)
! initialise global variables
call init0
call init1
! read Fermi energy from a file
call readfermi
! get the eigenvalues from file
do ik=1,nkpt
  call getevalsv(filext,vkl(:,ik),evalsv(:,ik))
end do
! generate the BSE state index arrays
call genidxbse
! allocate global BSE arrays
if (allocated(evalbse)) deallocate(evalbse)
allocate(evalbse(nmbse))
if (allocated(hmlbse)) deallocate(hmlbse)
allocate(hmlbse(nmbse,nmbse))
! read in BSE Hamiltonian matrix
open(50,file='HMLBSE.OUT',action='READ',form='UNFORMATTED',status='OLD', &
 iostat=iostat)
if (iostat.ne.0) then
  write(*,*)
  write(*,'("Error(writeevbse): error opening HMLBSE.OUT")')
  write(*,*)
  stop
end if
read(50) nmbse_
if (nmbse.ne.nmbse_) then
  write(*,*)
  write(*,'("Error(writeevbse): differing nmbse")')
  write(*,'(" current    : ",I6)') nmbse
  write(*,'(" HMLBSE.OUT : ",I6)') nmbse_
  write(*,*)
  stop
end if
read(50) hmlbse
close(50)
write(*,*)
write(*,'("Info(writeevbse): diagonalising the BSE Hamiltonian matrix")')
if (bsefull) then
! full non-Hermitian matrix
  allocate(w(nmbse))
  allocate(vl(1,1),vr(nmbse,nmbse))
  lwork=2*nmbse
  allocate(rwork(lwork),work(lwork))
  call zgeev('N','V',nmbse,hmlbse,nmbse,w,vl,1,vr,nmbse,work,lwork,rwork,info)
  if (info.ne.0) then
    write(*,*)
    write(*,'("Error(writeevbse): diagonalisation failed")')
    write(*,'(" ZGEEV returned INFO = ",I8)') info
    write(*,*)
    stop
  end if
  evalbse(:)=dble(w(:))
  hmlbse(:,:)=vr(:,:)
  deallocate(vl,vr,rwork,work)
else
! Hermitian block only
  allocate(rwork(3*nmbse))
  lwork=2*nmbse
  allocate(work(lwork))
  call zheev('V','U',nmbse,hmlbse,nmbse,evalbse,work,lwork,rwork,info)
  if (info.ne.0) then
    write(*,*)
    write(*,'("Error(writeevbse): diagonalisation failed")')
    write(*,'(" ZHEEV returned INFO = ",I8)') info
    write(*,*)
    stop
  end if
  deallocate(rwork,work)
end if
! write the BSE eigenvectors and eigenvalues to file
open(50,file='EVBSE.OUT',action='WRITE',form='UNFORMATTED')
write(50) nmbse
write(50) evalbse
write(50) hmlbse
close(50)
! write the BSE eigenvalues to file
open(50,file='EIGVAL_BSE.OUT',action='WRITE',form='FORMATTED')
write(50,'(I6," : nmbse")') nmbse
if (bsefull) then
  do a=1,nmbse
    write(50,'(I6,2G18.10)') a,dble(w(a)),aimag(w(a))
  end do
  deallocate(w)
else
  do a=1,nmbse
    write(50,'(I6,G18.10)') a,evalbse(a)
  end do
end if
close(50)
write(*,*)
write(*,'("Info(writeevbse):")')
write(*,'(" BSE eigenvectors and eigenvalues written to EVBSE.OUT")')
write(*,'(" BSE eigenvalues written to EIGVAL_BSE.OUT")')
! deallocate global BSE arrays
deallocate(evalbse,hmlbse)
return
end subroutine

