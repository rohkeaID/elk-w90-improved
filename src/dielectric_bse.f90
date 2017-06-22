
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dielectric_bse
use modmain
use modtest
implicit none
! local variables
integer a1,a2,ik1,jk1
integer i1,j1,ist1,jst1
integer iw,i,j,l
integer iostat,nmbse_
real(8) e,eji,t1,t2
complex(8) eta,zv(3),z1
character(256) fname
! allocatable arrays
real(8), allocatable :: w(:)
complex(8), allocatable :: pmat(:,:,:),sigma(:,:,:)
! initialise global variables
call init0
call init1
! read Fermi energy from a file
call readfermi
! get the eigenvalues from file
do ik1=1,nkpt
  call getevalsv(filext,vkl(:,ik1),evalsv(:,ik1))
end do
! generate the BSE state index arrays
call genidxbse
! allocate global BSE arrays
if (allocated(evalbse)) deallocate(evalbse)
allocate(evalbse(nmbse))
if (allocated(hmlbse)) deallocate(hmlbse)
allocate(hmlbse(nmbse,nmbse))
! read in the BSE eigenvectors and eigenvalues
open(50,file='EVBSE.OUT',action='READ',form='UNFORMATTED',status='OLD', &
 iostat=iostat)
if (iostat.ne.0) then
  write(*,*)
  write(*,'("Error(dielectric_bse): error opening EVBSE.OUT")')
  write(*,*)
  stop
end if
read(50) nmbse_
if (nmbse.ne.nmbse_) then
  write(*,*)
  write(*,'("Error(dielectric_bse): differing nmbse")')
  write(*,'(" current   : ",I6)') nmbse
  write(*,'(" EVBSE.OUT : ",I6)') nmbse_
  stop
end if
read(50) evalbse
read(50) hmlbse
close(50)
! allocate local arrays
allocate(w(nwplot))
allocate(pmat(nstsv,nstsv,3))
allocate(sigma(3,3,nwplot))
! generate energy grid (starting from zero)
t1=wplot(2)/dble(nwplot)
do iw=1,nwplot
  w(iw)=t1*dble(iw-1)
end do
! i divided by the complex relaxation time
eta=cmplx(0.d0,swidth,8)
sigma(:,:,:)=0.d0
do a2=1,nmbse
  e=evalbse(a2)
  zv(:)=0.d0
! loop over non-reduced k-points
  do ik1=1,nkptnr
! equivalent reduced k-point
    jk1=ivkik(ivk(1,ik1),ivk(2,ik1),ivk(3,ik1))
! read the momentum matrix elements from file
    call getpmat(.false.,vkl(:,ik1),pmat)
    do i1=1,nvbse
      ist1=istbse(i1,ik1)
      do j1=1,ncbse
        jst1=jstbse(j1,ik1)
        a1=ijkbse(i1,j1,ik1)
        eji=evalsv(jst1,jk1)-evalsv(ist1,jk1)
        z1=(e/eji)*hmlbse(a1,a2)
        zv(:)=zv(:)+z1*pmat(ist1,jst1,:)
      end do
    end do
  end do
  if (abs(e).gt.1.d-8) then
    do i=1,3
      do j=1,3
        z1=zv(i)*conjg(zv(j))/e
        sigma(i,j,:)=sigma(i,j,:)+z1/(w(:)-e+eta)+conjg(z1)/(w(:)+e+eta)
      end do
    end do
  end if
end do
z1=zi*occmax*wkptnr/omega
sigma(:,:,:)=z1*sigma(:,:,:)
! loop over tensor components
do l=1,noptcomp
  i=optcomp(1,l)
  j=optcomp(2,l)
  t1=0.d0
  if (i.eq.j) t1=1.d0
  write(fname,'("EPSILON_BSE_",2I1,".OUT")') i,j
  open(60,file=trim(fname),action='WRITE',form='FORMATTED')
  do iw=1,nwplot
    t2=t1-fourpi*aimag(sigma(i,j,iw)/(w(iw)+eta))
    write(60,'(2G18.10)') w(iw),t2
  end do
  write(60,'("     ")')
  do iw=1,nwplot
    t2=fourpi*dble(sigma(i,j,iw)/(w(iw)+eta))
    write(60,'(2G18.10)') w(iw),t2
  end do
  close(60)
end do
write(*,*)
write(*,'("Info(dielectric_bse):")')
write(*,'(" dielectric tensor written to EPSILON_BSE_ij.OUT")')
write(*,'(" for components")')
do l=1,noptcomp
  write(*,'("  i = ",I1,", j = ",I1)') optcomp(1:2,l)
end do
! write sigma to test file
call writetest(187,'BSE optical conductivity',nv=nwplot,tol=1.d-3,zva=sigma)
deallocate(w,pmat,sigma)
! deallocate global BSE arrays
deallocate(evalbse,hmlbse)
return
end subroutine

