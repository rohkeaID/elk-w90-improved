
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeexpmat
use modmain
implicit none
! local variables
integer nk,ik,jk,i,j
real(8) vgqc(3),gqc
real(8) a,b
! allocatable arrays
real(8), allocatable :: jlgqr(:,:)
complex(8), allocatable :: ylmgq(:),sfacgq(:)
complex(8), allocatable :: expmt(:,:),emat(:,:)
! initialise universal variables
call init0
call init1
call init2
! read in the density and potentials from file
call readstate
! read Fermi energy from file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! generate the phase factor function exp(iq.r) in the muffin-tins
allocate(jlgqr(njcmax,nspecies))
allocate(ylmgq(lmmaxo),sfacgq(natmtot))
allocate(expmt(npcmtmax,natmtot))
ngrf=1
call gengqrf(vecqc,vgqc,gqc,jlgqr,ylmgq,sfacgq)
call genexpmt(1,jlgqr,ylmgq,1,sfacgq,expmt)
expmt(:,:)=omega*expmt(:,:)
deallocate(jlgqr,ylmgq,sfacgq)
! number of k-points to write out
if (kstlist(1,1).le.0) then
  nk=nkpt
else
  nk=nkstlist
end if
open(50,file='EXPIQR.OUT',form='FORMATTED')
write(50,*)
write(50,'("q-vector (lattice coordinates) :")')
write(50,'(3G18.10)') vecql
write(50,'("q-vector (Cartesian coordinates) :")')
write(50,'(3G18.10)') vecqc
write(50,*)
write(50,'(I8," : number of k-points")') nk
write(50,'(I6," : number of states per k-point")') nstsv
allocate(emat(nstsv,nstsv))
do jk=1,nk
  if (kstlist(1,1).le.0) then
    ik=jk
  else
    ik=kstlist(1,jk)
  end if
  if ((ik.le.0).or.(ik.gt.nkpt)) then
    write(*,*)
    write(*,'("Error(writeexpiqr): k-point out of range : ",I8)') ik
    write(*,*)
    stop
  end if
  write(50,*)
  write(50,'(" k-point (lattice coordinates) :")')
  write(50,'(3G18.10)') vkl(:,ik)
  write(50,*)
  write(50,'(" k-point (Cartesian coordinates) :")')
  write(50,'(3G18.10)') vkc(:,ik)
  call genexpmat(vkl(:,ik),expmt,emat)
  do i=1,nstsv
    write(50,*)
    write(50,'(I6," : state i; state j, <...>, |<...>|^2 below")') i
    do j=1,nstsv
      a=dble(emat(i,j))
      b=aimag(emat(i,j))
      write(50,'(I6,3G18.10)') j,a,b,a**2+b**2
    end do
  end do
! end loop over k-points
end do
close(50)
write(*,*)
write(*,'("Info(writeexpmat)")')
write(*,'(" < i,k+q | exp(iq.r) | j,k > matrix elements written to &
 &EXPIQR.OUT")')
write(*,'(" for the q-vector in vecql and all k-points in kstlist")')
deallocate(expmt,emat)
return
end subroutine

