
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: findsymcrys
! !INTERFACE:
subroutine findsymcrys
! !USES:
use modmain
use modmpi
use modtest
! !DESCRIPTION:
!   Finds the complete set of symmetries which leave the crystal structure
!   (including the magnetic fields) invariant. A crystal symmetry is of the form
!   $\{\alpha_S|\alpha_R|{\bf t}\}$, where ${\bf t}$ is a translation vector,
!   $\alpha_R$ is a spatial rotation operation and $\alpha_S$ is a global spin
!   rotation. Note that the order of operations is important and defined to be
!   from right to left, i.e. translation followed by spatial rotation followed
!   by spin rotation. In the case of spin-orbit coupling $\alpha_S=\alpha_R$. In
!   order to determine the translation vectors, the entire atomic basis is
!   shifted so that the first atom in the smallest set of atoms of the same
!   species is at the origin. Then all displacement vectors between atoms in
!   this set are checked as possible symmetry translations. If the global
!   variable {\tt tshift} is set to {\tt .false.} then the shift is not
!   performed. See L. M. Sandratskii and P. G. Guletskii, {\it J. Phys. F: Met.
!   Phys.} {\bf 16}, L43 (1986) and the routine {\tt findsym}.
!
! !REVISION HISTORY:
!   Created April 2007 (JKD)
!EOP
!BOC
implicit none
! local variables
integer ia,ja,is,js
integer isym,nsym,i,n
integer lspl(48),lspn(48),ilspl
real(8) v0(3),v1(3),v2(3),t1
real(8) apl(3,maxatoms,maxspecies)
! allocatable arrays
integer, allocatable :: iea(:,:,:)
real(8), allocatable :: vtl(:,:)
! allocate local array
allocate(iea(natmmax,nspecies,48))
! allocate equivalent atom arrays
if (allocated(ieqatom)) deallocate(ieqatom)
allocate(ieqatom(natmmax,nspecies,maxsymcrys))
if (allocated(eqatoms)) deallocate(eqatoms)
allocate(eqatoms(natmmax,natmmax,nspecies))
! store position of first atom
if (natmtot.gt.0) v0(:)=atposl(:,1,1)
! find the smallest set of atoms
is=1
do js=1,nspecies
  if (natoms(js).lt.natoms(is)) is=js
end do
if ((tshift).and.(natmtot.gt.0)) then
! shift basis so that the first atom in the smallest atom set is at the origin
  v1(:)=atposl(:,1,is)
  do js=1,nspecies
    do ia=1,natoms(js)
! shift atom
      atposl(:,ia,js)=atposl(:,ia,js)-v1(:)
! map lattice coordinates back to [0,1)
      call r3frac(epslat,atposl(:,ia,js))
! determine the new Cartesian coordinates
      call r3mv(avec,atposl(:,ia,js),atposc(:,ia,js))
    end do
  end do
end if
! determine possible translation vectors from smallest set of atoms
n=max(natoms(is)*natoms(is),1)
allocate(vtl(3,n))
n=1
vtl(:,1)=0.d0
do ia=1,natoms(is)
  do ja=2,natoms(is)
! compute difference between two atom vectors
    v1(:)=atposl(:,ia,is)-atposl(:,ja,is)
! map lattice coordinates to [0,1)
    call r3frac(epslat,v1)
    do i=1,n
      t1=abs(vtl(1,i)-v1(1))+abs(vtl(2,i)-v1(2))+abs(vtl(3,i)-v1(3))
      if (t1.lt.epslat) goto 10
    end do
    n=n+1
    vtl(:,n)=v1(:)
10 continue
  end do
end do
! no translations required when symtype=0,2 (F. Cricchio)
if (symtype.ne.1) n=1
eqatoms(:,:,:)=.false.
nsymcrys=0
! loop over all possible translations
do i=1,n
! construct new array with translated positions
  do is=1,nspecies
    do ia=1,natoms(is)
      apl(:,ia,is)=atposl(:,ia,is)+vtl(:,i)
    end do
  end do
! find the symmetries for current translation
  call findsym(atposl,apl,nsym,lspl,lspn,iea)
  do isym=1,nsym
    nsymcrys=nsymcrys+1
    if (nsymcrys.gt.maxsymcrys) then
      write(*,*)
      write(*,'("Error(findsymcrys): too many crystal symmetries")')
      write(*,'(" Adjust maxsymcrys in modmain and recompile code")')
      write(*,*)
      stop
    end if
    vtlsymc(:,nsymcrys)=vtl(:,i)
    lsplsymc(nsymcrys)=lspl(isym)
    lspnsymc(nsymcrys)=lspn(isym)
    do is=1,nspecies
      do ia=1,natoms(is)
        ja=iea(ia,is,isym)
        ieqatom(ia,is,nsymcrys)=ja
        eqatoms(ia,ja,is)=.true.
        eqatoms(ja,ia,is)=.true.
      end do
    end do
  end do
end do
tsyminv=.false.
do isym=1,nsymcrys
! check if inversion symmetry is present
  i=lsplsymc(isym)
  if (all(symlat(:,:,i).eq.-symlat(:,:,1))) then
    tsyminv=.true.
! make inversion the second symmetry element (the identity is the first)
    v1(:)=vtlsymc(:,isym); vtlsymc(:,isym)=vtlsymc(:,2); vtlsymc(:,2)=v1(:)
    i=lsplsymc(isym); lsplsymc(isym)=lsplsymc(2); lsplsymc(2)=i
    i=lspnsymc(isym); lspnsymc(isym)=lspnsymc(2); lspnsymc(2)=i
    do is=1,nspecies
      do ia=1,natoms(is)
        i=ieqatom(ia,is,isym)
        ieqatom(ia,is,isym)=ieqatom(ia,is,2)
        ieqatom(ia,is,2)=i
      end do
    end do
    goto 20
  end if
end do
20 continue
! if inversion exists then shift basis so that inversion center is at origin
if (tsyminv.and.tshift) then
  v1(:)=v1(:)/2.d0
  do is=1,nspecies
    do ia=1,natoms(is)
! shift atom
      atposl(:,ia,is)=atposl(:,ia,is)+v1(:)
! map lattice coordinates back to [0,1)
      call r3frac(epslat,atposl(:,ia,is))
! map lattice coordinates to [-0.5,0.5)
      do i=1,3
        if (atposl(i,ia,is).gt.0.5d0) atposl(i,ia,is)=atposl(i,ia,is)-1.d0
      end do
! determine the new Cartesian coordinates
      call r3mv(avec,atposl(:,ia,is),atposc(:,ia,is))
    end do
  end do
! recalculate crystal symmetry translation vectors
  do isym=1,nsymcrys
    ilspl=isymlat(lsplsymc(isym))
    v2(:)=symlat(:,1,ilspl)*v1(1) &
         +symlat(:,2,ilspl)*v1(2) &
         +symlat(:,3,ilspl)*v1(3)
    vtlsymc(:,isym)=vtlsymc(:,isym)-v1(:)+v2(:)
    call r3frac(epslat,vtlsymc(:,isym))
  end do
end if
! set flag for zero translation vector
do isym=1,nsymcrys
  t1=sum(abs(vtlsymc(:,isym)))
  if (t1.lt.epslat) then
    tvzsymc(isym)=.true.
  else
    tvzsymc(isym)=.false.
  end if
end do
! check inversion does not include a translation
if (tsyminv) then
  if (.not.tvzsymc(2)) tsyminv=.false.
end if
if (natmtot.gt.0) then
  v1(:)=atposl(:,1,1)-v0(:)
  t1=abs(v1(1))+abs(v1(2))+abs(v1(3))
  if ((t1.gt.epslat).and.mp_mpi) then
    write(*,'("Info(findsymcrys): atomic basis shift (lattice) :")')
    write(*,'(3G18.10)') v1(:)
    write(*,'("See GEOMETRY.OUT for new atomic positions")')
  end if
end if
! write number of crystal symmetries to test file
call writetest(705,'number of crystal symmetries',iv=nsymcrys)
deallocate(iea,vtl)
return
end subroutine
!EOC

