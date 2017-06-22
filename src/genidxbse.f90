
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genidxbse
use modmain
implicit none
integer ik,jk,a,ntop
integer ist,jst,i,j,k
! allocatable arrays
integer, allocatable :: idx(:)
! check if the BSE extra valence or conduction states are in range
do i=1,nvxbse
  ist=istxbse(i)
  if ((ist.lt.1).or.(ist.gt.nstsv)) then
    write(*,*)
    write(*,'("Error(genidxbse): extra valence state out of range : ",I8)') ist
    write(*,*)
    stop
  end if
end do
do j=1,ncxbse
  jst=jstxbse(j)
  if ((jst.lt.1).or.(jst.gt.nstsv)) then
    write(*,*)
    write(*,'("Error(genidxbse): extra conduction state out of range : ",I8)') &
     jst
    write(*,*)
    stop
  end if
end do
! number of valence states for transitions
nvbse=nvbse0+nvxbse
! number of conduction states for transitions
ncbse=ncbse0+ncxbse
if ((nvbse.le.0).or.(ncbse.le.0)) then
  write(*,*)
  write(*,'("Error(genidxbse): invalid number of valence or conduction &
   &transition states : ",2I8)') nvbse,ncbse
  write(*,*)
  stop
end if
! total number of transitions
nvcbse=nvbse*ncbse
! block size in BSE matrix
nbbse=nvcbse*nkptnr
! BSE matrix size
if (bsefull) then
  nmbse=2*nbbse
else
  nmbse=nbbse
end if
allocate(idx(nstsv))
! allocate global BSE index arrays
if (allocated(istbse)) deallocate(istbse)
allocate(istbse(nvbse,nkptnr))
if (allocated(jstbse)) deallocate(jstbse)
allocate(jstbse(ncbse,nkptnr))
if (allocated(ijkbse)) deallocate(ijkbse)
allocate(ijkbse(nvbse,ncbse,nkptnr))
a=0
! loop over non-reduced k-points
do ik=1,nkptnr
! equivalent reduced k-point
  jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
! index for sorting the eigenvalues into ascending order
  call sortidx(nstsv,evalsv(:,jk),idx)
! find the topmost occupied band
  ntop=nstsv
  do ist=nstsv,1,-1
    if (evalsv(idx(ist),jk).lt.efermi) then
      ntop=ist
      exit
    end if
  end do
  if ((ntop-nvbse0+1).lt.1) then
    write(*,*)
    write(*,'("Error(genidxbse): not enough valence states, reduce nvbse")')
    write(*,*)
    stop
  end if
  if ((ntop+ncbse0).gt.nstsv) then
    write(*,*)
    write(*,'("Error(genidxbse): not enough conduction states")')
    write(*,'(" reduce ncbse or increase nempty")')
    write(*,*)
    stop
  end if
! index from BSE valence states to second-variational state numbers
  do i=1,nvbse0
    istbse(i,ik)=idx(ntop-nvbse0+i)
  end do
! index from BSE conduction states to second-variational state numbers
  do j=1,ncbse0
    jstbse(j,ik)=idx(ntop+j)
  end do
! add extra states to the list
  do i=1,nvxbse
    ist=istxbse(i)
    if (evalsv(ist,jk).gt.efermi) then
      write(*,*)
      write(*,'("Error(genidxbse): extra valence state above Fermi energy : ",&
       &I6)') ist
      write(*,'(" for k-point ",I8)') jk
      write(*,*)
      stop
    end if
    do k=1,nvbse0+i-1
      if (ist.eq.istbse(k,ik)) then
        write(*,*)
        write(*,'("Error(genidxbse): redundant extra valence state : ",I6)') ist
        write(*,'(" for k-point ",I8)') jk
        write(*,*)
        stop
      end if
    end do
    istbse(nvbse0+i,ik)=ist
  end do
  do j=1,ncxbse
    jst=jstxbse(j)
    if (evalsv(jst,jk).lt.efermi) then
      write(*,*)
      write(*,'("Error(genidxbse): extra conduction state below Fermi &
       &energy : ",I6)') jst
      write(*,'(" for k-point ",I8)') jk
      write(*,*)
      stop
    end if
    do k=1,ncbse0+j-1
      if (jst.eq.jstbse(k,ik)) then
        write(*,*)
        write(*,'("Error(genidxbse): redundant extra conduction state : ",&
         &I6)') jst
        write(*,'(" for k-point ",I8)') jk
        write(*,*)
        stop
      end if
    end do
    jstbse(ncbse0+j,ik)=jst
  end do
! index from BSE valence-conduction pair and k-point to location in BSE matrix
  do i=1,nvbse
    do j=1,ncbse
      a=a+1
      ijkbse(i,j,ik)=a
    end do
  end do
! end loop over non-reduced k-points
end do
deallocate(idx)
return
end subroutine

