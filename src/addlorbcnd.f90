
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine addlorbcnd
use modmain
implicit none
! local variables
integer is,nlo,l,io
if (.not.lorbcnd) return
! add conduction local-orbitals to each species
do is=1,nspecies
  nlo=nlorb(is)
  do l=0,lmaxmat
    nlo=nlo+1
    if (nlo.gt.maxlorb) then
      write(*,*)
      write(*,'("Error(addlorbcnd): nlorb too large : ",I8)') nlo
      write(*,'(" for species ",I4)') is
      write(*,'("Adjust maxlorb in modmain and recompile code")')
      write(*,*)
      stop
    end if
    lorbl(nlo,is)=l
    lorbord(nlo,is)=lorbordc
    do io=1,lorbordc
      lorbe0(io,nlo,is)=0.15d0
      lorbdm(io,nlo,is)=io-1
      lorbve(io,nlo,is)=.true.
    end do
  end do
  nlorb(is)=nlo
end do
return
end subroutine

