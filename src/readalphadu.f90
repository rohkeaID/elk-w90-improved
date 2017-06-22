
! Copyright (C) 2009 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine readalphadu
use modmain
use moddftu
implicit none
! local variables
integer is,ia,i,l
integer is_,ia_,l_
if (.not.readadu) return
open(50,file='ALPHADU'//trim(filext),action='READ',form='FORMATTED')
do i=1,ndftu
  is=idftu(1,i)
  l=idftu(2,i)
  do ia=1,natoms(is)
    read(50,*)
    read(50,*) is_,ia_,l_
    if ((is.ne.is_).or.(ia.ne.ia_).or.(l.ne.l_)) then
      write(*,*)
      write(*,'("Error(readalphadu): differing is, ia or l")')
      write(*,'(" current     : ",3I8)') is,ia,l
      write(*,'(" ALPHADU.OUT : ",3I8)') is_,ia_,l_
      write(*,*)
      stop
    end if
    read(50,*) alphadu(ia,i)
  end do
end do
close(50)
return
end subroutine

