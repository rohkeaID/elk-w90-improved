
! Copyright (C) 2008 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine readfdu
use modmain
use moddftu
implicit none
! local variables
integer is,i,k,l
integer is_,l_
fdu(:,:)=0.d0
! read Slater integrals from FDU.OUT
open(50,file='FDU'//trim(filext),action='READ',form='FORMATTED')
do i=1,ndftu
  is=idftu(1,i)
  l=idftu(2,i)
  read(50,*)
  read(50,*) is_,l_
  if ((is.ne.is_).or.(l.ne.l_)) then
    write(*,*)
    write(*,'("Error(readfdu): differing is or l")')
    write(*,'(" current : ",2I8)') is,l
    write(*,'(" FDU.OUT : ",2I8)') is_,l_
    write(*,*)
    stop
  end if
  do k=0,2*l,2
    read(50,*) fdu(k,i)
  end do
  read(50,*) ujdu(1,i)
  read(50,*) ujdu(2,i)
  if (inpdftu.ge.4) read(50,*) lambdadu(i)
end do
close(50)
return
end subroutine

