
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine energykncr
use modmain
implicit none
integer ist,is,ias,nr
! allocatable local arrays
real(8), allocatable :: rfmt(:,:)
! external functions
real(8) rfmtinp
external rfmtinp
! allocate local arrays
allocate(rfmt(lmmaxvr,nrmtmax))
! calculate the kinetic energy for core states
engykncr=0.d0
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
! sum of core eigenvalues
  do ist=1,nstsp(is)
    if (spcore(ist,is)) engykncr=engykncr+occcr(ist,ias)*evalcr(ist,ias)
  end do
! core density
  rfmt(:,:)=0.d0
  if (spincore) then
    rfmt(1,1:nr)=(rhocr(1:nr,ias,1)+rhocr(1:nr,ias,2))/y00
  else
    rfmt(1,1:nr)=rhocr(1:nr,ias,1)/y00
  end if
  engykncr=engykncr-rfmtinp(nr,nrmtinr(is),1,rsp(:,is),r2sp(:,is),rfmt, &
   vsmt(:,:,ias))
end do
deallocate(rfmt)
return
end subroutine

