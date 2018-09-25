
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rdmvaryc
! !INTERFACE:
subroutine rdmvaryc
! !USES:
use modmain
use modrdm
! !DESCRIPTION:
!   Calculates new {\tt evecsv} from old by using the derivatives of the total
!   energy w.r.t. {\tt evecsv}. A single step of steepest-descent is made.
!
! !REVISION HISTORY:
!   Created 2009 (Sharma)
!EOP
!BOC
implicit none
! local variables
integer ik,ist,jst
real(8) t1
complex(8) z1
! allocatable arrays
complex(8), allocatable :: dedc(:,:,:),evecsv(:,:),x(:)
! external functions
real(8) dznrm2
complex(8) zdotc
external dznrm2,zdotc
! compute the derivative w.r.t. evecsv
allocate(dedc(nstsv,nstsv,nkpt))
call rdmdedc(dedc)
allocate(evecsv(nstsv,nstsv),x(nstsv))
do ik=1,nkpt
! get the eigenvectors from file
  call getevecsv(filext,ik,vkl(:,ik),evecsv)
! calculate new evecsv
  evecsv(:,:)=evecsv(:,:)-taurdmc*dedc(:,:,ik)
! othogonalise evecsv (Gram-Schmidt)
  do ist=1,nstsv
    x(:)=evecsv(:,ist)
    do jst=1,ist-1
      z1=zdotc(nstsv,evecsv(:,jst),1,evecsv(:,ist),1)
      x(:)=x(:)-z1*evecsv(:,jst)
    end do
    t1=dznrm2(nstsv,x,1)
    t1=1.d0/t1
    evecsv(:,ist)=t1*x(:)
  end do
! write new evecsv to file
  call putevecsv(filext,ik,evecsv)
! end loop over k-points
end do
deallocate(dedc,evecsv,x)
return
end subroutine
!EOC

