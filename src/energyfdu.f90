
! Copyright (C) 2008 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: energyfdu
! !INTERFACE:
subroutine energyfdu
! !USES:
use modmain
use moddftu
use modmpi
! !DESCRIPTION:
!   Calculates the energies of radial functions to be used to calculate the
!   Slater integrals. By convention those energies are chosen to be the ones at
!   the center of the band.
!
! !REVISION HISTORY:
!   Created April 2008 (F. Cricchio)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ja,ias,jas
integer nnf,i,l
logical fnd
! automatic arrays
logical done(natmmax)
real(8) vr(nrmtmax)
nnf=0
! loop over DFT+U entries
do i=1,ndftu
  is=idftu(1,i)
  l=idftu(2,i)
  done(:)=.false.
  do ia=1,natoms(is)
    if (done(ia)) cycle
    ias=idxas(ia,is)
    vr(1:nrmt(is))=vsmt(1,1:nrmt(is),ias)*y00
! find the center of the band starting from -0.5 Ha
    fdue(l,ias)=-0.5d0
    call findband(solsc,l,nrmt(is),rsp(1,is),vr,epsband,demaxbnd,fdue(l,ias), &
     fnd)
    if (.not.fnd) nnf=nnf+1
    done(ia)=.true.
! copy to equivalent atoms
    do ja=1,natoms(is)
      if ((.not.done(ja)).and.(eqatoms(ia,ja,is))) then
        jas=idxas(ja,is)
        fdue(l,jas)=fdue(l,ias)
        done(ja)=.true.
      end if
    end do
! end loops over atoms and species
  end do
end do
if ((nnf.gt.0).and.mp_mpi) then
  write(*,*)
  write(*,'("Warning(energyfdu): could not find ",I3," energies")') nnf
end if
return
end subroutine
!EOC

