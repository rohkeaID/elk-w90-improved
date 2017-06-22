
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: allatoms
! !INTERFACE:
subroutine allatoms
! !USES:
use modmain
use modxcifc
! !DESCRIPTION:
!   Solves the Kohn-Sham-Dirac equations for each atom type in the solid and
!   finds the self-consistent radial wavefunctions, eigenvalues, charge
!   densities and potentials. The atomic densities can then be used to
!   initialise the crystal densities, and the atomic self-consistent potentials
!   can be appended to the muffin-tin potentials to solve for the core states.
!   Note that, irrespective of the value of {\tt xctype}, exchange-correlation
!   functional type 3 is used. See also {\tt atoms}, {\tt rhoinit},
!   {\tt gencore} and {\tt modxcifc}.
!
! !REVISION HISTORY:
!   Created September 2002 (JKD)
!   Modified for GGA, June 2007 (JKD)
!EOP
!BOC
implicit none
logical hybrid_
integer xcspin_,xcgrad_
integer is
real(8) hybridc_
character(512) xcdescr_
! allocatable arrays
real(8), allocatable :: rwf(:,:,:)
! allocate global species charge density and potential arrays
if (allocated(rhosp)) deallocate(rhosp)
allocate(rhosp(nrspmax,nspecies))
if (allocated(vrsp)) deallocate(vrsp)
allocate(vrsp(nrspmax,nspecies))
! get the exchange-correlation functional data
call getxcdata(xctsp,xcdescr_,xcspin_,xcgrad_,hybrid_,hybridc_)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(rwf)
!$OMP DO
do is=1,nspecies
  allocate(rwf(nrspmax,2,nstspmax))
  call atom(solsc,ptnucl,spzn(is),nstsp(is),nsp(:,is),lsp(:,is),ksp(:,is), &
   occsp(:,is),xctsp,xcgrad_,nrsp(is),rsp(:,is),evalsp(:,is),rhosp(:,is), &
   vrsp(:,is),rwf)
  deallocate(rwf)
end do
!$OMP END DO
!$OMP END PARALLEL
return
end subroutine
!EOC

