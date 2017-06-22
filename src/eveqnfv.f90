
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: eveqnfv
! !INTERFACE:
subroutine eveqnfv(nmatp,ngp,igpig,vpc,vgpc,apwalm,evalfv,evecfv)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   nmatp  : order of overlap and Hamiltonian matrices (in,integer)
!   ngp    : number of G+p-vectors (in,integer)
!   igpig  : index from G+p-vectors to G-vectors (in,integer(ngkmax))
!   vpc    : p-vector in Cartesian coordinates (in,real(3))
!   vgpc   : G+p-vectors in Cartesian coordinates (in,real(3,ngkmax))
!   apwalm : APW matching coefficients
!            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
!   evalfv : first-variational eigenvalues (out,real(nstfv))
!   evecfv : first-variational eigenvectors (out,complex(nmatmax,nstfv))
! !DESCRIPTION:
!   Solves the eigenvalue equation,
!   $$ (H-\epsilon O)b=0, $$
!   for the all the first-variational states of the input $p$-point.
!
! !REVISION HISTORY:
!   Created March 2004 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: nmatp,ngp
integer, intent(in) :: igpig(ngkmax)
real(8), intent(in) :: vpc(3),vgpc(3,ngkmax)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
real(8), intent(out) :: evalfv(nstfv)
complex(8), intent(out) :: evecfv(nmatmax,nstfv)
! local variables
integer ias,i
real(8) ts0,ts1
! allocatable arrays
complex(8), allocatable :: h(:,:),o(:,:)
!-----------------------------------------------!
!     Hamiltonian and overlap matrix set up     !
!-----------------------------------------------!
call timesec(ts0)
allocate(h(nmatp,nmatp),o(nmatp,nmatp))
!$OMP PARALLEL SECTIONS DEFAULT(SHARED) PRIVATE(ias,i)
!$OMP SECTION
! Hamiltonian
do i=1,nmatp
  h(1:i,i)=0.d0
end do
do ias=1,natmtot
  call hmlaa(ias,ngp,apwalm,nmatp,h)
  call hmlalo(ias,ngp,apwalm,nmatp,h)
  call hmllolo(ias,ngp,nmatp,h)
end do
call hmlistl(ngp,igpig,vgpc,nmatp,h)
!$OMP SECTION
! overlap
do i=1,nmatp
  o(1:i,i)=0.d0
end do
do ias=1,natmtot
  call olpaa(ias,ngp,apwalm,nmatp,o)
  call olpalo(ias,ngp,apwalm,nmatp,o)
  call olplolo(ias,ngp,nmatp,o)
end do
call olpistl(ngp,igpig,nmatp,o)
!$OMP END PARALLEL SECTIONS
call timesec(ts1)
!$OMP CRITICAL
timemat=timemat+ts1-ts0
!$OMP END CRITICAL
!---------------------------------------!
!     solve the eigenvalue equation     !
!---------------------------------------!
if (tefvr) then
! system has inversion symmetry: use real symmetric matrix eigen solver
  call eveqnfvr(nmatp,ngp,vpc,h,o,evalfv,evecfv)
else
! no inversion symmetry: use complex Hermitian matrix eigen solver
  call eveqnfvz(nmatp,h,o,evalfv,evecfv)
end if
deallocate(h,o)
return
end subroutine
!EOC

