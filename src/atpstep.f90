
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: atpstep
! !INTERFACE:
subroutine atpstep
! !USES:
use modmain
use modmpi
! !DESCRIPTION:
!   Makes a geometry optimisation step and updates the current atomic positions
!   according to the force on each atom. If ${\bf r}_{ij}^m$ is the position and
!   ${\bf F}_{ij}^m$ is the force acting on it for atom $j$ of species $i$ and
!   after time step $m$, then the new position is calculated by
!   $$ {\bf r}_{ij}^{m+1}={\bf r}_{ij}^m+\tau_{ij}^m\left({\bf F}_{ij}^m
!    +{\bf F}_{ij}^{m-1}\right), $$
!   where $\tau_{ij}^m$ is a parameter governing the size of the displacement.
!   If ${\bf F}_{ij}^m\cdot{\bf F}_{ij}^{m-1}>0$ then $\tau_{ij}^m$ is
!   increased, otherwise it is decreased.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias,n
real(8) t1
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! compute the dot-product between the current and previous total force
    t1=dot_product(forcetot(:,ias),forcetotp(:,ias))
! if the force is in the same direction then increase step size parameter
    if (t1.gt.0.d0) then
      tauatp(ias)=tauatp(ias)+tau0atp
    else
      tauatp(ias)=tau0atp
    end if
! make atomic position step
    atposc(:,ia,is)=atposc(:,ia,is)+tauatp(ias)*(forcetot(:,ias) &
     +forcetotp(:,ias))
  end do
end do
! each MPI process should have identical atomic positions
n=3*maxatoms*maxspecies
call mpi_bcast(atposc,n,mpi_double_precision,0,mpi_comm_kpt,ierror)
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! compute the lattice coordinates of the atomic positions
    call r3mv(ainv,atposc(:,ia,is),atposl(:,ia,is))
  end do
end do
return
end subroutine
!EOC
