
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: hmlbdg
! !INTERFACE:
subroutine hmlbdg(adelta,h)
! !USES:
use modmain
use modscdft
! !INPUT/OUTPUT PARAMETERS:
!   adelta : anomalous potential (in,complex(nbdg,nbdg))
!   h      : BdG Hamiltonian array (out,complex(nmbdg,nmbdg))
! !DESCRIPTION:
!   Sets up the Bogoliubov-de Gennes Hamiltonian matrix:
!   $$ H=\begin{pmatrix} h & \Delta \\
!    \Delta^{\dag} & -h^* \end{pmatrix} $$
!   where $h$ and $\Delta$ are the normal and anomalous parts, respectively.
!
! !REVISION HISTORY:
!   Created January 2012 (JKD)
!EOP
!BOC
implicit none
! arguments
complex(8), intent(in) :: adelta(nbdg,nbdg)
complex(8), intent(out) :: h(nmbdg,nmbdg)
! local variables
integer i,j,ik,ist
! zero the matrix
h(:,:)=0.d0
! diagonal part of matrix
do i=1,nbdg
  ik=idxbdg(1,i)
  ist=idxbdg(2,i)
  h(i,i)=evalsv(ist,ik)
  h(nbdg+i,nbdg+i)=-conjg(h(i,i))
end do
! off-diagonal anomalous part of matrix
do i=1,nbdg
  do j=1,nbdg
    h(i,nbdg+j)=adelta(i,j)
    h(nbdg+i,j)=conjg(adelta(j,i))
  end do
end do
return
end subroutine
!EOC

