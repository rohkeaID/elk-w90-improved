
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genscph(p,dph)
use modmain
use modphonon
use modstore
implicit none
! arguments
integer, intent(in) :: p
real(8), intent(in) :: dph
! local variables
integer is,ia,na,i
real(8) vc(3),t1
if ((p.ne.0).and.(p.ne.1)) then
  write(*,*)
  write(*,'("Error(genscph): phase (p) should be 0 or 1 : ",I8)') p
  write(*,*)
  stop
end if
! find the smallest supercell which contains the q-vector
call findscq(iqph,avec0,nscph,vscph)
! construct supercell atomic positions and magnetic fields
do is=1,nspecies
  na=0
  do ia=1,natoms0(is)
    do i=1,nscph
      na=na+1
      if (na.gt.maxatoms) then
        write(*,*)
        write(*,'("Error(genscph): too many atoms in supercell : ",I8)') na
        write(*,'(" for species ",I4)') is
        write(*,'("Adjust maxatoms in modmain and recompile code")')
        write(*,*)
        stop
      end if
      vc(:)=vscph(:,i)+atposc0(:,ia,is)
! add small periodic displacement
      if ((isph.eq.is).and.(iaph.eq.ia)) then
        t1=dot_product(vqc(:,iqph),vscph(:,i))
        if (p.eq.0) then
          vc(ipph)=vc(ipph)+dph*cos(t1)
        else
          vc(ipph)=vc(ipph)+dph*sin(t1)
        end if
      end if
! convert to new lattice coordinates
      call r3mv(ainv,vc,atposl(:,na,is))
      call r3frac(epslat,atposl(:,na,is))
! set muffin-tin fields and fixed spin moments if required
      if (spinpol) then
        bfcmt0(:,na,is)=bfcmt00(:,ia,is)
        mommtfix(:,na,is)=mommtfix0(:,ia,is)
      end if
    end do
  end do
  natoms(is)=na
end do
return
end subroutine

