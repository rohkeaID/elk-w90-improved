
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genscss
use modmain
use modstore
implicit none
! local variables
integer is,ia,na,i
real(8) vc(3),cs,sn,t1
! automatic arrays
real(8) vsc(3,nqptnr)
! find the smallest supercell which contains q-vector
call findscq(iqss,avec0,nscss,vsc)
! construct supercell atomic positions and magnetic fields
do is=1,nspecies
  na=0
  do ia=1,natoms0(is)
    do i=1,nscss
      na=na+1
      if (na.gt.maxatoms) then
        write(*,*)
        write(*,'("Error(genscss): too many atoms in supercell : ",I8)') na
        write(*,'(" for species ",I4)') is
        write(*,'("Adjust maxatoms in modmain and recompile code")')
        write(*,*)
        stop
      end if
      vc(:)=vsc(:,i)+atposc0(:,ia,is)
! new atomic position in lattice coordinates
      call r3mv(ainv,vc,atposl(:,na,is))
! rotate external B-field and fixed spin moment vector by angle q.r
      t1=dot_product(vqc(:,iqss),vc(:))
      cs=cos(t1); sn=sin(t1)
      bfcmt0(1,na,is)=cs*bfcmt00(1,ia,is)-sn*bfcmt00(2,ia,is)
      bfcmt0(2,na,is)=sn*bfcmt00(1,ia,is)+cs*bfcmt00(2,ia,is)
      bfcmt0(3,na,is)=bfcmt00(3,ia,is)
      mommtfix(1,na,is)=cs*mommtfix0(1,ia,is)-sn*mommtfix0(2,ia,is)
      mommtfix(2,na,is)=sn*mommtfix0(1,ia,is)+cs*mommtfix0(2,ia,is)
      mommtfix(3,na,is)=mommtfix0(3,ia,is)
    end do
  end do
  natoms(is)=na
end do
return
end subroutine

