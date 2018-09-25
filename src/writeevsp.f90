
! Copyright (C) 2015 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeevsp
use modmain
implicit none
! local variables
integer is,ist
! solve the atomic Dirac-Kohn-Sham ground state for all species
call init0
! write out the atomic eigenvalues for each species
open(50,file='EVALSP.OUT',form='FORMATTED')
write(50,*)
write(50,'("Kohn-Sham-Dirac eigenvalues for all atomic species")')
write(50,*)
write(50,'("Exchange-correlation functional : ",3I6)') xctsp(:)
do is=1,nspecies
  write(50,*)
  write(50,'("Species : ",I4," (",A,")",I4)') is,trim(spsymb(is))
  do ist=1,nstsp(is)
    write(50,'(" n = ",I2,", l = ",I2,", k = ",I2," : ",G18.10)') nsp(ist,is), &
     lsp(ist,is),ksp(ist,is),evalsp(ist,is)
  end do
end do
close(50)
write(*,*)
write(*,'("Info(writeevsp)")')
write(*,'(" Kohn-Sham-Dirac eigenvalues written to EVALSP.OUT for all atomic &
 &species")')
return
end subroutine

