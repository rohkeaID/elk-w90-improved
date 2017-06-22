
! Copyright (C) 2009 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writetmdu
use modmain
use moddftu
implicit none
if (dftu.eq.0) then
  write(*,*)
  write(*,'("Error(writetmdu): dftu = 0")')
  write(*,*)
  stop
end if
! initialize universal variables
call init0
! read density matrix from file DMATMT.OUT
call readdmatmt
! read Slater integrals from FDU.OUT
call readfdu
! open TMDFTU.OUT
open(50,file='TMDFTU.OUT',action='WRITE',form='FORMATTED')
if (spinorb) then
  call writetm3du(50)
else
  call writetm2du(50)
end if
close(50)
write(*,*)
write(*,'("Info(writetmdu): Tensor moment decomposition of density matrix")')
write(*,'(" and DFT+U Hartree-Fock energy written to TMDFTU.OUT")')
return
end subroutine

