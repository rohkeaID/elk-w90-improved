
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine wfcrplot
use modmain
implicit none
! local variables
integer ist,is,ia,ias,ir
character(256) fname
! initialise universal variables
call init0
! read density and potentials from file
call readstate
! generate the core wavefunctions
call gencore
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(fname,'("WFCORE_S",I2.2,"_A",I4.4,".OUT")') is,ia
    open(50,file=trim(fname),form='FORMATTED')
    do ist=1,nstsp(is)
      if (spcore(ist,is)) then
        do ir=1,nrsp(is)
          write(50,'(2G18.10)') rsp(ir,is),rwfcr(ir,1,ist,ias)
        end do
        write(50,'("     ")')
      end if
    end do
    close(50)
  end do
end do
write(*,*)
write(*,'("Info(wfcrplot):")')
write(*,'(" Core state wavefunctions written to WFCORE_Sss_Aaaaa.OUT")')
write(*,'(" for all species and atoms")')
return
end subroutine

