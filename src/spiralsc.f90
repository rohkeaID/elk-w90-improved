
! Copyright (C) 2012 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine spiralsc
use modmain
use modstore
implicit none
! local variables
integer nq,iq,jq
real(8) qc
! store original parameters
natoms0(:)=natoms(:)
avec0(:,:)=avec(:,:)
atposl0(:,:,:)=atposl(:,:,:)
bfcmt00(:,:,:)=bfcmt0(:,:,:)
mommtfix0(:,:,:)=mommtfix(:,:,:)
autokpt0=autokpt
ngridk0(:)=ngridk
! initialise universal variables
call init0
! initialise q-point dependent variables
call init2
! store original parameters
atposc0(:,:,:)=atposc(:,:,:)
10 continue
call sstask(80,filext)
! if nothing more to do then restore input parameters and return
if (iqss.eq.0) then
  filext='.OUT'
  natoms(:)=natoms0(:)
  avec(:,:)=avec0(:,:)
  atposl(:,:,:)=atposl0(:,:,:)
  bfcmt0(:,:,:)=bfcmt00(:,:,:)
  mommtfix(:,:,:)=mommtfix0(:,:,:)
  autokpt=autokpt0
  ngridk(:)=ngridk0(:)
  return
end if
! spiral dry run: just generate empty SS files
if (task.eq.352) goto 10
write(*,'("Info(spiralsc): working on ",A)') 'SS'//trim(filext)
! determine k-point grid size from radkpt
autokpt=.true.
! generate the spin-spiral supercell
call genscss
! initialise or read the charge density and potentials from file
if (task.eq.350) then
  trdstate=.false.
else
  trdstate=.true.
end if
! run the ground-state calculation
call gndstate
write(80,'(I6,T20," : number of unit cells in supercell")') nscss
write(80,'(G18.10,T20," : total energy per unit cell")') engytot/dble(nscss)
write(80,*)
write(80,'("q-point in lattice and Cartesian coordinates :")')
write(80,'(3G18.10)') vql(:,iqss)
write(80,'(3G18.10)') vqc(:,iqss)
qc=sqrt(vqc(1,iqss)**2+vqc(2,iqss)**2+vqc(3,iqss)**2)
write(80,'(G18.10,T20," : length of q-vector")') qc
write(80,*)
nq=nint(dble(nqptnr)*wqpt(iqss))
write(80,'(I6,T20," : number of equivalent q-points")') nq
write(80,'("Equivalent q-points in lattice and Cartesian coordinates :")')
do iq=1,nqptnr
  jq=iqmap(ivq(1,iq),ivq(2,iq),ivq(3,iq))
  if (jq.eq.iqss) then
    write(80,'(3G18.10)') vql(:,iq)
    write(80,'(3G18.10)') vqc(:,iq)
    write(80,*)
  end if
end do
close(80)
! delete the eigenvector files
call delevec
goto 10
return
end subroutine

