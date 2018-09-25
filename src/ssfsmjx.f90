
! Copyright (C) 2016 A. Jakobsson.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine ssfsmjx
use modmain
use modjx
use modrandom
use modstore
use modmpi
implicit none
! local variables
integer iq,i
integer is,ia,ias
real(8) tp(2),t1
! initialise universal variables
call init0
! store original variables
bfcmt0_(:,:,:)=bfcmt0(:,:,:)
spinsprl_=spinsprl
vqlss_(:)=vqlss(:)
! enable spin-spirals
spinsprl=.true.
! enable fixed spin direction
fsmtype=-2
! open SSFSMJX.OUT
if (mp_mpi) then
  open(71,file='SSFSMJX.OUT',form='FORMATTED')
end if
! loop over spin-spiral q-vectors
do iq=1,nqssjx
  if (mp_mpi) then
    write(*,'("Info(ssfsmjx): spin-spiral q-vector ",I6," of ",I6)') iq,nqssjx
  end if
! generate random q-vector
  do i=1,3
    vqlss(i)=2.d0*randomu()-1.d0
  end do
! generate random fixed spin moment directions for those atoms with non-zero
! external magnetic fields
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      tp(1)=thssjx(1)+(thssjx(2)-thssjx(1))*randomu()
      tp(2)=twopi*randomu()
      t1=bfcmt0_(3,ia,is)
      mommtfix(1,ia,is)=t1*sin(tp(1))*cos(tp(2))
      mommtfix(2,ia,is)=t1*sin(tp(1))*sin(tp(2))
      mommtfix(3,ia,is)=t1*cos(tp(1))
      bfcmt0(:,ia,is)=-mommtfix(:,ia,is)
    end do
  end do
! synchronise MPI processes
  call mpi_barrier(mpicom,ierror)
! run the ground-state spin-spiral calculation
  call gndstate
! write data to file
  if (mp_mpi) then
    write(71,*)
    write(71,'(3G18.10)') vqlss(:)
    do is=1,nspecies
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        t1=sqrt(sum(mommt(:,ias)**2))
        t1=-sign(t1,bfcmt0(3,ia,is))
        write(71,'(6G18.10)') tp(:),bfsmcmt(:,ias),t1
      end do
    end do
    flush(71)
  end if
end do
close(71)
if (mp_mpi) then
  write(*,*)
  write(*,'("Info(ssfsmjx):")')
  write(*,'(" Spin-spiral fixed spin moment data written to SSFSMJX.OUT")')
end if
! restore original input parameters
bfcmt0(:,:,:)=bfcmt0_(:,:,:)
spinsprl=spinsprl_
vqlss(:)=vqlss_(:)
return
end subroutine

