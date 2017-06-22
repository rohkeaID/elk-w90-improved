
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine geomopt
use modmain
use modmpi
use modstore
implicit none
! local variables
integer istp,jstp,i
real(8) ds
! store original rmtdelta (minimum distance between muffin-tin surfaces)
rmtdelta0=rmtdelta
! make rmtdelta large enough to accommodate changes in atomic positions
rmtdelta=0.2d0
! initialise global variables (and the muffin-tin radii)
call init0
call init1
! trim vsig for |G| > gmaxvr/2
trimvg0=trimvg
trimvg=.true.
! make rmtdelta small so the muffin-tin radii are not subsequently adjusted
rmtdelta=0.01d0
! store orginal volume
omega0=omega
! atomic forces are required
tforce=.true.
if (task.eq.3) then
  trdstate=.true.
else
  trdstate=.false.
end if
! initial atomic step sizes
if (allocated(tauatp)) deallocate(tauatp)
allocate(tauatp(natmtot))
tauatp(:)=tau0atp
! initialise the previous total force on each atom
if (allocated(forcetotp)) deallocate(forcetotp)
allocate(forcetotp(3,natmtot))
forcetotp(:,:)=0.d0
! initial lattice vector step size
taulatv(:)=tau0latv
! initialise current and previous stress tensor
stressp(:)=0.d0
if (mp_mpi) then
! open TOTENERGY.OUT
  open(71,file='TOTENERGY_OPT.OUT',action='WRITE',form='FORMATTED')
! open FORCEMAX.OUT
  open(72,file='FORCEMAX.OUT',action='WRITE',form='FORMATTED')
! open GEOMETRY_OPT.OUT
  open(73,file='GEOMETRY_OPT.OUT',action='WRITE',form='FORMATTED')
! open IADIST_OPT.OUT
  open(74,file='IADIST_OPT.OUT',action='WRITE',form='FORMATTED')
! open FORCES_OPT.OUT
  open(75,file='FORCES_OPT.OUT',action='WRITE',form='FORMATTED')
! open MOMENTM_OPT.OUT
  if (spinpol) then
    open(78,file='MOMENTM_OPT.OUT',action='WRITE',form='FORMATTED')
  end if
! open STRESSMAX.OUT and STRESS_OPT.OUT if required
  if (latvopt.ne.0) then
    open(86,file='STRESSMAX.OUT',action='WRITE',form='FORMATTED')
    open(87,file='STRESS_OPT.OUT',action='WRITE',form='FORMATTED')
    open(88,file='OMEGA_OPT.OUT',action='WRITE',form='FORMATTED')
  end if
end if
if (mp_mpi) write(*,*)
do istp=1,maxlatvstp
  do jstp=1,maxatpstp
    if (mp_mpi) then
      write(*,'("Info(geomopt): Lattice and atomic position optimisation &
       &steps : ",2I6)') istp,jstp
    end if
! ground-state and forces calculation
    call gndstate
! check for stop signal
    if (tstop) goto 10
! subsequent calculations will read in the potential from STATE.OUT
    trdstate=.true.
! update the atomic positions
    call atpstep
! write total energy, forces, atomic positions, interatomic distances to file
    if (mp_mpi) then
      write(71,'(G22.12)') engytot
      call flushifc(71)
      write(72,'(G18.10)') forcemax
      call flushifc(72)
      write(73,*); write(73,*)
      write(73,'("! Lattice and atomic position optimisation steps : ",2I6)') &
       istp,jstp
      call writegeom(73)
      call flushifc(73)
      write(74,*); write(74,*)
      write(74,'("Lattice and atomic position optimisation steps : ",2I6)') &
       istp,jstp
      call writeiad(74)
      call flushifc(74)
      write(75,*); write(75,*)
      write(75,'("Lattice and atomic position optimisation steps : ",2I6)') &
       istp,jstp
      call writeforces(75)
      write(75,*)
      write(75,'("Maximum force magnitude over all atoms (target) : ",G18.10,&
       &" (",G18.10,")")') forcemax,epsforce
      call flushifc(75)
      if (spinpol) then
        write(78,'(G22.12)') momtotm
        call flushifc(78)
      end if
    end if
! broadcast forcemax from master process to all other processes
    call mpi_bcast(forcemax,1,mpi_double_precision,0,mpi_comm_kpt,ierror)
! check force convergence
    if (forcemax.le.epsforce) then
      if (mp_mpi) then
        write(75,*)
        write(75,'("Force convergence target achieved")')
      end if
      exit
    end if
    if ((jstp.eq.maxatpstp).and.mp_mpi) then
      write(*,*)
      write(*,'("Warning(geomopt): atomic position optimisation failed to &
       &converge in ",I6," steps")') maxatpstp
    end if
! store the current forces array
    forcetotp(:,:)=forcetot(:,:)
! end loop over atomic position optimisation
  end do
! exit lattice optimisation loop if required
  if (latvopt.eq.0) exit
! generate the stress tensor
  call genstress
! take average of current and previous stress tensors
  stress(:)=0.5d0*(stress(:)+stressp(:))
! check for stop signal
  if (tstop) goto 10
! update the lattice vectors
  call latvstep
! write stress tensor components and maximum magnitude to file
  if (mp_mpi) then
    write(86,'(G18.10)') stressmax
    call flushifc(86)
    write(87,*)
    write(87,'("Lattice vector optimisation step : ",I6)') istp
    write(87,'("Derivative of total energy w.r.t. strain tensors :")')
    do i=1,nstrain
      write(87,'(G18.10)') stress(i)
    end do
    call flushifc(87)
    write(88,'(G18.10)') omega
    call flushifc(88)
  end if
! check for stress convergence
  if (latvopt.eq.1) then
    ds=sum(abs(stress(:)))
  else
! stress may be non-zero because of volume constraint; check change in stress
! tensor instead
    ds=sum(abs(stress(:)-stressp(:)))
  end if
! broadcase ds from master process to all other processes
  call mpi_bcast(ds,1,mpi_double_precision,0,mpi_comm_kpt,ierror)
  if ((istp.ge.3).and.(ds.le.epsstress*tau0latv)) then
    if (mp_mpi) then
      write(87,*)
      write(87,'("Stress convergence target achieved")')
    end if
    exit
  end if
  if ((istp.eq.maxlatvstp).and.mp_mpi) then
    write(*,*)
    write(*,'("Warning(geomopt): lattice vector optimisation failed to &
     &converge in ",I6," steps")') maxlatvstp
  end if
  stressp(1:nstrain)=stress(1:nstrain)
! end loop over lattice optimisation
end do
10 continue
if (mp_mpi) then
  close(71); close(72); close(73); close(74); close(75)
  if (spinpol) close(78)
  if (latvopt.ne.0) then
    close(86); close(87); close(88)
  end if
end if
! ground-state should be run again after lattice vector optimisation
if (latvopt.ne.0) call gndstate
! restore original parameters
trimvg=trimvg0
rmtdelta=rmtdelta0
return
end subroutine

