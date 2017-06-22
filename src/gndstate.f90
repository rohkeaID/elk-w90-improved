
! Copyright (C) 2002-2013 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: gndstate
! !INTERFACE:
subroutine gndstate
! !USES:
use modmain
use modmpi
use moddftu
! !DESCRIPTION:
!   Computes the self-consistent Kohn-Sham ground-state. General information is
!   written to the file {\tt INFO.OUT}. First- and second-variational
!   eigenvalues, eigenvectors and occupancies are written to the unformatted
!   files {\tt EVALFV.OUT}, {\tt EVALSV.OUT}, {\tt EVECFV.OUT}, {\tt EVECSV.OUT}
!   and {\tt OCCSV.OUT}. The density, magnetisation, Kohn-Sham potential and
!   magnetic field are written to {\tt STATE.OUT}.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!   Added MPI, August 2010 (JKD)
!EOP
!BOC
implicit none
! local variables
logical exist
integer ik,nwork,n
real(8) dv,etp,de,timetot
! allocatable arrays
real(8), allocatable :: v(:),work(:)
! initialise global variables
call init0
call init1
! set the stop signal to .false.
tstop=.false.
! initialise OEP variables if required
if (xctype(1).lt.0) call init2
if (task.eq.0) trdstate=.false.
if (task.eq.1) trdstate=.true.
! only the MPI master process should write files
if (mp_mpi) then
! write the real and reciprocal lattice vectors to file
  call writelat
! write symmetry matrices to file
  call writesym
! output the k-point set to file
  call writekpts
! write lattice vectors and atomic positions to file
  open(50,file='GEOMETRY'//trim(filext),action='WRITE',form='FORMATTED')
  call writegeom(50)
  close(50)
! write interatomic distances to file
  open(50,file='IADIST'//trim(filext),action='WRITE',form='FORMATTED')
  call writeiad(50)
  close(50)
! open INFO.OUT file
  open(60,file='INFO'//trim(filext),action='WRITE',form='FORMATTED')
! open TOTENERGY.OUT
  open(61,file='TOTENERGY'//trim(filext),action='WRITE',form='FORMATTED')
! open FERMIDOS.OUT
  open(62,file='FERMIDOS'//trim(filext),action='WRITE',form='FORMATTED')
! open MOMENT.OUT if required
  if (spinpol) open(63,file='MOMENT'//trim(filext),action='WRITE', &
   form='FORMATTED')
! open GAP.OUT
  open(64,file='GAP'//trim(filext),action='WRITE',form='FORMATTED')
! open RMSDVS.OUT
  open(65,file='RMSDVS'//trim(filext),action='WRITE',form='FORMATTED')
! open DTOTENERGY.OUT
  open(66,file='DTOTENERGY'//trim(filext),action='WRITE',form='FORMATTED')
! open TMDFTU.OUT
  if (tmwrite) open(67,file='TMDFTU'//trim(filext),action='WRITE', &
   form='FORMATTED')
! open MOMENTM.OUT
  if (spinpol) open(68,file='MOMENTM'//trim(filext),action='WRITE', &
   form='FORMATTED')
! write out general information to INFO.OUT
  call writeinfo(60)
  write(60,*)
end if
! initialise or read the charge density and potentials from file
iscl=0
if (trdstate) then
  call readstate
  if (mp_mpi) then
    write(60,'("Potential read in from STATE.OUT")')
  end if
  if (autolinengy) call readfermi
else
  call rhoinit
  call potks
  call genvsig
  if (mp_mpi) write(60,'("Density and potential initialised from atomic data")')
end if
if (mp_mpi) call flushifc(60)
! size of mixing vector
n=lmmaxvr*nrmtmax*natmtot+ngtot
if (spinpol) n=n+ndmag*(lmmaxvr*nrcmtmax*natmtot+ngtot)
if (tvmatmt) n=n+2*((lmmaxdm*nspinor)**2)*natmtot
! allocate mixing array
allocate(v(n))
! determine the size of the mixer work array
nwork=-1
call mixerifc(mixtype,n,v,dv,nwork,v)
allocate(work(nwork))
! initialise the mixer
iscl=0
call mixpack(.true.,n,v)
call mixerifc(mixtype,n,v,dv,nwork,work)
! set last self-consistent loop flag
tlast=.false.
etp=0.d0
! begin the self-consistent loop
if (mp_mpi) then
  write(60,*)
  write(60,'("+------------------------------+")')
  write(60,'("| Self-consistent loop started |")')
  write(60,'("+------------------------------+")')
end if
do iscl=1,maxscl
  if (mp_mpi) then
    write(60,*)
    write(60,'("+--------------------+")')
    write(60,'("| Loop number : ",I4," |")') iscl
    write(60,'("+--------------------+")')
  end if
  if (iscl.ge.maxscl) then
    if (mp_mpi) then
      write(60,*)
      write(60,'("Reached self-consistent loops maximum")')
    end if
    write(*,*)
    write(*,'("Warning(gndstate): failed to reach self-consistency after ",I4,&
     &" loops")') iscl
    tlast=.true.
  end if
  if (mp_mpi) call flushifc(60)
! generate the core wavefunctions and densities
  call gencore
! find the new linearisation energies
  call linengy
! write out the linearisation energies
  if (mp_mpi) call writelinen
! generate the APW radial functions
  call genapwfr
! generate the local-orbital radial functions
  call genlofr
! compute the overlap radial integrals
  call olprad
! compute the Hamiltonian radial integrals
  call hmlrad
! generate the spin-orbit coupling radial functions
  call gensocfr
! generate the first- and second-variational eigenvectors and eigenvalues
  call genevfsv
! find the occupation numbers and Fermi energy
  call occupy
  if (autoswidth.and.mp_mpi) then
    write(60,*)
    write(60,'("New smearing width : ",G18.10)') swidth
  end if
  if (mp_mpi) then
! write the occupation numbers to file
    do ik=1,nkpt
      call putoccsv(filext,ik,occsv(:,ik))
    end do
! write eigenvalues to file
    call writeeval
! write the Fermi energy to file
    call writefermi
  end if
! synchronise MPI processes
  call mpi_barrier(mpi_comm_kpt,ierror)
! generate the density and magnetisation
  call rhomag
! DFT+U or fixed tensor moment calculation
  if ((dftu.ne.0).or.(ftmtype.ne.0)) then
! generate the muffin-tin density matrix used for computing the potential matrix
    call gendmatmt
! write the FTM tensor moments to file
    if (ftmtype.ne.0) call writeftm
! generate the DFT+U or FTM muffin-tin potential matrices
    call genvmatmt
  end if
  if (dftu.ne.0) then
    if (mp_mpi) then
! write the DFT+U matrices to file
      call writedftu
! calculate and write tensor moments to file
      if (tmwrite) then
        if (spinorb) then
          call writetm3du(67)
        else
          call writetm2du(67)
        end if
      end if
    end if
  end if
! compute the Kohn-Sham potentials and magnetic fields
  call potks
  if (mp_mpi) then
    if ((xcgrad.eq.3).and.(c_tb09.ne.0.d0)) then
      write(60,*)
      write(60,'("Tran-Blaha ''09 constant c : ",G18.10)') c_tb09
    end if
  end if
! pack interstitial and muffin-tin potential and field into one array
  call mixpack(.true.,n,v)
! mix in the old potential and field with the new
  call mixerifc(mixtype,n,v,dv,nwork,work)
! make sure every MPI process has a numerically identical potential
  if (np_mpi.gt.1) then
    call mpi_bcast(v,n,mpi_double_precision,0,mpi_comm_kpt,ierror)
  end if
! unpack potential and field
  call mixpack(.false.,n,v)
! add the fixed spin moment effective field (after mixing)
  call fsmfield
! Fourier transform Kohn-Sham potential to G-space
  call genvsig
! reduce the external magnetic fields if required
  if (reducebf.lt.1.d0) then
    bfieldc(:)=bfieldc(:)*reducebf
    bfcmt(:,:,:)=bfcmt(:,:,:)*reducebf
  end if
! compute the energy components
  call energy
  if (mp_mpi) then
! output energy components
    call writeengy(60)
    write(60,*)
    write(60,'("Density of states at Fermi energy : ",G18.10)') fermidos
    write(60,'(" (states/Hartree/unit cell)")')
    write(60,*)
    write(60,'("Estimated indirect band gap : ",G18.10)') bandgap(1)
    write(60,'(" from k-point ",I6," to k-point ",I6)') ikgap(1),ikgap(2)
    write(60,'("Estimated direct band gap   : ",G18.10)') bandgap(2)
    write(60,'(" at k-point ",I6)') ikgap(3)
! write total energy to TOTENERGY.OUT
    write(61,'(G22.12)') engytot
    call flushifc(61)
! write DOS at Fermi energy to FERMIDOS.OUT
    write(62,'(G18.10)') fermidos
    call flushifc(62)
! output charges and moments
    call writechg(60)
    if (spinpol) then
! write total moment to MOMENT.OUT
      write(63,'(3G18.10)') momtot(1:ndmag)
      call flushifc(63)
! write total moment magnitude to MOMENTM.OUT
      write(68,'(G18.10)') momtotm
      call flushifc(68)
    end if
! write estimated Kohn-Sham indirect band gap
    write(64,'(G22.12)') bandgap(1)
    call flushifc(64)
! output effective fields for fixed spin moment calculations
    if (fsmtype.ne.0) call writefsm(60)
! check for WRITE file
    inquire(file='WRITE',exist=exist)
    if (exist) then
      write(60,*)
      write(60,'("WRITE file exists - writing STATE.OUT")')
      call writestate
      open(50,file='WRITE')
      close(50,status='DELETE')
    end if
! write STATE.OUT file if required
    if (nwrite.ge.1) then
      if (mod(iscl,nwrite).eq.0) then
        call writestate
        write(60,*)
        write(60,'("Wrote STATE.OUT")')
      end if
    end if
  end if
! exit self-consistent loop if required
  if (tlast) goto 10
! check for convergence
  if (iscl.ge.2) then
    if (mp_mpi) then
      write(60,*)
      write(60,'("RMS change in Kohn-Sham potential (target) : ",G18.10," (",&
       &G18.10,")")') dv,epspot
      write(65,'(G18.10)') dv
      call flushifc(65)
    end if
    de=abs(engytot-etp)
    if (mp_mpi) then
      write(60,'("Absolute change in total energy (target)   : ",G18.10," (",&
       &G18.10,")")') de,epsengy
      write(66,'(G18.10)') de
      call flushifc(66)
    end if
    if ((dv.lt.epspot).and.(de.lt.epsengy)) then
      if (mp_mpi) then
        write(60,*)
        write(60,'("Convergence targets achieved")')
      end if
      tlast=.true.
    end if
  end if
  if ((xctype(1).lt.0).and.mp_mpi) then
    write(60,*)
    write(60,'("Magnitude of OEP residual : ",G18.10)') resoep
  end if
! average the current and previous total energies and store
  if (iscl.gt.1) then
    etp=0.75d0*engytot+0.25d0*etp
  else
    etp=engytot
  end if
! check for STOP file (only master process)
  if (mp_mpi) then
    inquire(file='STOP',exist=exist)
    if (exist) then
      write(60,*)
      write(60,'("STOP file exists - stopping self-consistent loop")')
      open(50,file='STOP')
      close(50,status='DELETE')
      tstop=.true.
      tlast=.true.
    end if
  end if
! broadcast tlast and tstop from master process to all other processes
  call mpi_bcast(tlast,1,mpi_logical,0,mpi_comm_kpt,ierror)
  call mpi_bcast(tstop,1,mpi_logical,0,mpi_comm_kpt,ierror)
! output the current total CPU time
  timetot=timeinit+timemat+timefv+timesv+timerho+timepot+timefor
  if (mp_mpi) then
    write(60,*)
    write(60,'("Time (CPU seconds) : ",F12.2)') timetot
  end if
! end the self-consistent loop
end do
10 continue
! synchronise MPI processes
call mpi_barrier(mpi_comm_kpt,ierror)
if (mp_mpi) then
  write(60,*)
  write(60,'("+------------------------------+")')
  write(60,'("| Self-consistent loop stopped |")')
  write(60,'("+------------------------------+")')
! write density and potentials to file only if maxscl > 1
  if (maxscl.gt.1) then
    call writestate
    write(60,*)
    write(60,'("Wrote STATE.OUT")')
  end if
end if
! compute forces if required
if (tforce) then
  call force
! output forces to INFO.OUT
  if (mp_mpi) call writeforces(60)
end if
! total time used
timetot=timeinit+timemat+timefv+timesv+timerho+timepot+timefor
! output timing information
if (mp_mpi) then
  write(60,*)
  write(60,'("Timings (CPU seconds) :")')
  write(60,'(" initialisation",T40,": ",F12.2)') timeinit
  write(60,'(" Hamiltonian and overlap matrix set up",T40,": ",F12.2)') timemat
  write(60,'(" first-variational eigenvalue equation",T40,": ",F12.2)') timefv
  if (spinpol) then
    write(60,'(" second-variational calculation",T40,": ",F12.2)') timesv
  end if
  write(60,'(" charge density calculation",T40,": ",F12.2)') timerho
  write(60,'(" potential calculation",T40,": ",F12.2)') timepot
  if (tforce) then
    write(60,'(" force calculation",T40,": ",F12.2)') timefor
  end if
  write(60,'(" total",T40,": ",F12.2)') timetot
  write(60,*)
  write(60,'("+----------------------------+")')
  write(60,'("| Elk version ",I1.1,".",I1.1,".",I2.2," stopped |")') version
  write(60,'("+----------------------------+")')
! close the INFO.OUT file
  close(60)
! close the TOTENERGY.OUT file
  close(61)
! close the FERMIDOS.OUT file
  close(62)
! close the MOMENT.OUT and MOMENTM.OUT files
  if (spinpol) then
    close(63); close(68)
  end if
! close the GAP.OUT file
  close(64)
! close the RMSDVS.OUT file
  close(65)
! close the DTOTENERGY.OUT file
  close(66)
! close TMDFTU.OUT file
  if (tmwrite) close(67)
end if
deallocate(v,work)
! synchronise MPI processes
call mpi_barrier(mpi_comm_kpt,ierror)
return
end subroutine
!EOC
