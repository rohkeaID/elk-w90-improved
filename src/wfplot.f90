
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine wfplot
use modmain
implicit none
! local variables
integer ik,ist,ispn
real(8) x,t1
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: evecfv(:,:),evecsv(:,:)
! external functions
real(8) sdelta
external sdelta
! initialise universal variables
call init0
call init1
! read the density and potentials from file
call readstate
! read Fermi energy from file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! set the occupancies
if ((task.eq.61).or.(task.eq.62).or.(task.eq.63)) then
  ik=kstlist(1,1)
  ist=kstlist(2,1)
  if ((ik.lt.1).or.(ik.gt.nkpt)) then
    write(*,*)
    write(*,'("Error(wfplot): k-point out of range : ",I8)') ik
    write(*,*)
    stop
  end if
  if ((ist.lt.1).or.(ist.gt.nstsv)) then
    write(*,*)
    write(*,'("Error(wfplot): state out of range : ",I8)') ist
    write(*,*)
    stop
  end if
! plotting a single wavefunction
  occsv(:,:)=0.d0
  occsv(ist,ik)=1.d0/wkpt(ik)
else
! plotting an STM image by setting occupancies to be a delta function at the
! Fermi energy
  t1=1.d0/swidth
  do ik=1,nkpt
! get the eigenvalues from file
    call getevalsv(filext,vkl(:,ik),evalsv(:,ik))
    do ist=1,nstsv
      x=(efermi-evalsv(ist,ik))*t1
      occsv(ist,ik)=occmax*wkpt(ik)*sdelta(stype,x)*t1
    end do
  end do
end if
! set the charge density to zero
rhomt(:,:,:)=0.d0
rhoir(:)=0.d0
! compute the charge density with the new occupancies
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
do ik=1,nkpt
! get the eigenvectors from file
  call getevecfv(filext,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
  call getevecsv(filext,vkl(:,ik),evecsv)
! find the matching coefficients
  do ispn=1,nspnfv
    call match(ngk(ispn,ik),gkc(:,ispn,ik),tpgkc(:,:,ispn,ik), &
     sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
  end do
! add to the density
  call rhomagk(ngk(:,ik),igkig(:,:,ik),wkpt(ik),occsv(:,ik),apwalm,evecfv, &
   evecsv)
end do
deallocate(apwalm,evecfv,evecsv)
! convert muffin-tin density/magnetisation to spherical harmonics
call rhomagsh
! symmetrise the density for the STM plot
if (task.eq.162) then
  call symrf(lradstp,rhomt,rhoir)
end if
! convert the density from a coarse to a fine radial mesh
call rfmtctof(rhomt)
! write the wavefunction modulus squared plot to file
select case(task)
case(61)
  open(50,file='WF1D.OUT',action='WRITE',form='FORMATTED')
  open(51,file='WFLINES.OUT',action='WRITE',form='FORMATTED')
  call plot1d(50,51,1,rhomt,rhoir)
  close(50)
  close(51)
  write(*,*)
  write(*,'("Info(wfplot):")')
  write(*,'(" 1D wavefunction modulus squared written to WF1D.OUT")')
  write(*,'(" vertex location lines written to WFLINES.OUT")')
case(62)
  open(50,file='WF2D.OUT',action='WRITE',form='FORMATTED')
  call plot2d(50,1,rhomt,rhoir)
  close(50)
  write(*,*)
  write(*,'("Info(wfplot):")')
  write(*,'(" 2D wavefunction modulus squared written to WF2D.OUT")')
case(162)
  open(50,file='STM2D.OUT',action='WRITE',form='FORMATTED')
  call plot2d(50,1,rhomt,rhoir)
  close(50)
  write(*,*)
  write(*,'("Info(wfplot):")')
  write(*,'(" 2D STM image written to STM2D.OUT")')
case(63)
  open(50,file='WF3D.OUT',action='WRITE',form='FORMATTED')
  call plot3d(50,1,rhomt,rhoir)
  close(50)
  write(*,*)
  write(*,'("Info(wfplot):")')
  write(*,'(" 3D wavefunction modulus squared written to WF3D.OUT")')
end select
if (task.ne.162) then
  write(*,'(" for k-point ",I8," and state ",I6)') kstlist(1,1),kstlist(2,1)
end if
return
end subroutine

