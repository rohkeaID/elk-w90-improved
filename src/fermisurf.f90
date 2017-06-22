
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine fermisurf
use modmain
implicit none
! local variables
integer ik,jk,nf,f
integer nst,ist
integer ist0,ist1
integer jst0,jst1
integer i1,i2,i3
real(8) e0,e1,prd
! allocatable arrays
integer, allocatable :: idx(:)
real(8), allocatable :: evalfv(:,:),e(:)
complex(8), allocatable :: evecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:)
! initialise universal variables
call init0
call init1
! read density and potentials from file
call readstate
! read Fermi energy from file
call readfermi
! find the new linearisation energies
call linengy
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
! begin parallel loop over reduced k-points set
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evalfv,evecfv,evecsv)
!$OMP DO
do ik=1,nkpt
  allocate(evalfv(nstfv,nspnfv))
  allocate(evecfv(nmatmax,nstfv,nspnfv))
  allocate(evecsv(nstsv,nstsv))
!$OMP CRITICAL
  write(*,'("Info(fermisurf): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL
! solve the first- and second-variational eigenvalue equations
  call eveqn(ik,evalfv,evecfv,evecsv)
  deallocate(evalfv,evecfv,evecsv)
! end loop over reduced k-points set
end do
!$OMP END DO
!$OMP END PARALLEL
! if iterative diagonalisation is used the eigenvalues must be reordered
if (tefvit.and.(.not.spinpol)) then
  allocate(idx(nstsv),e(nstsv))
  do ik=1,nkpt
    e(:)=evalsv(:,ik)
    call sortidx(nstsv,e,idx)
    do ist=1,nstsv
      evalsv(ist,ik)=e(idx(ist))
    end do
  end do
  deallocate(idx,e)
end if
! number of files to plot (2 for collinear magnetism, 1 otherwise)
if (ndmag.eq.1) then
  nf=2
else
  nf=1
end if
do f=1,nf
  if (nf.eq.2) then
    if (f.eq.1) then
      open(50,file='FERMISURF_UP.OUT',action='WRITE',form='FORMATTED')
      jst0=1; jst1=nstfv
    else
      open(50,file='FERMISURF_DN.OUT',action='WRITE',form='FORMATTED')
      jst0=nstfv+1; jst1=2*nstfv
    end if
  else
    open(50,file='FERMISURF.OUT',action='WRITE',form='FORMATTED')
    jst0=1; jst1=nstsv
  end if
! find the range of eigenvalues which contribute to the Fermi surface (Lars)
  ist0=jst1; ist1=jst0
  do ist=jst0,jst1
    e0=minval(evalsv(ist,:)); e1=maxval(evalsv(ist,:))
! determine if the band crosses the Fermi energy
    if ((e0.lt.efermi).and.(e1.gt.efermi)) then
      ist0=min(ist0,ist); ist1=max(ist1,ist)
    end if
  end do
  nst=ist1-ist0+1
  if (task.eq.100) then
! write product of eigenstates minus the Fermi energy
    write(50,'(3I6," : grid size")') np3d(:)
    do i3=0,ngridk(3)-1
      do i2=0,ngridk(2)-1
        do i1=0,ngridk(1)-1
          ik=ivkiknr(i1,i2,i3)
          jk=ivkik(i1,i2,i3)
          prd=1.d0
          do ist=ist0,ist1
            prd=prd*(evalsv(ist,jk)-efermi)
          end do
          write(50,'(4G18.10)') vkc(:,ik),prd
        end do
      end do
    end do
  else
! write the eigenvalues minus the Fermi energy separately
    write(50,'(4I6," : grid size, number of states")') np3d(:),nst
    do i3=0,ngridk(3)-1
      do i2=0,ngridk(2)-1
        do i1=0,ngridk(1)-1
          ik=ivkiknr(i1,i2,i3)
          jk=ivkik(i1,i2,i3)
          write(50,'(3G18.10,40F14.8)') vkc(:,ik),evalsv(ist0:ist1,jk)-efermi
        end do
      end do
    end do
  end if
  close(50)
end do
write(*,*)
write(*,'("Info(fermisurf):")')
if (ndmag.eq.1) then
  write(*,'(" 3D Fermi surface data written to FERMISURF_UP.OUT and &
   &FERMISURF_DN.OUT")')
else
  write(*,'(" 3D Fermi surface data written to FERMISURF.OUT")')
end if
if (task.eq.100) then
  write(*,'(" in terms of the product of eigenvalues minus the Fermi energy")')
else
  write(*,'(" in terms of separate eigenvalues minus the Fermi energy")')
end if
return
end subroutine

