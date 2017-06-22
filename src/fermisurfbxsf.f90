
! Copyright (C) 2009 F. Cricchio, F. Bultmark and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine fermisurfbxsf
use modmain
implicit none
! local variables
integer ik,nf,f
integer nst,ist,i
integer ist0,ist1
integer jst0,jst1
integer i1,i2,i3
integer j1,j2,j3
real(8) vc(3,4),e0,e1
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
  write(*,'("Info(fermisurfbxsf): ",I6," of ",I6," k-points")') ik,nkpt
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
! plotting box in Cartesian coordinates
do i=1,4
  vc(:,i)=bvec(:,1)*kptboxl(1,i)+bvec(:,2)*kptboxl(2,i)+bvec(:,3)*kptboxl(3,i)
end do
! number of files to plot (2 for collinear magnetism, 1 otherwise)
if (ndmag.eq.1) then
  nf=2
else
  nf=1
end if
do f=1,nf
  if (nf.eq.2) then
    if (f.eq.1) then
      open(50,file='FERMISURF_UP.bxsf',action='WRITE',form='FORMATTED')
      jst0=1; jst1=nstfv
    else
      open(50,file='FERMISURF_DN.bxsf',action='WRITE',form='FORMATTED')
      jst0=nstfv+1; jst1=2*nstfv
    end if
  else
    open(50,file='FERMISURF.bxsf',action='WRITE',form='FORMATTED')
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
  write(50,'(" BEGIN_INFO")')
  write(50,'(" # Band-XCRYSDEN-Structure-File for Fermi surface plotting")')
  write(50,'(" # created by Elk version ",I1.1,".",I1.1,".",I2.2)') version
  write(50,'(" # Launch as: xcrysden --bxsf FERMISURF(_UP/_DN).bxsf")')
  write(50,'("   Fermi Energy: ",G18.10)') 0.d0
  write(50,'(" END_INFO")')
  write(50,'(" BEGIN_BLOCK_BANDGRID_3D")')
  write(50, '(" band_energies")')
  write(50,'(" BANDGRID_3D_BANDS")')
  write(50,'(I4)') nst
  write(50,'(3I6)') ngridk(:)+1
  do i=1,4
    write(50,'(3G18.10)') vc(:,i)
  end do
  do ist=ist0,ist1
    write(50,'(" BAND: ",I4)') ist
    do i1=0,ngridk(1)
      j1=mod(i1,ngridk(1))
      do i2=0,ngridk(2)
        j2=mod(i2,ngridk(2))
        do i3=0,ngridk(3)
          j3=mod(i3,ngridk(3))
          ik=ivkik(j1,j2,j3)
          write(50,'(G18.10)') evalsv(ist,ik)-efermi
        end do
      end do
    end do
  end do
  write(50,'(" END_BANDGRID_3D")')
  write(50,'(" END_BLOCK_BANDGRID_3D")')
  close(50)
end do
write(*,*)
write(*,'("Info(fermisurfbxsf):")')
if (ndmag.eq.1) then
  write(*,'(" 3D Fermi surface data written to FERMISURF_UP.bxsf and &
   &FERMISURF_DN.bxsf")')
else
  write(*,'(" 3D Fermi surface data written to FERMISURF.bxsf")')
end if
write(*,'(" for plotting with XCrysDen (Fermi energy set to zero)")')
write(*,*)
write(*,'(" Launch as: xcrysden --bxsf FERMISURF(_UP/_DN).bxsf")')
return
end subroutine

