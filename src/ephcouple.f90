
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine ephcouple
use modmain
use modphonon
use modmpi
use modomp
use modstore
implicit none
! local variables
integer iq,ik,jk,jkq
integer ist,jst,isym,ip
integer is,ia,ias,js,jas
integer nr,nri,np,i,n,nthd
real(8) vl(3),de,x
real(8) t1,t2,t3,t4
! allocatable arrays
real(8), allocatable :: wq(:,:),gq(:,:)
complex(8), allocatable :: dynq(:,:,:),ev(:,:),a(:,:)
complex(8), allocatable :: dvmt(:,:,:),dvir(:,:)
complex(8), allocatable :: zfmt(:),gzfmt(:,:,:)
complex(8), allocatable :: ephmat(:,:,:)
! external functions
real(8) sdelta
external sdelta
! increase the angular momentum cut-off on the inner part of the muffin-tin
lmaxi_=lmaxi
lmaxi=max(lmaxi,4)
! initialise universal variables
call init0
call init1
call init2
! allocate global arrays
if (allocated(dvsmt)) deallocate(dvsmt)
allocate(dvsmt(npmtmax,natmtot))
if (allocated(dvsir)) deallocate(dvsir)
allocate(dvsir(ngtot))
! allocate local arrays
allocate(wq(nbph,nqpt),gq(nbph,nqpt))
allocate(dynq(nbph,nbph,nqpt),ev(nbph,nbph),a(nbph,nbph))
allocate(dvmt(npcmtmax,natmtot,nbph),dvir(ngtot,nbph))
allocate(zfmt(npmtmax),gzfmt(npmtmax,3,natmtot))
! read in the density and potentials from file
call readstate
! Fourier transform Kohn-Sham potential to G-space
call genvsig
! read Fermi energy from file
call readfermi
! set the speed of light >> 1 (non-relativistic approximation)
solsc=sol*100.d0
! new file extension for eigenvector files with c >> 1
filext='_EPH.OUT'
! generate the first- and second-variational eigenvectors and eigenvalues
call linengy
call genapwfr
call genlofr
call olprad
call hmlrad
call gensocfr
call genevfsv
! precise determine of the Fermi energy
t1=swidth
swidth=1.d-5
call occupy
swidth=t1
! restore the speed of light
solsc=sol
! read in the dynamical matrices
call readdyn(dynq)
! apply the acoustic sum rule
call sumrule(dynq)
! compute the gradients of the Kohn-Sham potential for the rigid-ion term
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  call rtozfmt(nr,nri,vsmt(:,ias),zfmt)
  call gradzfmt(nr,nri,rsp(:,is),zfmt,npmtmax,gzfmt(:,:,ias))
end do
! loop over phonon q-points
do iq=1,nqpt
  if (mp_mpi) write(*,'("Info(ephcouple): ",I6," of ",I6," q-points")') iq,nqpt
! diagonalise the dynamical matrix
  call dynev(dynq(:,:,iq),wq(:,iq),ev)
! generate the matrix for converting between Cartesian and phonon coordinates
  call genmcph(wq(:,iq),ev,a)
  i=0
  do is=1,nspecies
    nr=nrmt(is)
    nri=nrmti(is)
    np=npmt(is)
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      do ip=1,3
        i=i+1
! read in the Cartesian change in Kohn-Sham potential
        call readdvs(iq,is,ia,ip)
! add the rigid-ion term
        dvsmt(1:np,ias)=dvsmt(1:np,ias)-gzfmt(1:np,ip,ias)
        do jas=1,natmtot
          js=idxis(jas)
          call zfmtftoc(nrmt(js),nrmti(js),dvsmt(:,jas),dvmt(:,jas,i))
        end do
! multiply the interstitial potential with the characteristic function
        dvir(:,i)=dvsir(:)*cfunir(:)
      end do
    end do
  end do
! energy window for calculating the electron-phonon matrix elements averaged
! over the Fermi surface
  de=4.d0*swidth
! zero the phonon linewidths array
  gq(:,iq)=0.d0
  call omp_hold(nkptnr/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ephmat,jk,vl,isym,jkq) &
!$OMP PRIVATE(t1,t2,t3,t4,ist,jst,x,i) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
  do ik=1,nkptnr
! distribute among MPI processes
    if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
    allocate(ephmat(nstsv,nstsv,nbph))
! equivalent reduced k-point
    jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
! compute the electron-phonon coupling matrix elements
    call genephmat(iq,ik,de,a,dvmt,dvir,ephmat)
! k+q-vector in lattice coordinates
    vl(:)=vkl(:,ik)+vql(:,iq)
! index to k+q-vector
    call findkpt(vl,isym,jkq)
    t1=pi*wkptnr*occmax
! loop over second-variational states
    do ist=1,nstsv
      x=(evalsv(ist,jkq)-efermi)/swidth
      t2=t1*sdelta(stype,x)/swidth
      do jst=1,nstsv
! loop over phonon branches
        do i=1,nbph
          x=(evalsv(ist,jkq)-evalsv(jst,jk)-wq(i,iq))/swidth
          t3=t2*sdelta(stype,x)/swidth
          t4=dble(ephmat(ist,jst,i))**2+aimag(ephmat(ist,jst,i))**2
!$OMP CRITICAL(ephcouple_)
          gq(i,iq)=gq(i,iq)+wq(i,iq)*t3*t4
!$OMP END CRITICAL(ephcouple_)
        end do
      end do
    end do
    deallocate(ephmat)
! end loop over k-points
  end do
!$OMP END DO
!$OMP END PARALLEL
  call omp_free(nthd)
! end loop over phonon q-points
end do
! add gq from each process
if (np_mpi.gt.1) then
  n=nbph*nqpt
  call mpi_allreduce(mpi_in_place,gq,n,mpi_double_precision,mpi_sum,mpicom, &
   ierror)
end if
! restore the default file extension
filext='.OUT'
if (mp_mpi) then
! write the phonon linewidths to file
  call writegamma(gq)
! write electron-phonon coupling constants to file
  call writelambda(wq,gq)
end if
deallocate(wq,gq,dynq,ev,a)
deallocate(dvmt,dvir,zfmt,gzfmt)
! restore original input parameters
lmaxi=lmaxi_
return
end subroutine

