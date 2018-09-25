
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine deveqnfv(n,nq,igpig,igpqig,vgpc,vgpqc,evalfv,apwalm, &
 apwalmq,dapwalm,dapwalmq,evecfv,devalfvp,devecfv)
use modmain
use modphonon
use modomp
implicit none
! arguments
integer, intent(in) :: n,nq
integer, intent(in) :: igpig(ngkmax),igpqig(ngkmax)
real(8), intent(in) :: vgpc(3,ngkmax),vgpqc(3,ngkmax)
real(8), intent(in) :: evalfv(nstfv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(in) :: apwalmq(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(in) :: dapwalm(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: dapwalmq(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
real(8), intent(out) :: devalfvp(nstfv)
complex(8), intent(out) :: devecfv(nmatmax,nstfv)
! local variables
integer nm,nmq,ias,ist,i
integer lwork,info,nthd
real(8) t1
complex(8) z1
! allocatable arrays
real(8), allocatable :: w(:),rwork(:)
complex(8), allocatable :: h(:,:),o(:,:),dh(:,:),od(:,:)
complex(8), allocatable :: x(:),y(:),work(:)
! external functions
complex(8) zdotc
external zdotc
! matrix sizes for k and k+q
nm=n+nlotot
nmq=nq+nlotot
allocate(h(nmq,nmq),o(nmq,nmq))
! compute the Hamiltonian and overlap matrices at p+q
call omp_hold(2,nthd)
!$OMP PARALLEL SECTIONS DEFAULT(SHARED) &
!$OMP PRIVATE(ias) &
!$OMP NUM_THREADS(nthd)
!$OMP SECTION
! Hamiltonian
h(:,:)=0.d0
do ias=1,natmtot
  call hmlaa(ias,nq,apwalmq(:,:,:,ias),nmq,h)
  call hmlalo(ias,nq,apwalmq(:,:,:,ias),nmq,h)
  call hmllolo(ias,nq,nmq,h)
end do
call hmlistl(nq,igpqig,vgpqc,nmq,h)
!$OMP SECTION
! overlap
o(:,:)=0.d0
do ias=1,natmtot
  call olpaa(ias,nq,apwalmq(:,:,:,ias),nmq,o)
  call olpalo(ias,nq,apwalmq(:,:,:,ias),nmq,o)
  call olplolo(ias,nq,nmq,o)
end do
call olpistl(nq,igpqig,nmq,o)
!$OMP END PARALLEL SECTIONS
call omp_free(nthd)
! solve the generalised eigenvalue problem (H - e_i O)|v_i> = 0
! (note: these are also the eigenvalues/vectors of O^(-1)H )
lwork=2*nmq
allocate(w(nmq),rwork(3*nmq),work(lwork))
call zhegv(1,'V','U',nmq,h,nmq,o,nmq,w,work,lwork,rwork,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(deveqnfv): diagonalisation failed")')
  write(*,'(" ZHEGV returned INFO = ",I8)') info
  write(*,*)
  stop
end if
deallocate(rwork,o,work)
! compute the Hamiltonian and overlap matrix derivatives
allocate(dh(nmq,nm),od(nmq,nm))
call omp_hold(2,nthd)
!$OMP PARALLEL SECTIONS DEFAULT(SHARED) &
!$OMP PRIVATE(ias) &
!$OMP NUM_THREADS(nthd)
!$OMP SECTION
dh(:,:)=0.d0
do ias=1,natmtot
  call dhmlaa(ias,n,nq,apwalm(:,:,:,ias),apwalmq(:,:,:,ias),dapwalm,dapwalmq, &
   nmq,dh)
  call dhmlalo(ias,n,nq,apwalm(:,:,:,ias),apwalmq(:,:,:,ias),dapwalm,dapwalmq, &
   nmq,dh)
  call dhmllolo(ias,n,nq,nmq,dh)
end do
call dhmlistl(n,nq,igpig,igpqig,vgpc,vgpqc,nmq,dh)
!$OMP SECTION
od(:,:)=0.d0
do ias=1,natmtot
  call dolpaa(ias,n,nq,apwalm(:,:,:,ias),apwalmq(:,:,:,ias),dapwalm,dapwalmq, &
   nmq,od)
  call dolpalo(ias,n,nq,dapwalm,dapwalmq,nmq,od)
end do
call dolpistl(n,nq,igpig,igpqig,nmq,od)
!$OMP END PARALLEL SECTIONS
call omp_free(nthd)
allocate(x(nmq),y(nmq))
! loop over states
do ist=1,nstfv
! compute |dv_i> = V(e_i - D)^(-1)V^(*t)(dH - e_i dO)|v_i>
  z1=-evalfv(ist)
  call zgemv('N',nmq,nm,z1,od,nmq,evecfv(:,ist),1,zzero,x,1)
  call zgemv('N',nmq,nm,zone,dh,nmq,evecfv(:,ist),1,zone,x,1)
! compute the first-order change in eigenvalue
  if (tphq0) then
    z1=zdotc(nmq,evecfv(:,ist),1,x,1)
    devalfvp(ist)=dble(z1)
  else
    devalfvp(ist)=0.d0
  end if
  call zgemv('C',nmq,nmq,zone,h,nmq,x,1,zzero,y,1)
  do i=1,nmq
    t1=evalfv(ist)-w(i)
    if (abs(t1).gt.epsdev) then
      y(i)=y(i)/t1
    else
      y(i)=0.d0
    end if
  end do
  call zgemv('N',nmq,nmq,zone,h,nmq,y,1,zzero,devecfv(:,ist),1)
end do
deallocate(w,h,dh,od,x,y)
return
end subroutine

