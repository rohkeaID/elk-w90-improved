
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine deveqnfv(n,nq,igpig,igpqig,vgpc,vgpqc,evalfv,apwalm, &
 apwalmq,dapwalm,dapwalmq,evecfv,devalfv,devecfv)
use modmain
use modphonon
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
complex(8), intent(out) :: devalfv(nstfv)
complex(8), intent(out) :: devecfv(nmatmax,nstfv)
! local variables
integer nm,nmq,ias,ist,i
integer lwork,info
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
!$OMP PARALLEL SECTIONS DEFAULT(SHARED) PRIVATE(ias)
!$OMP SECTION
! Hamiltonian
h(:,:)=0.d0
do ias=1,natmtot
  call hmlaa(ias,nq,apwalmq,nmq,h)
  call hmlalo(ias,nq,apwalmq,nmq,h)
  call hmllolo(ias,nq,nmq,h)
end do
call hmlistl(nq,igpqig,vgpqc,nmq,h)
!$OMP SECTION
! overlap
o(:,:)=0.d0
do ias=1,natmtot
  call olpaa(ias,nq,apwalmq,nmq,o)
  call olpalo(ias,nq,apwalmq,nmq,o)
  call olplolo(ias,nq,nmq,o)
end do
call olpistl(nq,igpqig,nmq,o)
!$OMP END PARALLEL SECTIONS
! solve the generalised eigenvalue problem (H - e_i)|v_i> = 0
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
!$OMP PARALLEL SECTIONS DEFAULT(SHARED) PRIVATE(ias)
!$OMP SECTION
dh(:,:)=0.d0
do ias=1,natmtot
  call dhmlaa(ias,n,nq,apwalm,apwalmq,dapwalm,dapwalmq,nmq,dh)
  call dhmlalo(ias,n,nq,apwalm,apwalmq,dapwalm,dapwalmq,nmq,dh)
  call dhmllolo(ias,n,nq,nmq,dh)
end do
call dhmlistl(n,nq,igpig,igpqig,vgpc,vgpqc,nmq,dh)
!$OMP SECTION
od(:,:)=0.d0
do ias=1,natmtot
  call dolpaa(ias,n,nq,apwalm,apwalmq,dapwalm,dapwalmq,nmq,od)
  call dolpalo(ias,n,nq,dapwalm,dapwalmq,nmq,od)
end do
call dolpistl(n,nq,igpig,igpqig,nmq,od)
!$OMP END PARALLEL SECTIONS
allocate(x(nmq),y(nmq))
! loop over states
do ist=1,nstfv
! compute |dv_i> = V(e_i - D)^(-1)V^(*t)(dH - e_i dO)|v_i>
  z1=-evalfv(ist)
  call zgemv('N',nmq,nm,z1,od,nmq,evecfv(:,ist),1,zzero,x,1)
  call zgemv('N',nmq,nm,zone,dh,nmq,evecfv(:,ist),1,zone,x,1)
! compute the first-order change in eigenvalue
  if (iqph.eq.iq0) then
    devalfv(ist)=zdotc(nmq,evecfv(:,ist),1,x,1)
  else
    devalfv(ist)=0.d0
  end if
  call zgemv('C',nmq,nmq,zone,h,nmq,x,1,zzero,y,1)
  do i=1,nmq
    t1=evalfv(ist)-w(i)
    if (abs(t1).gt.epsph) then
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

