
! Copyright (C) 2012 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genspchi0
! !INTERFACE:
subroutine genspchi0(ik,lock,scsr,vqpl,jlgqr,ylmgq,sfacgq,chi0)
! !USES:
use modmain
use modomp
! !INPUT/OUTPUT PARAMETERS:
!   ik     : k-point from non-reduced set (in,integer)
!   lock   : OpenMP locks for frequency index of chi0 (in,integer(nwrf))
!   scsr   : scissor correction (in,real)
!   vqpl   : input q-point in lattice coordinates (in,real(3))
!   jlgqr  : spherical Bessel functions evaluated on the coarse radial mesh for
!            all species and G+q-vectors (in,real(njcmax,nspecies,ngrf))
!   ylmgq  : spherical harmonics of the G+q-vectors (in,complex(lmmaxo,ngrf))
!   sfacgq : structure factors of G+q-vectors (in,complex(ngrf,natmtot))
!   chi0   : spin-dependent Kohn-Sham response function in G-space
!            (out,complex(ngrf,4,ngrf,4,nwrf))
! !DESCRIPTION:
!   Computes the spin-dependent Kohn-Sham response function:
!   \begin{align*}
!    \chi_{\alpha\beta,\alpha'\beta'}({\bf r},{\bf r}',\omega)
!    & \equiv\frac{\delta\rho_{\alpha\beta}({\bf r},\omega)}
!    {\delta v_{\alpha'\beta'}({\bf r}',\omega)} \\
!    & =\frac{1}{N_k}\sum_{i{\bf k},j{\bf k}'}(f_{i{\bf k}}-f_{j{\bf k}'})
!    \frac{\langle i{\bf k}|\hat{\rho}_{\beta\alpha}({\bf r})|j{\bf k}'\rangle
!    \langle j{\bf k}'|\hat{\rho}_{\alpha'\beta'}({\bf r}')|i{\bf k}\rangle}
!    {w+(\varepsilon_{i{\bf k}}-\varepsilon_{j{\bf k}'})+i\eta},
!   \end{align*}
!   where $\alpha$ and $\beta$ are spin-coordinates, $N_k$ is the number of
!   $k$-points, $f_{i{\bf k}}$ are the occupancies, $v$ is the Kohn-Sham
!   potential and $\hat{\rho}$ is the spin-density operator. With translational
!   symmetry in mind, we adopt the following convention for its Fourier
!   transform:
!   $$ \chi_{\alpha\beta,\alpha'\beta'}({\bf G},{\bf G}',{\bf q},\omega)=
!    \frac{1}{\Omega}\int d^3r\,d^3r'\,e^{-i({\bf G}+{\bf q})\cdot{\bf r}}
!    e^{i({\bf G}'+{\bf q})\cdot{\bf r}'}
!    \chi_{\alpha\beta,\alpha'\beta'}({\bf r},{\bf r}',\omega). $$
!   Let
!   $$ Z_{i{\bf k},j{\bf k}+{\bf q}}^{\alpha\beta}({\bf G})\equiv
!    \int d^3r\,e^{i({\bf G}+{\bf q})\cdot{\bf r}}
!    \varphi_{j{\bf k}+{\bf q},\alpha}^*({\bf r})
!    \varphi_{i{\bf k},\beta}({\bf r}) $$
!   then the response function in $G$-space can be written
!   $$ \chi_{\alpha\beta,\alpha'\beta'}({\bf G},{\bf G}',{\bf q},\omega)=
!    \frac{1}{N_k\Omega}\sum_{i{\bf k},j{\bf k}+{\bf q}}
!    (f_{i{\bf k}}-f_{j{\bf k}})
!    \frac{\left[Z_{i{\bf k},j{\bf k}+{\bf q}}^{\alpha\beta}({\bf G})\right]^*
!    Z_{i{\bf k},j{\bf k}+{\bf q}}^{\alpha'\beta'}({\bf G}')}
!    {w+(\varepsilon_{i{\bf k}}-\varepsilon_{j{\bf k}+{\bf q}})+i\eta}. $$
!
! !REVISION HISTORY:
!   Created March 2012 (SS and JKD)
!EOP
!BOC
implicit none
! local variables
integer, intent(in) :: ik
integer(8), intent(in) :: lock(nwrf)
real(8), intent(in) :: scsr,vqpl(3),jlgqr(njcmax,nspecies,ngrf)
complex(8), intent(in) :: ylmgq(lmmaxo,ngrf),sfacgq(ngrf,natmtot)
complex(8), intent(inout) :: chi0(ngrf,4,ngrf,4,nwrf)
! local variables
logical tz(4)
integer isym,jk,jkq,iw
integer nst,nstq,ist,jst,kst,lst
integer ig,jg,a,b,i,j,nthd
real(8) vkql(3),eij,t1
complex(8) z1
! automatic arrays
integer idx(nstsv),idxq(nstsv)
integer ngp(nspnfv),ngpq(nspnfv)
! allocatable arrays
integer, allocatable :: igpig(:,:),igpqig(:,:)
complex(8), allocatable :: wfmt(:,:,:,:),wfir(:,:,:)
complex(8), allocatable :: wfmtq(:,:,:,:),wfirq(:,:,:)
complex(8), allocatable :: zrhomt(:,:),zrhoir(:),zrhoig(:,:)
complex(8), allocatable :: zw(:),c(:,:,:,:)
if (.not.spinpol) then
  write(*,*)
  write(*,'("Error(genspchi0): spin-unpolarised calculation")')
  write(*,*)
  stop
end if
! k+q-vector in lattice coordinates
vkql(:)=vkl(:,ik)+vqpl(:)
! equivalent reduced k-points for k and k+q
jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
call findkpt(vkql,isym,jkq)
! count and index states at k and k+q in energy window
nst=0
do ist=1,nstsv
  if (abs(evalsv(ist,jk)-efermi).lt.emaxrf) then
    nst=nst+1
    idx(nst)=ist
  end if
end do
nstq=0
do jst=1,nstsv
  if (abs(evalsv(jst,jkq)-efermi).lt.emaxrf) then
    nstq=nstq+1
    idxq(nstq)=jst
  end if
end do
! generate the wavefunctions for all states at k and k+q in energy window
allocate(igpig(ngkmax,nspnfv))
allocate(wfmt(npcmtmax,natmtot,nspinor,nst),wfir(ngtc,nspinor,nst))
call genwfsvp(.false.,.false.,nst,idx,ngdc,igfc,vkl(:,ik),ngp,igpig,wfmt,ngtc, &
 wfir)
deallocate(igpig)
allocate(igpqig(ngkmax,nspnfv))
allocate(wfmtq(npcmtmax,natmtot,nspinor,nstq),wfirq(ngtc,nspinor,nstq))
call genwfsvp(.false.,.false.,nstq,idxq,ngdc,igfc,vkql,ngpq,igpqig,wfmtq,ngtc, &
 wfirq)
deallocate(igpqig)
call omp_hold(nst,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zrhomt,zrhoir,zrhoig,zw,c) &
!$OMP PRIVATE(jst,kst,lst,t1,eij,iw,i,j) &
!$OMP PRIVATE(a,b,tz,ig,jg,z1) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ist=1,nst
  allocate(zrhomt(npcmtmax,natmtot),zrhoir(ngtc))
  allocate(zrhoig(ngrf,4),zw(nwrf),c(ngrf,4,ngrf,4))
  kst=idx(ist)
  do jst=1,nstq
    lst=idxq(jst)
    t1=wkptnr*omega*(occsv(kst,jk)-occsv(lst,jkq))
    if (abs(t1).lt.1.d-8) cycle
    eij=evalsv(kst,jk)-evalsv(lst,jkq)
! scissor operator
    if (abs(scsr).gt.1.d-8) then
      if (eij.gt.0.d0) then
        eij=eij+scsr
      else
        eij=eij-scsr
      end if
    end if
! frequency-dependent part in response function formula for all frequencies
    do iw=1,nwrf
      zw(iw)=t1/(eij+wrf(iw))
    end do
! compute the complex density in G+q-space
    i=0
    do a=1,2
      do b=1,2
        i=i+1
! find which contributions are zero for collinear case
        tz(i)=.false.
        if (.not.ncmag) then
          if (((a.eq.1).and.(kst.gt.nstfv)).or. &
              ((a.eq.2).and.(kst.le.nstfv)).or. &
              ((b.eq.1).and.(lst.gt.nstfv)).or. &
              ((b.eq.2).and.(lst.le.nstfv))) then
            tz(i)=.true.
            cycle
          end if
        end if
        call genzrho(.true.,.false.,ngtc,wfmt(:,:,a,ist),wfir(:,a,ist), &
         wfmtq(:,:,b,jst),wfirq(:,b,jst),zrhomt,zrhoir)
        call zftzf(ngrf,jlgqr,ylmgq,ngrf,sfacgq,zrhomt,zrhoir,zrhoig(:,i))
      end do
    end do
! Hermitian part of matrix
    do j=1,4
      if (tz(j)) cycle
      do jg=1,ngrf
        z1=conjg(zrhoig(jg,j))
        do i=1,4
          if (tz(i)) cycle
          do ig=1,ngrf
            c(ig,i,jg,j)=zrhoig(ig,i)*z1
          end do
        end do
      end do
    end do
    do iw=1,nwrf
      z1=zw(iw)
      call omp_set_lock(lock(iw))
      do j=1,4
        if (tz(j)) cycle
        do jg=1,ngrf
          do i=1,4
            if (tz(i)) cycle
            call zaxpy(ngrf,z1,c(:,i,jg,j),1,chi0(:,i,jg,j,iw),1)
          end do
        end do
      end do
      call omp_unset_lock(lock(iw))
    end do
! end loop over jst
  end do
  deallocate(zrhomt,zrhoir,zrhoig,zw,c)
! end loop over ist
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
deallocate(wfmt,wfmtq,wfir,wfirq)
return
end subroutine
!EOC

