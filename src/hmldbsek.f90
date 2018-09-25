
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hmldbsek(ik2)
use modmain
use modomp
implicit none
! arguments
integer, intent(in) :: ik2
! local variables
integer ik1,ist1,ist2,jst1,jst2
integer i1,i2,j1,j2,a1,a2,b1,b2
integer iv(3),iq,ig,jg,nthd
real(8) vl(3),vc(3),t0,t1,t2
complex(8) zsum
! automatic arrays
integer idx(nstsv),ngp(nspnfv)
! allocatable arrays
integer, allocatable :: igpig(:,:)
real(8), allocatable :: vgqc(:,:),gqc(:),gclgq(:),jlgqr(:,:,:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:)
complex(8), allocatable :: wfmt1(:,:,:,:),wfir1(:,:,:)
complex(8), allocatable :: wfmt2(:,:,:,:),wfir2(:,:,:)
complex(8), allocatable :: zrhomt(:,:),zrhoir(:)
complex(8), allocatable :: zvv(:,:,:),zcc(:,:,:)
complex(8), allocatable :: zvc(:,:,:),zcv(:,:,:)
complex(8), allocatable :: epsi(:,:,:)
allocate(igpig(ngkmax,nspnfv))
allocate(vgqc(3,ngrf),gqc(ngrf),gclgq(ngrf))
allocate(jlgqr(njcmax,nspecies,ngrf))
allocate(ylmgq(lmmaxo,ngrf),sfacgq(ngrf,natmtot))
allocate(wfmt1(npcmtmax,natmtot,nspinor,nstsv),wfir1(ngtc,nspinor,nstsv))
allocate(wfmt2(npcmtmax,natmtot,nspinor,nstsv),wfir2(ngtc,nspinor,nstsv))
allocate(zvv(ngrf,nvbse,nvbse),zcc(ngrf,ncbse,ncbse))
allocate(epsi(ngrf,ngrf,nwrf))
if (bsefull) then
  allocate(zvc(ngrf,nvbse,ncbse))
  allocate(zcv(ngrf,ncbse,nvbse))
end if
! index to all states
do ist1=1,nstsv
  idx(ist1)=ist1
end do
! generate the wavefunctions for all states of k-point ik2
call genwfsvp(.false.,.false.,nstsv,idx,ngdc,igfc,vkl(:,ik2),ngp,igpig,wfmt2, &
 ngtc,wfir2)
! begin loop over ik1
do ik1=1,nkptnr
! generate the wavefunctions for all states of k-point ik1
  call genwfsvp(.false.,.false.,nstsv,idx,ngdc,igfc,vkl(:,ik1),ngp,igpig, &
   wfmt1,ngtc,wfir1)
! determine equivalent q-vector in first Brillouin zone
  iv(:)=ivk(:,ik1)-ivk(:,ik2)
  iv(:)=modulo(iv(:),ngridk(:))
  iq=iqmap(iv(1),iv(2),iv(3))
! q-vector in lattice and Cartesian coordinates
  vl(:)=vkl(:,ik1)-vkl(:,ik2)
  vc(:)=vkc(:,ik1)-vkc(:,ik2)
! generate the G+q-vectors and related quantities
  call gengqrf(vc,vgqc,gqc,jlgqr,ylmgq,sfacgq)
! generate the regularised Coulomb Green's function in G+q-space
  call gengclgq(.true.,iq,ngrf,gqc,gclgq)
! symmetrise the Coulomb Green's function
  gclgq(:)=sqrt(gclgq(:))
! compute the <v|exp(-i(G+q).r)|v'> matrix elements
  call omp_hold(nvbse,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zrhomt,zrhoir,ist1,ist2,i2) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
  do i1=1,nvbse
    allocate(zrhomt(npcmtmax,natmtot),zrhoir(ngtc))
    ist1=istbse(i1,ik1)
    do i2=1,nvbse
      ist2=istbse(i2,ik2)
      call genzrho(.true.,.true.,ngtc,wfmt2(:,:,:,ist2),wfir2(:,:,ist2), &
       wfmt1(:,:,:,ist1),wfir1(:,:,ist1),zrhomt,zrhoir)
      call zftzf(ngrf,jlgqr,ylmgq,ngrf,sfacgq,zrhomt,zrhoir,zvv(:,i1,i2))
    end do
    deallocate(zrhomt,zrhoir)
  end do
!$OMP END DO
!$OMP END PARALLEL
  call omp_free(nthd)
! compute the <c|exp(-i(G+q).r)|c'> matrix elements
  call omp_hold(ncbse,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zrhomt,zrhoir,jst1,jst2,j2) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
  do j1=1,ncbse
    allocate(zrhomt(npcmtmax,natmtot),zrhoir(ngtc))
    jst1=jstbse(j1,ik1)
    do j2=1,ncbse
      jst2=jstbse(j2,ik2)
      call genzrho(.true.,.true.,ngtc,wfmt2(:,:,:,jst2),wfir2(:,:,jst2), &
       wfmt1(:,:,:,jst1),wfir1(:,:,jst1),zrhomt,zrhoir)
      call zftzf(ngrf,jlgqr,ylmgq,ngrf,sfacgq,zrhomt,zrhoir,zcc(:,j1,j2))
    end do
    deallocate(zrhomt,zrhoir)
  end do
!$OMP END DO
!$OMP END PARALLEL
  call omp_free(nthd)
! matrix elements for full BSE kernel if required
  if (bsefull) then
! compute the <v|exp(-i(G+q).r)|c'> matrix elements
    call omp_hold(nvbse,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zrhomt,zrhoir,ist1,jst2,j2) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
    do i1=1,nvbse
      allocate(zrhomt(npcmtmax,natmtot),zrhoir(ngtc))
      ist1=istbse(i1,ik1)
      do j2=1,ncbse
        jst2=jstbse(j2,ik2)
        call genzrho(.true.,.true.,ngtc,wfmt2(:,:,:,jst2),wfir2(:,:,jst2), &
         wfmt1(:,:,:,ist1),wfir1(:,:,ist1),zrhomt,zrhoir)
        call zftzf(ngrf,jlgqr,ylmgq,ngrf,sfacgq,zrhomt,zrhoir,zvc(:,i1,j2))
      end do
      deallocate(zrhomt,zrhoir)
    end do
!$OMP END DO
!$OMP END PARALLEL
    call omp_free(nthd)
! compute the <c|exp(-i(G+q).r)|v'> matrix elements
    call omp_hold(ncbse,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zrhomt,zrhoir,jst1,ist2,i2) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
    do j1=1,ncbse
      allocate(zrhomt(npcmtmax,natmtot),zrhoir(ngtc))
      jst1=jstbse(j1,ik1)
      do i2=1,nvbse
        ist2=istbse(i2,ik2)
        call genzrho(.true.,.true.,ngtc,wfmt2(:,:,:,ist2),wfir2(:,:,ist2), &
         wfmt1(:,:,:,jst1),wfir1(:,:,jst1),zrhomt,zrhoir)
        call zftzf(ngrf,jlgqr,ylmgq,ngrf,sfacgq,zrhomt,zrhoir,zcv(:,j1,i2))
      end do
      deallocate(zrhomt,zrhoir)
    end do
!$OMP END DO
!$OMP END PARALLEL
    call omp_free(nthd)
  end if
! get RPA inverse epsilon from file
  call getcfgq('EPSINV.OUT',vl,ngrf,nwrf,epsi)
  t0=wkptnr*omega
  do i1=1,nvbse
    do j1=1,ncbse
      a1=ijkbse(i1,j1,ik1)
      do i2=1,nvbse
        do j2=1,ncbse
          a2=ijkbse(i2,j2,ik2)
          zsum=0.d0
          do ig=1,ngrf
            t1=t0*gclgq(ig)
            do jg=1,ngrf
              t2=t1*gclgq(jg)
              zsum=zsum+t2*epsi(ig,jg,1)*conjg(zcc(ig,j1,j2))*zvv(jg,i1,i2)
            end do
          end do
          hmlbse(a1,a2)=hmlbse(a1,a2)-zsum
! compute off-diagonal blocks if required
          if (bsefull) then
            b1=a1+nbbse
            b2=a2+nbbse
            hmlbse(b1,b2)=hmlbse(b1,b2)+conjg(zsum)
            zsum=0.d0
            do ig=1,ngrf
              t1=t0*gclgq(ig)
              do jg=1,ngrf
                t2=t1*gclgq(jg)
                zsum=zsum+t2*epsi(ig,jg,1)*conjg(zcv(ig,j1,i2))*zvc(jg,i1,j2)
              end do
            end do
            hmlbse(a1,b2)=hmlbse(a1,b2)-zsum
            hmlbse(b1,a2)=hmlbse(b1,a2)+conjg(zsum)
          end if
! end loop over i2 and j2
        end do
      end do
! end loop over i1 and j1
    end do
  end do
! end loop over ik1
end do
deallocate(igpig,vgqc,gqc,gclgq,jlgqr)
deallocate(ylmgq,sfacgq)
deallocate(wfmt1,wfmt2,wfir1,wfir2)
deallocate(zvv,zcc,epsi)
if (bsefull) deallocate(zvc,zcv)
return
end subroutine

