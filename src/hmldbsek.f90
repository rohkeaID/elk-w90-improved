
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hmldbsek(ik2)
use modmain
implicit none
! arguments
integer, intent(in) :: ik2
! local variables
integer ik1,iv(3),iq,igq0,ig,jg
integer i1,i2,j1,j2,a1,a2,b1,b2
integer ist1,ist2,jst1,jst2
real(8) vl(3),vc(3),t0,t1,t2
complex(8) zsum
character(256) fname
! automatic arrays
integer idx(nstsv)
real(8) vcl(ngrf)
! allocatable arrays
real(8), allocatable :: vgqc(:,:),gqc(:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:)
complex(8), allocatable :: wfmt1(:,:,:,:,:),wfir1(:,:,:)
complex(8), allocatable :: wfmt2(:,:,:,:,:),wfir2(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:),zrhoir(:)
complex(8), allocatable :: zvv(:,:,:),zcc(:,:,:)
complex(8), allocatable :: zvc(:,:,:),zcv(:,:,:)
complex(8), allocatable :: epsinv(:,:,:)
allocate(vgqc(3,ngvec),gqc(ngvec))
allocate(ylmgq(lmmaxvr,ngvec),sfacgq(ngvec,natmtot))
allocate(wfmt1(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir1(ngtot,nspinor,nstsv))
allocate(wfmt2(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir2(ngtot,nspinor,nstsv))
allocate(zvv(ngrf,nvbse,nvbse),zcc(ngrf,ncbse,ncbse))
allocate(epsinv(ngrf,ngrf,nwrf))
if (bsefull) then
  allocate(zvc(ngrf,nvbse,ncbse))
  allocate(zcv(ngrf,ncbse,nvbse))
end if
! index to all states
do ist1=1,nstsv
  idx(ist1)=ist1
end do
! generate the wavefunctions for all states of k-point ik2
call genwfsvp(.false.,.false.,nstsv,idx,vkl(:,ik2),wfmt2,ngtot,wfir2)
! filename for inverse dielectric function
fname='EPSINV_RPA.OUT'
! begin loop over ik1
do ik1=1,nkptnr
! generate the wavefunctions for all states of k-point ik1
  call genwfsvp(.false.,.false.,nstsv,idx,vkl(:,ik1),wfmt1,ngtot,wfir1)
! determine equivalent q-vector in first Brillouin zone
  iv(:)=ivk(:,ik1)-ivk(:,ik2)
  iv(:)=modulo(iv(:),ngridk(:))
  iq=iqmap(iv(1),iv(2),iv(3))
! q-vector in lattice and Cartesian coordinates
  vl(:)=vkl(:,ik1)-vkl(:,ik2)
  vc(:)=vkc(:,ik1)-vkc(:,ik2)
! generate the G+q-vectors and related quantities
  call gengqrf(vc,igq0,vgqc,gqc,ylmgq,sfacgq)
! compute the regularised Coulomb interaction
  do ig=1,ngrf
    if (ig.eq.igq0) then
! volume of small parallelepiped around q-point (see genwiq2)
      t2=omegabz*wkptnr
! average symmetrised interaction over volume
      vcl(ig)=sqrt(fourpi*wiq2(iq)/t2)
    else
! G+q-vector is outside FBZ so use symmetrised 4 pi/(G+q)^2 interaction
      vcl(ig)=sqrt(fourpi)/gqc(ig)
    end if
  end do
! compute the <v|exp(-i(G+q).r)|v'> matrix elements
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zrhomt,zrhoir,ist1,ist2,i2)
!$OMP DO
  do i1=1,nvbse
    allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot),zrhoir(ngtot))
    ist1=istbse(i1,ik1)
    do i2=1,nvbse
      ist2=istbse(i2,ik2)
      call genzrho(.true.,.true.,wfmt2(:,:,:,:,ist2),wfir2(:,:,ist2), &
       wfmt1(:,:,:,:,ist1),wfir1(:,:,ist1),zrhomt,zrhoir)
      call zftzf(ngrf,gqc,ylmgq,ngrf,sfacgq,zrhomt,zrhoir,zvv(:,i1,i2))
    end do
    deallocate(zrhomt,zrhoir)
  end do
!$OMP END DO
!$OMP END PARALLEL
! compute the <c|exp(-i(G+q).r)|c'> matrix elements
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zrhomt,zrhoir,jst1,jst2,j2)
!$OMP DO
  do j1=1,ncbse
    allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot),zrhoir(ngtot))
    jst1=jstbse(j1,ik1)
    do j2=1,ncbse
      jst2=jstbse(j2,ik2)
      call genzrho(.true.,.true.,wfmt2(:,:,:,:,jst2),wfir2(:,:,jst2), &
       wfmt1(:,:,:,:,jst1),wfir1(:,:,jst1),zrhomt,zrhoir)
      call zftzf(ngrf,gqc,ylmgq,ngrf,sfacgq,zrhomt,zrhoir,zcc(:,j1,j2))
    end do
    deallocate(zrhomt,zrhoir)
  end do
!$OMP END DO
!$OMP END PARALLEL
! matrix elements for full BSE kernel if required
  if (bsefull) then
! compute the <v|exp(-i(G+q).r)|c'> matrix elements
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zrhomt,zrhoir,ist1,jst2,j2)
!$OMP DO
    do i1=1,nvbse
      allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot),zrhoir(ngtot))
      ist1=istbse(i1,ik1)
      do j2=1,ncbse
        jst2=jstbse(j2,ik2)
        call genzrho(.true.,.true.,wfmt2(:,:,:,:,jst2),wfir2(:,:,jst2), &
         wfmt1(:,:,:,:,ist1),wfir1(:,:,ist1),zrhomt,zrhoir)
        call zftzf(ngrf,gqc,ylmgq,ngrf,sfacgq,zrhomt,zrhoir,zvc(:,i1,j2))
      end do
      deallocate(zrhomt,zrhoir)
    end do
!$OMP END DO
!$OMP END PARALLEL
! compute the <c|exp(-i(G+q).r)|v'> matrix elements
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zrhomt,zrhoir,jst1,ist2,i2)
!$OMP DO
    do j1=1,ncbse
      allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot),zrhoir(ngtot))
      jst1=jstbse(j1,ik1)
      do i2=1,nvbse
        ist2=istbse(i2,ik2)
        call genzrho(.true.,.true.,wfmt2(:,:,:,:,ist2),wfir2(:,:,ist2), &
         wfmt1(:,:,:,:,jst1),wfir1(:,:,jst1),zrhomt,zrhoir)
        call zftzf(ngrf,gqc,ylmgq,ngrf,sfacgq,zrhomt,zrhoir,zcv(:,j1,i2))
      end do
      deallocate(zrhomt,zrhoir)
    end do
!$OMP END DO
!$OMP END PARALLEL
  end if
! get RPA inverse epsilon from file
  call getcf2pt(fname,vl,ngrf,nwrf,epsinv)
  t0=wkptnr*omega
  do i1=1,nvbse
    do j1=1,ncbse
      a1=ijkbse(i1,j1,ik1)
      do i2=1,nvbse
        do j2=1,ncbse
          a2=ijkbse(i2,j2,ik2)
          zsum=0.d0
          do ig=1,ngrf
            t1=t0*vcl(ig)
            do jg=1,ngrf
              t2=t1*vcl(jg)
              zsum=zsum+t2*epsinv(ig,jg,1)*conjg(zcc(ig,j1,j2))*zvv(jg,i1,i2)
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
              t1=t0*vcl(ig)
              do jg=1,ngrf
                t2=t1*vcl(jg)
                zsum=zsum+t2*epsinv(ig,jg,1)*conjg(zcv(ig,j1,i2))*zvc(jg,i1,j2)
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
deallocate(vgqc,gqc,ylmgq,sfacgq)
deallocate(wfmt1,wfmt2,wfir1,wfir2)
deallocate(zvv,zcc,epsinv)
if (bsefull) deallocate(zvc,zcv)
return
end subroutine
