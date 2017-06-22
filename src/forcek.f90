
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: forcek
! !INTERFACE:
subroutine forcek(ik)
! !USES:
use modmain
! !DESCRIPTION:
!   Computes the {\bf k}-dependent contribution to the incomplete basis set
!   (IBS) force. See the calling routine {\tt force} for a full description.
!
! !REVISION HISTORY:
!   Created June 2006 (JKD)
!   Updated for spin-spiral case, May 2007 (Francesco Cricchio and JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ik
! local variables
integer ispn0,ispn1,ispn,jspn
integer n,nm,nm2,is,ias,ist,jst
integer iv(3),ig,i,j,k,l
real(8) sum,t1
complex(8) z1,z2
! allocatable arrays
integer, allocatable :: ijg(:)
real(8), allocatable :: dp(:),evalfv(:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
complex(8), allocatable :: h(:),o(:)
complex(8), allocatable :: dlh(:),dlo(:)
complex(8), allocatable :: vh(:),vo(:)
complex(8), allocatable :: ffv(:,:),y(:)
! external functions
complex(8) zdotc
external zdotc
nm2=nmatmax**2
! allocate local arrays
allocate(ijg(nm2),dp(nm2))
allocate(evalfv(nstfv,nspnfv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv,nspnfv))
allocate(evecsv(nstsv,nstsv))
allocate(h(nm2),o(nm2),dlh(nm2),dlo(nm2))
allocate(vh(nmatmax),vo(nmatmax))
allocate(ffv(nstfv,nstfv),y(nstfv))
! get the eigenvalues/vectors from file
call getevalfv(filext,vkl(:,ik),evalfv)
call getevecfv(filext,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
call getevecsv(filext,vkl(:,ik),evecsv)
! loop over first-variational spin components
do jspn=1,nspnfv
  if (spinsprl) then
    ispn0=jspn; ispn1=jspn
  else
    ispn0=1; ispn1=nspinor
  end if
  n=ngk(jspn,ik)
  nm=nmat(jspn,ik)
  do j=1,n
    k=(j-1)*nm
    do i=1,j
      k=k+1
      iv(:)=ivg(:,igkig(i,jspn,ik))-ivg(:,igkig(j,jspn,ik))
      iv(:)=modulo(iv(:)-intgv(1,:),ngridg(:))+intgv(1,:)
      ijg(k)=ivgig(iv(1),iv(2),iv(3))
      dp(k)=0.5d0*dot_product(vgkc(:,i,jspn,ik),vgkc(:,j,jspn,ik))
    end do
  end do
! find the matching coefficients
  call match(n,gkc(:,jspn,ik),tpgkc(:,:,jspn,ik),sfacgk(:,:,jspn,ik),apwalm)
! loop over species and atoms
  do ias=1,natmtot
    is=idxis(ias)
! Hamiltonian and overlap matrices
    h(:)=0.d0
    call hmlaa(ias,n,apwalm,nm,h)
    call hmlalo(ias,n,apwalm,nm,h)
    o(:)=0.d0
    call olpaa(ias,n,apwalm,nm,o)
    call olpalo(ias,n,apwalm,nm,o)
! loop over Cartesian directions
    do l=1,3
! APW-APW contribution
      do j=1,n
        k=(j-1)*nm
        do i=1,j
          k=k+1
          ig=ijg(k)
          t1=vgc(l,ig)
          z1=-ffacg(ig,is)*conjg(sfacg(ig,ias))
          z2=t1*(dp(k)*z1+h(k))
          dlh(k)=cmplx(-aimag(z2),dble(z2),8)
          z2=t1*(z1+o(k))
          dlo(k)=cmplx(-aimag(z2),dble(z2),8)
        end do
      end do
      do j=n+1,nm
        k=(j-1)*nm
! APW-local-orbital contribution
        do i=1,n
          k=k+1
          t1=vgkc(l,i,jspn,ik)
          z1=t1*h(k)
          dlh(k)=cmplx(-aimag(z1),dble(z1),8)
          z1=t1*o(k)
          dlo(k)=cmplx(-aimag(z1),dble(z1),8)
        end do
! zero the local-orbital-local-orbital contribution
        do i=n+1,j
          k=k+1
          dlh(k)=0.d0
          dlo(k)=0.d0
        end do
      end do
! compute the force matrix elements in the first-variational basis
      do jst=1,nstfv
        call zhemv('U',nm,zone,dlh,nm,evecfv(:,jst,jspn),1,zzero,vh,1)
        call zhemv('U',nm,zone,dlo,nm,evecfv(:,jst,jspn),1,zzero,vo,1)
        t1=evalfv(jst,jspn)
        do ist=1,nstfv
          z1=zdotc(nm,evecfv(:,ist,jspn),1,vh,1)
          z2=zdotc(nm,evecfv(:,ist,jspn),1,vo,1)
          ffv(ist,jst)=z1-t1*z2
        end do
      end do
! compute the force using the second-variational coefficients if required
      sum=0.d0
      if (tevecsv) then
! spin-polarised case
        do j=1,nstsv
          do ispn=ispn0,ispn1
            i=(ispn-1)*nstfv+1
            call zgemv('N',nstfv,nstfv,zone,ffv,nstfv,evecsv(i,j),1,zzero,y,1)
            z1=zdotc(nstfv,evecsv(i,j),1,y,1)
            sum=sum+occsv(j,ik)*dble(z1)
          end do
        end do
      else
! spin-unpolarised case
        do j=1,nstsv
          sum=sum+occsv(j,ik)*dble(ffv(j,j))
        end do
      end if
!$OMP CRITICAL
      forceibs(l,ias)=forceibs(l,ias)+wkpt(ik)*sum
!$OMP END CRITICAL
! end loop over Cartesian components
    end do
! end loop over atoms and species
  end do
! end loop over first-variational spins
end do
deallocate(ijg,dp,evalfv,apwalm,evecfv,evecsv)
deallocate(h,o,dlh,dlo,vh,vo,ffv,y)
return
end subroutine
!EOC

