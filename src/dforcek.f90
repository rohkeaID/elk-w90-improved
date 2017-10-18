
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dforcek(ik,dyn)
use modmain
use modphonon
use modpw
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(inout) :: dyn(3,natmtot)
! local variables
integer ispn0,ispn1,ispn,jspn
integer n,nq,nm,nmq,nm2
integer is,ias,ist,jst,jk
integer iv(3),ig,i,j,k,l
real(8) t1
complex(8) zsum,z1,z2
complex(8) dt1,dz1,dz2
! allocatable arrays
integer, allocatable :: ijg(:,:),ijgq(:,:)
real(8), allocatable :: dp(:,:),dpq(:,:),evalfv(:,:)
complex(8), allocatable :: apwalm(:,:,:,:),apwalmq(:,:,:,:),dapwalm(:,:,:)
complex(8), allocatable :: evecfv(:,:,:),devecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:),devecsv(:,:)
complex(8), allocatable :: h(:),o(:),dlh(:),dlo(:)
complex(8), allocatable :: hq(:,:),oq(:,:),dh(:,:),od(:,:)
complex(8), allocatable :: dlhq(:,:),dloq(:,:),ddlh(:,:),ddlo(:,:)
complex(8), allocatable :: vh(:),vo(:),dvh(:),dvo(:)
complex(8), allocatable :: ffv(:,:),dffv(:,:),y(:),dy(:)
! external functions
complex(8) zdotc
external zdotc
nm2=nmatmax**2
! allocate local arrays
allocate(ijg(nmatmax,nmatmax),ijgq(nmatmax,nmatmax))
allocate(dp(nmatmax,nmatmax),dpq(nmatmax,nmatmax))
allocate(evalfv(nstfv,nspnfv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(apwalmq(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(dapwalm(ngkmax,apwordmax,lmmaxapw))
allocate(evecfv(nmatmax,nstfv,nspnfv))
allocate(devecfv(nmatmax,nstfv,nspnfv))
allocate(h(nm2),o(nm2),dlh(nm2),dlo(nm2))
allocate(hq(nmatmax,nmatmax),oq(nmatmax,nmatmax))
allocate(dh(nmatmax,nmatmax),od(nmatmax,nmatmax))
allocate(dlhq(nmatmax,nmatmax),dloq(nmatmax,nmatmax))
allocate(ddlh(nmatmax,nmatmax),ddlo(nmatmax,nmatmax))
allocate(vh(nmatmax),vo(nmatmax),dvh(nmatmax),dvo(nmatmax))
allocate(ffv(nstfv,nstfv),dffv(nstfv,nstfv),y(nstfv),dy(nstfv))
! equivalent reduced k-point
jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
! get the eigenvalues/vectors from file
call getevalfv(filext,0,vkl(:,ik),evalfv)
call getevecfv(filext,0,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
! get the eigenvalue/vector derivatives from file
call getdevecfv(ik,iqph,isph,iaph,ipph,devecfv)
if (tevecsv) then
  allocate(evecsv(nstsv,nstsv),devecsv(nstsv,nstsv))
  call getevecsv(filext,0,vkl(:,ik),evecsv)
  call getdevecsv(ik,iqph,isph,iaph,ipph,devecsv)
end if
! loop over first-variational spin components
do jspn=1,nspnfv
  if (spinsprl) then
    ispn0=jspn; ispn1=jspn
  else
    ispn0=1; ispn1=nspinor
  end if
  n=ngk(jspn,ik)
  nq=ngkq(jspn,ik)
  nm=n+nlotot
  nmq=nq+nlotot
  do j=1,n
    do i=1,n
      iv(:)=ivg(:,igkig(i,jspn,ik))-ivg(:,igkig(j,jspn,ik))
      iv(:)=modulo(iv(:)-intgv(1,:),ngridg(:))+intgv(1,:)
      ijg(i,j)=ivgig(iv(1),iv(2),iv(3))
      dp(i,j)=0.5d0*dot_product(vgkc(:,i,jspn,ik),vgkc(:,j,jspn,ik))
    end do
  end do
  do j=1,n
    do i=1,nq
      iv(:)=ivg(:,igkqig(i,jspn,ik))-ivg(:,igkig(j,jspn,ik))
      iv(:)=modulo(iv(:)-intgv(1,:),ngridg(:))+intgv(1,:)
      ijgq(i,j)=ivgig(iv(1),iv(2),iv(3))
      dpq(i,j)=0.5d0*dot_product(vgkqc(:,i,jspn,ik),vgkc(:,j,jspn,ik))
    end do
  end do
! find the matching coefficients
  call match(n,gkc(:,jspn,ik),tpgkc(:,:,jspn,ik),sfacgk(:,:,jspn,ik),apwalm)
  call match(nq,gkqc(:,jspn,ik),tpgkqc(:,:,jspn,ik),sfacgkq(:,:,jspn,ik), &
   apwalmq)
! find the matching coefficient derivatives
  call dmatch(iasph,ipph,n,vgkc(:,:,jspn,ik),apwalm,dapwalm)
! loop over species and atoms
  do ias=1,natmtot
    is=idxis(ias)
! Hamiltonian and overlap matrices
    h(:)=0.d0
    call hmlaa(ias,n,apwalm(:,:,:,ias),nm,h)
    call hmlalo(ias,n,apwalm(:,:,:,ias),nm,h)
    o(:)=0.d0
    call olpaa(ias,n,apwalm(:,:,:,ias),nm,o)
    call olpalo(ias,n,apwalm(:,:,:,ias),nm,o)
    hq(:,:)=0.d0
    call hmlaaq(ias,n,nq,apwalm(:,:,:,ias),apwalmq(:,:,:,ias),nmatmax,hq)
    call hmlaloq(ias,n,nq,apwalm(:,:,:,ias),apwalmq(:,:,:,ias),nmatmax,hq)
    oq(:,:)=0.d0
    call olpaaq(ias,n,nq,apwalm(:,:,:,ias),apwalmq(:,:,:,ias),nmatmax,oq)
    call olpaloq(ias,n,nq,apwalm(:,:,:,ias),apwalmq(:,:,:,ias),nmatmax,oq)
! Hamiltonian and overlap derivatives
    dh(:,:)=0.d0
    call dhmlaa(ias,n,n,apwalm(:,:,:,ias),apwalm(:,:,:,ias),dapwalm,dapwalm, &
     nmatmax,dh)
    call dhmlalo(ias,n,n,apwalm(:,:,:,ias),apwalm(:,:,:,ias),dapwalm,dapwalm, &
     nmatmax,dh)
    od(:,:)=0.d0
    call dolpaa(ias,n,n,apwalm(:,:,:,ias),apwalm(:,:,:,ias),dapwalm,dapwalm, &
     nmatmax,od)
    call dolpalo(ias,n,n,dapwalm,dapwalm,nmatmax,od)
! loop over Cartesian directions
    do l=1,3
! APW-APW contribution
      do j=1,n
        k=(j-1)*nm
        do i=1,j
          k=k+1
          ig=ijg(i,j)
          t1=vgc(l,ig)
          z1=-ffacg(ig,is)*conjg(sfacg(ig,ias))
          z2=t1*(dp(i,j)*z1+h(k))
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
! non-square H/O(G+k+q,G'+k) matrices
! APW-APW contribution
      do j=1,n
        do i=1,nq
          ig=ijgq(i,j)
          t1=vgqc(l,ig)
          z1=-ffacgq(ig,is)*conjg(sfacgq(ig,ias))
          z2=t1*(dpq(i,j)*z1+hq(i,j))
          dlhq(i,j)=cmplx(-aimag(z2),dble(z2),8)
          z2=t1*(z1+oq(i,j))
          dloq(i,j)=cmplx(-aimag(z2),dble(z2),8)
        end do
! local-orbital-APW derivative
        t1=-vgkc(l,j,jspn,ik)
        do i=nq+1,nmq
          z1=t1*hq(i,j)
          dlhq(i,j)=cmplx(-aimag(z1),dble(z1),8)
          z1=t1*oq(i,j)
          dloq(i,j)=cmplx(-aimag(z1),dble(z1),8)
        end do
      end do
      do j=n+1,nm
! APW-local-orbital contribution
        do i=1,nq
          t1=vgkqc(l,i,jspn,ik)
          z1=t1*hq(i,j)
          dlhq(i,j)=cmplx(-aimag(z1),dble(z1),8)
          z1=t1*oq(i,j)
          dloq(i,j)=cmplx(-aimag(z1),dble(z1),8)
        end do
! zero the local-orbital-local-orbital contribution
        do i=nq+1,nmq
          dlhq(i,j)=0.d0
          dloq(i,j)=0.d0
        end do
      end do
! APW-APW derivative
      do j=1,n
        do i=1,n
          ig=ijg(i,j)
          t1=vgc(l,ig)
          if (ias.eq.iasph) then
            z1=-ffacg(ig,is)*conjg(sfacg(ig,ias))
            dz1=vgc(ipph,ig)*cmplx(aimag(z1),-dble(z1),8)
          else
            dz1=0.d0
          end if
          z2=t1*(dp(i,j)*dz1+dh(i,j))
          ddlh(i,j)=cmplx(-aimag(z2),dble(z2),8)
          z2=t1*(dz1+od(i,j))
          ddlo(i,j)=cmplx(-aimag(z2),dble(z2),8)
        end do
! local-orbital-APW derivative
        t1=-vgkc(l,j,jspn,ik)
        do i=n+1,nm
          z1=t1*dh(i,j)
          ddlh(i,j)=cmplx(-aimag(z1),dble(z1),8)
          z1=t1*od(i,j)
          ddlo(i,j)=cmplx(-aimag(z1),dble(z1),8)
        end do
      end do
! APW-local-orbital derivative
      do j=n+1,nm
        do i=1,n
          t1=vgkc(l,i,jspn,ik)
          z1=t1*dh(i,j)
          ddlh(i,j)=cmplx(-aimag(z1),dble(z1),8)
          z1=t1*od(i,j)
          ddlo(i,j)=cmplx(-aimag(z1),dble(z1),8)
        end do
! zero the local-orbital-local-orbital derivative
        do i=n+1,nm
          ddlh(i,j)=0.d0
          ddlo(i,j)=0.d0
        end do
      end do
      if (tphq0) then
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
      end if
! compute the force derivative matrix elements in the first-variational basis
      dffv(:,:)=0.d0
      do jst=1,nstfv
        call zhemv('U',nm,zone,dlo,nm,evecfv(:,jst,jspn),1,zzero,vo,1)
        call zgemv('N',nm,nm,zone,ddlh,nmatmax,evecfv(:,jst,jspn),1,zzero,dvh,1)
        call zgemv('N',nm,nm,zone,ddlo,nmatmax,evecfv(:,jst,jspn),1,zzero,dvo,1)
        t1=evalfv(jst,jspn)
        dt1=devalfv(jst,jspn,ik)
        do ist=1,nstfv
          z2=zdotc(nm,evecfv(:,ist,jspn),1,vo,1)
          dz1=zdotc(nm,evecfv(:,ist,jspn),1,dvh,1)
          dz2=zdotc(nm,evecfv(:,ist,jspn),1,dvo,1)
          dffv(ist,jst)=dffv(ist,jst)+dz1-dt1*z2-t1*dz2
        end do
        call zgemv('C',nmq,nm,zone,dlhq,nmatmax,devecfv(:,jst,jspn),1,zzero, &
         dvh,1)
        call zgemv('C',nmq,nm,zone,dloq,nmatmax,devecfv(:,jst,jspn),1,zzero, &
         dvo,1)
        do ist=1,nstfv
          dz1=2.d0*zdotc(nm,evecfv(:,ist,jspn),1,dvh,1)
          dz2=2.d0*zdotc(nm,evecfv(:,ist,jspn),1,dvo,1)
          dffv(ist,jst)=dffv(ist,jst)+dz1-t1*dz2
        end do
      end do
      zsum=0.d0
      if (tevecsv) then
! spin-polarised case
        do j=1,nstsv
          do ispn=ispn0,ispn1
            i=(ispn-1)*nstfv+1
            call zgemv('N',nstfv,nstfv,zone,ffv,nstfv,evecsv(i,j),1,zzero,y,1)
            call zgemv('N',nstfv,nstfv,zone,dffv,nstfv,evecsv(i,j),1,zzero,dy,1)
            call zgemv('N',nstfv,nstfv,zone,ffv,nstfv,devecsv(i,j),1,zone,dy,1)
            dz1=zdotc(nstfv,evecsv(i,j),1,dy,1)
            dz1=dz1+zdotc(nstfv,devecsv(i,j),1,y,1)
            zsum=zsum+occsv(j,jk)*dz1
!******** doccsv
          end do
        end do
      else
! spin-unpolarised case
        do j=1,nstsv
          zsum=zsum+occsv(j,jk)*dffv(j,j)
          if (tphq0) then
            zsum=zsum+doccsv(j,ik)*dble(ffv(j,j))
          end if
        end do
      end if
!$OMP CRITICAL
      dyn(l,ias)=dyn(l,ias)-wkptnr*zsum
!$OMP END CRITICAL
! end loop over Cartesian components
    end do
! end loop over atoms and species
  end do
! end loop over first-variational spins
end do
deallocate(ijg,ijgq,dp,dpq,evalfv)
deallocate(apwalm,apwalmq,dapwalm)
deallocate(evecfv,devecfv)
if (tevecsv) deallocate(evecsv,devecsv)
deallocate(h,o,dlh,dlo,hq,oq,dh,od)
deallocate(dlhq,dloq,ddlh,ddlo)
deallocate(vh,vo,dvh,dvo,ffv,dffv,y,dy)
return
end subroutine

