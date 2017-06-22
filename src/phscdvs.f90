
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phscdvs(p,vsmt0,vsir0)
use modmain
use modphonon
use modstore
implicit none
! arguments
integer, intent(in) :: p
real(8), intent(in) :: vsmt0(lmmaxvr,nrmtmax,natmtot),vsir0(ngtot)
! local variables
integer is,ia,ja,ias,jas
integer nr,nri,iv(3),i
integer ig0,ifg0,ifg
real(8) vl(3),vc(3),t1
complex(8) z0,z1,z2
! allocatable arrays
real(8), allocatable :: rfmt(:,:)
complex(8), allocatable :: zfmt(:,:),zfir(:)
! prefactor
z0=1.d0/(2.d0*deltaph)
! multiply by i for sin-like displacement
if (p.eq.1) z0=z0*zi
!------------------------------!
!     muffin-tin potential     !
!------------------------------!
allocate(rfmt(lmmaxvr,nrmtmax),zfmt(lmmaxvr,nrmtmax))
z1=z0/dble(nscph)
ias=0
jas=0
do is=1,nspecies
  nr=nrmt(is)
  nri=nrmtinr(is)
  ja=0
  do ia=1,natoms0(is)
    ias=ias+1
    do i=1,nscph
      ja=ja+1
      jas=jas+1
! compute the difference between the perturbed and unperturbed potentials
      rfmt(:,1:nr)=vsmt(:,1:nr,jas)-vsmt0(:,1:nr,jas)
! convert real potential difference to a complex spherical harmonic expansion
      call rtozfmt(nr,nri,1,rfmt,1,zfmt)
! the muffin-tin potential should have an *explicit* phase exp(iq.r)
      t1=-dot_product(vqc(:,iqph),vscph(:,i))
      z2=z1*cmplx(cos(t1),sin(t1),8)
! add to total
      call zfmtadd(nr,nri,z2,zfmt,dvsmt(:,:,ias))
    end do
! end loop over atoms and species
  end do
end do
deallocate(rfmt,zfmt)
!--------------------------------!
!     interstitial potential     !
!--------------------------------!
! Fourier transform interstitial potential derivative to G-space
allocate(zfir(ngtot))
zfir(:)=z0*(vsir(:)-vsir0(:))
call zfftifc(3,ngridg,-1,zfir)
! convert to G+q-space
do ig0=1,ngtot0
  ifg0=igfft0(ig0)
  vl(:)=dble(ivg0(:,ig0))+vql(:,iqph)
  call r3mv(bvec0,vl,vc)
  call r3mv(binv,vc,vl)
  iv(:)=nint(vl(:))
  if ((iv(1).ge.intgv(1,1)).and.(iv(1).le.intgv(2,1)).and. &
      (iv(2).ge.intgv(1,2)).and.(iv(2).le.intgv(2,2)).and. &
      (iv(3).ge.intgv(1,3)).and.(iv(3).le.intgv(2,3))) then
    ifg=igfft(ivgig(iv(1),iv(2),iv(3)))
    dvsir(ifg0)=dvsir(ifg0)+zfir(ifg)
  else
    dvsir(ifg0)=0.d0
  end if
end do
! Fourier transform back to real-space
if (p.eq.1) call zfftifc(3,ngridg0,1,dvsir)
deallocate(zfir)
return
end subroutine

