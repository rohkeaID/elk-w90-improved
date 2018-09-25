
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phscdvs(p,vsmt_,vsir_)
use modmain
use modphonon
use modstore
implicit none
! arguments
integer, intent(in) :: p
real(8), intent(in) :: vsmt_(npmtmax,natmtot),vsir_(ngtot)
! local variables
integer is,ia,ja,ias,jas
integer nr,nri,np,i
integer iv(3),ig_,ifg_,ifg
real(8) vl(3),vc(3),t1
complex(8) z0,z1,z2
! allocatable arrays
real(8), allocatable :: rfmt(:)
complex(8), allocatable :: zfmt(:),zfir(:)
! prefactor
z0=1.d0/(2.d0*deltaph)
! multiply by i for sin-like displacement
if (p.eq.1) z0=z0*zi
!------------------------------!
!     muffin-tin potential     !
!------------------------------!
allocate(rfmt(npmtmax),zfmt(npmtmax))
z1=z0/dble(nscph)
ias=0
jas=0
do is=1,nspecies
  nr=nrmt(is)
  nri=nrmti(is)
  np=npmt(is)
  ja=0
  do ia=1,natoms_(is)
    ias=ias+1
    do i=1,nscph
      ja=ja+1
      jas=jas+1
! compute the difference between the perturbed and unperturbed potentials
      rfmt(1:np)=vsmt(1:np,jas)-vsmt_(1:np,jas)
! convert real potential difference to a complex spherical harmonic expansion
      call rtozfmt(nr,nri,rfmt,zfmt)
! the muffin-tin potential should have an *explicit* phase exp(iq.r)
      t1=-dot_product(vqc(:,iqph),vscph(:,i))
      z2=z1*cmplx(cos(t1),sin(t1),8)
! add to total
      dvsmt(1:np,ias)=dvsmt(1:np,ias)+z2*zfmt(1:np)
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
zfir(:)=z0*(vsir(:)-vsir_(:))
call zfftifc(3,ngridg,-1,zfir)
! convert to G+q-space
do ig_=1,ngtot_
  ifg_=igfft_(ig_)
  vl(:)=dble(ivg_(:,ig_))+vql(:,iqph)
  call r3mv(bvec_,vl,vc)
  call r3mv(binv,vc,vl)
  iv(:)=nint(vl(:))
  if ((iv(1).ge.intgv(1,1)).and.(iv(1).le.intgv(2,1)).and. &
      (iv(2).ge.intgv(1,2)).and.(iv(2).le.intgv(2,2)).and. &
      (iv(3).ge.intgv(1,3)).and.(iv(3).le.intgv(2,3))) then
    ifg=igfft(ivgig(iv(1),iv(2),iv(3)))
    dvsir(ifg_)=dvsir(ifg_)+zfir(ifg)
  else
    dvsir(ifg_)=0.d0
  end if
end do
! Fourier transform back to real-space
if (p.eq.1) call zfftifc(3,ngridg_,1,dvsir)
deallocate(zfir)
return
end subroutine

