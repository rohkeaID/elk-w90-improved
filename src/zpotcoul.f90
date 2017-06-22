
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: zpotcoul
! !INTERFACE:
subroutine zpotcoul(nr,nri,ld1,r,igp0,gpc,jlgpr,ylmgp,sfacgp,zrhoir,ld2, &
 zvclmt,zvclir,zrho0)
! !USES:
use modmain
use modphonon
! !INPUT/OUTPUT PARAMETERS:
!   nr     : number of radial points for each species (in,integer(nspecies))
!   nri    : number of points on inner part of muffin-tin (in,integer(nspecies))
!   ld1    : leading dimension (in,integer)
!   r      : radial mesh for each species (in,real(ld1,nspecies))
!   igp0   : index of the shortest G+p-vector (in,integer)
!   gpc    : G+p-vector lengths (in,real(ngvec))
!   jlgpr  : spherical Bessel functions for evergy G+p-vector and muffin-tin
!            radius (in,real(0:lnpsd,ngvec,nspecies))
!   ylmgp  : spherical harmonics of the G+p-vectors (in,complex(lmmaxvr,ngvec))
!   sfacgp : structure factors of the G+p-vectors (in,complex(ngvec,natmtot))
!   zrhoir : interstitial charge density (in,complex(ngtot))
!   ld2    : leading dimension (in,integer)
!   zvclmt : muffin-tin Coulomb potential, with the contribution from the
!            isolated muffin-tin density precalculated and passed in
!            (inout,complex(lmmaxvr,ld2,natmtot))
!   zvclir : interstitial Coulomb potential (out,complex(ngtot))
!   zrho0  : G+p=0 term of the pseudocharge density (out,complex)
! !DESCRIPTION:
!   Calculates the Coulomb potential of a complex charge density by solving
!   Poisson's equation using the method of M. Weinert, {\it J. Math. Phys.}
!   {\bf 22}, 2433 (1981). First, the multipole moments of the muffin-tin charge
!   are determined for the $j$th atom of the $i$th species by
!   $$ q_{ij;lm}^{\rm MT}=\int_0^{R_i}r^{l+2}\rho_{ij;lm}(r)dr+z_{ij}Y_{00}
!   \,\delta_{l,0}\;, $$
!   where $R_i$ is the muffin-tin radius and $z_{ij}$ is a point charge located
!   at the atom center (usually the nuclear charge, which should be taken as
!   {\bf negative}). Next, the multipole moments of the continuation of the
!   interstitial density, $\rho^{\rm I}$, into the muffin-tin are found with
!   $$ q_{ij;lm}^{\rm I}=4\pi i^l R_i^{l+3}\sum_{\bf G}\frac{j_{l+1}(GR_i)}
!    {GR_i}\rho^{\rm I}({\bf G})\exp(i{\bf G}\cdot{\bf r}_{ij})Y_{lm}^*
!    (\hat{\bf G}), $$
!   remembering that
!   $$ \lim_{x\rightarrow 0}\frac{j_{l+n}(x)}{x^n}=\frac{1}{(2n+1)!!}
!    \delta_{l,0} $$
!   should be used for the case ${\bf G}=0$. A pseudocharge is now constructed
!   which is equal to the real density in the interstitial region and whose
!   multipoles are the difference between the real and interstitial muffin-tin
!   multipoles. This pseudocharge density is smooth in the sense that it can be
!   expanded in terms of the finite set of ${\bf G}$-vectors. In each muffin-tin
!   the pseudocharge has the form
!   $$ \rho_{ij}^{\rm P}({\bf r})=\rho^{\rm I}({\bf r}-{\bf r}_{ij})+\sum_{lm}
!    \rho_{ij;lm}^{\rm P}\frac{1}{R_i^{l+3}}\left(\frac{r}{R_i}\right)^l\left(1-
!    \frac{r^2}{R_i^2}\right)^{N_i}Y_{lm}(\hat{\bf r}) $$
!   where
!   $$ \rho_{ij;lm}^{\rm P}=\frac{(2l+2N_i+3)!!}{2^N_iN_i!(2l+1)!!}\left(
!    q_{ij;lm}^{\rm MT}-q_{ij;lm}^{\rm I}\right) $$
!   and $N_i\approx\frac{1}{2}R_iG_{\rm max}$ is generally a good choice.
!   The pseudocharge in reciprocal space is given by
!   $$ \rho^{\rm P}({\bf G})=\rho^{\rm I}({\bf G})+\sum_{ij;lm}2^{N_i}N_i!
!    \frac{4\pi(-i)^l}{\Omega R_i^l}\frac{j_{l+N_i+1}(GR_i)}{(GR_i)^{N_i+1}}
!    \rho_{ij;lm}^{\rm P}\exp(-i{\bf G}\cdot{\bf r}_{ij})Y_{lm}(\hat{\bf G}) $$
!   which may be used for solving Poisson's equation directly
!   $$ V^{\rm P}({\bf G})=\begin{cases}
!     4\pi\frac{\rho^{\rm P}({\bf G})}{G^2} & G>0 \\
!     0 & G=0 \end{cases}\;. $$
!   The usual Green's function approach is then employed to determine the
!   potential in the muffin-tin sphere due to charge in the sphere. In other
!   words
!   $$ V_{ij;lm}^{\rm MT}(r)=\frac{4\pi}{2l+1}\left(\frac{1}{r^{l+1}}\int_0^r
!    \rho_{ij;lm}^{\rm MT}(r'){r'}^{l+2}dr'+r^l\int_r^{R_i}\frac{
!    \rho_{ij;lm}^{\rm MT}(r')}{{r'}^{l-1}}dr'\right)+\frac{1}{Y_{00}}
!    \frac{z_{ij}}{r}\delta_{l,0} $$
!   where the last term is the monopole arising from the point charge. All that
!   remains is to add the homogenous solution of Poisson's equation,
!   $$ V_{ij}^{\rm H}({\bf r})=\sum_{lm}V_{ij;lm}^{\rm H}\left(\frac{r}
!    {R_i}\right)^lY_{lm}(\hat{\bf r}), $$
!   to the muffin-tin potential so that it is continuous at the muffin-tin
!   boundary. Therefore the coefficients, $\rho_{ij;lm}^{\rm H}$, are given by
!   $$ V_{ij;lm}^{\rm H}=4\pi i^l\sum_{\bf G}j_{l}(Gr)V^{\rm P}({\bf G})
!    \exp(i{\bf G}\cdot{\bf r}_{ij})Y_{lm}^*(\hat{\bf G})-V_{ij;lm}^{\rm MT}
!    (R_i). $$
!   Finally note that the ${\bf G}$-vectors passed to the routine can represent
!   vectors with a non-zero offset, ${\bf G}+{\bf p}$ say, which is required for
!   calculating Coulomb matrix elements.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: nr(nspecies),nri(nspecies)
integer, intent(in) :: ld1
real(8), intent(in) :: r(ld1,nspecies)
integer, intent(in) :: igp0
real(8), intent(in) :: gpc(ngvec)
real(8), intent(in) :: jlgpr(0:lnpsd,ngvec,nspecies)
complex(8), intent(in) :: ylmgp(lmmaxvr,ngvec)
complex(8), intent(in) :: sfacgp(ngvec,natmtot)
complex(8), intent(in) :: zrhoir(ngtot)
integer, intent(in) :: ld2
complex(8), intent(inout) :: zvclmt(lmmaxvr,ld2,natmtot)
complex(8), intent(out) :: zvclir(ngtot)
complex(8), intent(out) :: zrho0
! local variables
integer is,ia,ias,iro,ir
integer l,m,lm,ig,ifg
real(8) t0,t1,t2,t3
complex(8) zsum1,zsum2,z1,z2
! automatic arrays
real(8) rmtl(0:lmaxvr+3,nspecies)
real(8) rl(ld1,0:lmaxvr)
complex(8) qlm(lmmaxvr,natmtot)
complex(8) zl(0:lmaxvr),zlm(lmmaxvr)
! external functions
real(8) factnm
external factnm
! compute (R_mt)^l
do is=1,nspecies
  rmtl(0,is)=1.d0
  do l=1,lmaxvr+3
    rmtl(l,is)=rmtl(l-1,is)*rmt(is)
  end do
end do
! compute the multipole moments from the muffin-tin potentials
do ias=1,natmtot
  is=idxis(ias)
  lm=0
  do l=0,lmaxvr
    t1=dble(2*l+1)*rmtl(l+1,is)/fourpi
    do m=-l,l
      lm=lm+1
      qlm(lm,ias)=t1*zvclmt(lm,nr(is),ias)
    end do
  end do
end do
! Fourier transform density to G-space and store in zvclir
call zcopy(ngtot,zrhoir,1,zvclir,1)
call zfftifc(3,ngridg,-1,zvclir)
! subtract the multipole moments of the interstitial charge density
do is=1,nspecies
  do l=0,lmaxvr
    zl(l)=fourpi*zil(l)*rmtl(l+2,is)
  end do
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ig=1,ngvec
      ifg=igfft(ig)
      if (gpc(ig).gt.epslat) then
        z1=zvclir(ifg)*sfacgp(ig,ias)/gpc(ig)
        lm=0
        do l=0,lmaxvr
          z2=jlgpr(l+1,ig,is)*z1*zl(l)
          do m=-l,l
            lm=lm+1
            qlm(lm,ias)=qlm(lm,ias)-z2*conjg(ylmgp(lm,ig))
          end do
        end do
      else
        t1=(fourpi/3.d0)*rmtl(3,is)*y00
        qlm(1,ias)=qlm(1,ias)-t1*zvclir(ifg)
      end if
    end do
  end do
end do
! find the smooth pseudocharge within the muffin-tin whose multipoles are the
! difference between the real muffin-tin and interstitial multipoles
t0=fourpi/omega
t1=factnm(2*lnpsd+1,2)
do ias=1,natmtot
  is=idxis(ias)
  lm=0
  do l=0,lmaxvr
    t2=t1/(factnm(2*l+1,2)*rmtl(l,is))
    z1=t2*zilc(l)
    do m=-l,l
      lm=lm+1
      zlm(lm)=z1*qlm(lm,ias)
    end do
  end do
! add the pseudocharge and real interstitial densities in G-space
  do ig=1,ngvec
    ifg=igfft(ig)
    if (gpc(ig).gt.epslat) then
      t2=gpc(ig)*rmt(is)
      t3=1.d0/t2**lnpsd
      zsum1=t3*zlm(1)*ylmgp(1,ig)
      lm=1
      do l=1,lmaxvr
        lm=lm+1
        zsum2=zlm(lm)*ylmgp(lm,ig)
        do m=1-l,l
          lm=lm+1
          zsum2=zsum2+zlm(lm)*ylmgp(lm,ig)
        end do
        t3=t3*t2
        zsum1=zsum1+t3*zsum2
      end do
      z1=t0*jlgpr(lnpsd,ig,is)*conjg(sfacgp(ig,ias))
      zvclir(ifg)=zvclir(ifg)+z1*zsum1
    else
      t2=t0*y00/factnm(2*lnpsd+1,2)
      zvclir(ifg)=zvclir(ifg)+t2*zlm(1)
    end if
  end do
end do
! solve Poisson's equation in G+p-space for the pseudocharge
do ig=1,igp0-1
  ifg=igfft(ig)
  zvclir(ifg)=(fourpi/gpc(ig)**2)*zvclir(ifg)
end do
do ig=igp0+1,ngvec
  ifg=igfft(ig)
  zvclir(ifg)=(fourpi/gpc(ig)**2)*zvclir(ifg)
end do
! set zrho0 (pseudocharge density coefficient of the smallest G+p-vector)
if (igp0.gt.0) then
  ifg=igfft(igp0)
  zrho0=zvclir(ifg)
  zvclir(ifg)=0.d0
end if
! match potentials at muffin-tin boundary by adding homogeneous solution
do is=1,nspecies
! compute (r/R_mt)^l
  t1=1.d0/rmt(is)
  do ir=1,nr(is)
    t2=t1*r(ir,is)
    rl(ir,0)=1.d0
    do l=1,lmaxvr
      rl(ir,l)=rl(ir,l-1)*t2
    end do
  end do
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! find the spherical harmonic expansion of the interstitial potential at the
! muffin-tin radius
    zlm(:)=0.d0
    do ig=1,ngvec
      ifg=igfft(ig)
      z1=fourpi*zvclir(ifg)*sfacgp(ig,ias)
      lm=0
      do l=0,lmaxvr
        z2=jlgpr(l,ig,is)*z1*zil(l)
        do m=-l,l
          lm=lm+1
          zlm(lm)=zlm(lm)+z2*conjg(ylmgp(lm,ig))
        end do
      end do
    end do
! add homogenous solution
    lm=0
    do l=0,lmaxvr
      if (l.le.lmaxinr) then
        iro=1
      else
        iro=nri(is)+1
      end if
      do m=-l,l
        lm=lm+1
        z1=zlm(lm)-zvclmt(lm,nr(is),ias)
! store the nuclear potential without the self-term for the phonon dynamical
! matrix calculation
        if (tphdyn) then
          if (ias.eq.iasph) then
            zvnmt(lm,iro:nr(is))=z1*rl(iro:nr(is),l)
          end if
        end if
        do ir=iro,nr(is)
          zvclmt(lm,ir,ias)=zvclmt(lm,ir,ias)+z1*rl(ir,l)
        end do
      end do
    end do
  end do
end do
! Fourier transform interstitial potential to real-space
call zfftifc(3,ngridg,1,zvclir)
return
end subroutine
!EOC

