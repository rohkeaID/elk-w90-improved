
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: match
! !INTERFACE:
subroutine match(ngp,gpc,tpgpc,sfacgp,apwalm)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ngp    : number of G+p-vectors (in,integer)
!   gpc    : length of G+p-vectors (in,real(ngkmax))
!   tpgpc  : (theta, phi) coordinates of G+p-vectors (in,real(2,ngkmax))
!   sfacgp : structure factors of G+p-vectors (in,complex(ngkmax,natmtot))
!   apwalm : APW matching coefficients
!            (out,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
! !DESCRIPTION:
!   Computes the $({\bf G+p})$-dependent matching coefficients for the APW basis
!   functions. Inside muffin-tin $\alpha$, the APW functions are given by
!   $$ \phi^{\alpha}_{\bf G+p}({\bf r})=\sum_{l=0}^{l_{\rm max}}
!    \sum_{m=-l}^{l}\sum_{j=1}^{M^{\alpha}_l}A^{\alpha}_{jlm}({\bf G+p})
!    u^{\alpha}_{jl}(r)Y_{lm}(\hat{{\bf r}}), $$
!   where $A^{\alpha}_{jlm}({\bf G+p})$ is the matching coefficient,
!   $M^{\alpha}_l$ is the order of the APW and $u^{\alpha}_{jl}$ is the radial
!   function. In the interstitial region, an APW function is a plane wave,
!   $\exp(i({\bf G+p})\cdot{\bf r})/\sqrt{\Omega}$, where $\Omega$ is the unit
!   cell volume. Ensuring continuity up to the $(M^{\alpha}_l-1)$th derivative
!   across the muffin-tin boundary therefore requires that the matching
!   coefficients satisfy
!   $$ \sum_{j=1}^{M^{\alpha}_l}D_{ij}A^{\alpha}_{jlm}({\bf G+p})=b_i\;, $$
!   where
!   $$ D_{ij}=\left.\frac{d^{i-1}u^{\alpha}_{jl}(r)}{dr^{i-1}}
!    \right|_{r=R_{\alpha}} $$
!   and
!   $$ b_i=\frac{4\pi i^l}{\sqrt{\Omega}}|{\bf G+p}|^{i-1}j^{(i-1)}_l
!    (|{\bf G+p}|R_{\alpha})\exp(i({\bf G+p})\cdot{\bf r}_{\alpha})Y^*_{lm}
!    (\widehat{{\bf G+p}}), $$
!   with ${\bf r}_{\alpha}$ the atomic position and $R_{\alpha}$ the muffin-tin
!   radius. See routine {\tt wavefmt}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!   Fixed documentation, June 2006 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ngp
real(8), intent(in) :: gpc(ngkmax),tpgpc(2,ngkmax)
complex(8), intent(in) :: sfacgp(ngkmax,natmtot)
complex(8), intent(out) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
! local variables
integer np,is,ia,ias,omax
integer l,m,lm,io,jo
integer i,ir,igp,info
real(8) t0,t1
complex(8) z1,z2
! allocatable arrays
integer, allocatable :: ipiv(:)
real(8), allocatable :: c(:),djl(:,:,:)
complex(8), allocatable :: ylmgp(:,:),a(:,:),b(:,:)
! external functions
real(8) polynom
external polynom
! polynomial order
np=max(apwordmax+1,4)
allocate(ipiv(np))
allocate(c(np))
allocate(djl(0:lmaxapw,apwordmax,ngp))
allocate(ylmgp(lmmaxapw,ngp))
allocate(a(apwordmax,apwordmax))
allocate(b(apwordmax,ngp*(2*lmaxapw+1)))
! compute the spherical harmonics of the G+p-vectors
do igp=1,ngp
  call genylm(lmaxapw,tpgpc(:,igp),ylmgp(:,igp))
end do
t0=fourpi/sqrt(omega)
! loop over species
do is=1,nspecies
! maximum APW order for this species
  omax=maxval(apword(1:lmaxapw,is))
! starting point on radial mesh for fitting polynomial of order np
  ir=nrmt(is)-np+1
! evaluate the spherical Bessel function derivatives for all G+p-vectors
  do igp=1,ngp
    t1=gpc(igp)*rmt(is)
    do io=1,omax
      call sbesseldm(io-1,lmaxapw,t1,djl(:,io,igp))
    end do
    t1=1.d0
    do io=2,omax
      t1=t1*gpc(igp)
      djl(:,io,igp)=t1*djl(:,io,igp)
    end do
  end do
! loop over atoms
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! begin loop over l
    do l=0,lmaxapw
      z1=t0*zil(l)
! set up matrix of derivatives
      do jo=1,apword(l,is)
        do io=1,apword(l,is)
          a(io,jo)=polynom(io-1,np,rsp(ir,is),apwfr(ir,1,jo,l,ias),c,rmt(is))
        end do
      end do
! set up target vectors
      i=0
      do igp=1,ngp
        z2=z1*sfacgp(igp,ias)
        do m=-l,l
          lm=idxlm(l,m)
          i=i+1
          do io=1,apword(l,is)
            b(io,i)=djl(l,io,igp)*z2*conjg(ylmgp(lm,igp))
          end do
        end do
      end do
! solve the general complex linear systems
      call zgesv(apword(l,is),i,a,apwordmax,ipiv,b,apwordmax,info)
      if (info.ne.0) then
        write(*,*)
        write(*,'("Error(match): could not find APW matching coefficients")')
        write(*,'(" for species ",I4," and atom ",I4)') is,ia
        write(*,'(" ZGESV returned INFO = ",I8)') info
        write(*,*)
        stop
      end if
      i=0
      do igp=1,ngp
        do m=-l,l
          lm=idxlm(l,m)
          i=i+1
          do io=1,apword(l,is)
            apwalm(igp,io,lm,ias)=b(io,i)
          end do
        end do
      end do
! end loop over l
    end do
! end loops over atoms and species
  end do
end do
deallocate(ipiv,c,djl,ylmgp,a,b)
return
end subroutine
!EOC
