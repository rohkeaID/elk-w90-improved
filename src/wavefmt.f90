
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: wavefmt
! !INTERFACE:
subroutine wavefmt(lrstp,ias,ngp,apwalm,evecfv,wfmt)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   lrstp  : radial step length (in,integer)
!   ias    : joint atom and species number (in,integer)
!   ngp    : number of G+p-vectors (in,integer)
!   apwalm : APW matching coefficients (in,complex(ngkmax,apwordmax,lmmaxapw))
!   evecfv : first-variational eigenvector (in,complex(nmatmax))
!   wfmt   : complex muffin-tin wavefunction passed in as real array
!            (out,real(2,*))
! !DESCRIPTION:
!   Calculates the first-variational wavefunction in the muffin-tin in terms of
!   a spherical harmonic expansion. For atom $\alpha$ and a particular $k$-point
!   ${\bf p}$, the $r$-dependent $(l,m)$-coefficients of the wavefunction for
!   the $i$th state are given by
!   $$ \Phi^{i{\bf p}}_{\alpha lm}(r)=\sum_{\bf G}b^{i{\bf p}}_{\bf G}
!    \sum_{j=1}^{M^{\alpha}_l}A^{\alpha}_{jlm}({\bf G+p})u^{\alpha}_{jl}(r)
!    +\sum_{j=1}^{N^{\alpha}}b^{i{\bf p}}_{(\alpha,j,m)}v^{\alpha}_j(r)
!    \delta_{l,l_j}, $$
!   where $b^{i{\bf p}}$ is the $i$th eigenvector returned from routine
!   {\tt eveqn}; $A^{\alpha}_{jlm}({\bf G+p})$ is the matching coefficient;
!   $M^{\alpha}_l$ is the order of the APW; $u^{\alpha}_{jl}$ is the APW radial
!   function; $N^{\alpha}$ is the number of local-orbitals; $v^{\alpha}_j$ is
!   the $j$th local-orbital radial function; and $(\alpha,j,m)$ is a compound
!   index for the location of the local-orbital in the eigenvector. See routines
!   {\tt genapwfr}, {\tt genlofr}, {\tt match} and {\tt eveqn}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!   Fixed description, October 2004 (C. Brouder)
!   Removed argument ist, November 2006 (JKD)
!   Changed arguments and optimised, December 2014 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: lrstp,ias,ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: evecfv(nmatmax)
real(8), intent(out) :: wfmt(2,*)
! local variables
integer is,ldi,ldo,io,ilo
integer nrc,nrci,nrco,iro
integer l,m,lm,npc,npci,i
complex(8) z1
! external functions
complex(8) zdotu
external zdotu
is=idxis(ias)
ldi=2*lmmaxi
ldo=2*lmmaxo
iro=nrmti(is)+lrstp
if (lrstp.eq.1) then
  nrc=nrmt(is)
  nrci=nrmti(is)
  npc=npmt(is)
  npci=npmti(is)
else if (lrstp.eq.lradstp) then
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
  npci=npcmti(is)
else
  write(*,*)
  write(*,'("Error(wavefmt): invalid lrstp : ",I8)') lrstp
  write(*,*)
  stop
end if
nrco=nrc-nrci
! zero the wavefunction
wfmt(:,1:npc)=0.d0
!-----------------------!
!     APW functions     !
!-----------------------!
lm=0
do l=0,lmaxo
  do m=-l,l
    lm=lm+1
    i=npci+lm
    do io=1,apword(l,is)
      z1=zdotu(ngp,evecfv,1,apwalm(:,io,lm),1)
      if (abs(dble(z1)).gt.1.d-14) then
        if (l.le.lmaxi) then
          call daxpy(nrci,dble(z1),apwfr(1,1,io,l,ias),lrstp,wfmt(1,lm),ldi)
        end if
        call daxpy(nrco,dble(z1),apwfr(iro,1,io,l,ias),lrstp,wfmt(1,i),ldo)
      end if
      if (abs(aimag(z1)).gt.1.d-14) then
        if (l.le.lmaxi) then
          call daxpy(nrci,aimag(z1),apwfr(1,1,io,l,ias),lrstp,wfmt(2,lm),ldi)
        end if
        call daxpy(nrco,aimag(z1),apwfr(iro,1,io,l,ias),lrstp,wfmt(2,i),ldo)
      end if
    end do
  end do
end do
!---------------------------------!
!     local-orbital functions     !
!---------------------------------!
do ilo=1,nlorb(is)
  l=lorbl(ilo,is)
  do m=-l,l
    lm=idxlm(l,m)
    i=npci+lm
    z1=evecfv(ngp+idxlo(lm,ilo,ias))
    if (abs(dble(z1)).gt.1.d-14) then
      if (l.le.lmaxi) then
        call daxpy(nrci,dble(z1),lofr(1,1,ilo,ias),lrstp,wfmt(1,lm),ldi)
      end if
      call daxpy(nrco,dble(z1),lofr(iro,1,ilo,ias),lrstp,wfmt(1,i),ldo)
    end if
    if (abs(aimag(z1)).gt.1.d-14) then
      if (l.le.lmaxi) then
        call daxpy(nrci,aimag(z1),lofr(1,1,ilo,ias),lrstp,wfmt(2,lm),ldi)
      end if
      call daxpy(nrco,aimag(z1),lofr(iro,1,ilo,ias),lrstp,wfmt(2,i),ldo)
    end if
  end do
end do
return
end subroutine
!EOC

