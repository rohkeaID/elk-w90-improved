
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: olprad
! !INTERFACE:
subroutine olprad
! !USES:
use modmain
! !DESCRIPTION:
!   Calculates the radial overlap integrals of the APW and local-orbital basis
!   functions. In other words, for atom $\alpha$, it computes integrals of the
!   form
!   $$ o^{\alpha}_{qp}=\int_0^{R_i}u^{\alpha}_{q;l_p}(r)v^{\alpha}_p(r)r^2dr $$
!   and
!   $$ o^{\alpha}_{pp'}=\int_0^{R_i}v^{\alpha}_p(r)v^{\alpha}_{p'}(r)r^2dr,
!    \quad l_p=l_{p'} $$
!   where $u^{\alpha}_{q;l}$ is the $q$th APW radial function for angular
!   momentum $l$; and $v^{\alpha}_p$ is the $p$th local-orbital radial function
!   and has angular momentum $l_p$.
!
! !REVISION HISTORY:
!   Created November 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ias,nr
integer ilo,jlo,l,io
real(8) t1
! automatic arrays
real(8) fr(nrmtmax)
! external functions
real(8) fintgt
external fintgt
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
!--------------------------------------!
!     APW-local-orbital integtrals     !
!--------------------------------------!
  do ilo=1,nlorb(is)
    l=lorbl(ilo,is)
    do io=1,apword(l,is)
      fr(1:nr)=apwfr(1:nr,1,io,l,ias)*lofr(1:nr,1,ilo,ias)*r2sp(1:nr,is)
      t1=fintgt(-1,nr,rsp(:,is),fr)
      oalo(io,ilo,ias)=t1
    end do
  end do
!-----------------------------------------------!
!     local-orbital-local-orbital integrals     !
!-----------------------------------------------!
  do ilo=1,nlorb(is)
    l=lorbl(ilo,is)
    do jlo=1,nlorb(is)
      if (lorbl(jlo,is).eq.l) then
        fr(1:nr)=lofr(1:nr,1,ilo,ias)*lofr(1:nr,1,jlo,ias)*r2sp(1:nr,is)
        t1=fintgt(-1,nr,rsp(:,is),fr)
        ololo(ilo,jlo,ias)=t1
      end if
    end do
  end do
end do
return
end subroutine
!EOC

