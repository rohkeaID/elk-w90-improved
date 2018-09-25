
! Copyright (C) 2002-2016 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: hmlrad
! !INTERFACE:
subroutine hmlrad
! !USES:
use modmain
use modomp
! !DESCRIPTION:
!   Calculates the radial Hamiltonian integrals of the APW and local-orbital
!   basis functions. In other words, for atom $\alpha$, it computes integrals of
!   the form
!   $$ h^{\alpha}_{qq';ll'l''m''}=\begin{cases}
!    \int_0^{R_i}u^{\alpha}_{q;l}(r)H u^{\alpha}_{q';l'}(r)r^2dr & l''=0 \\
!    \int_0^{R_i}u^{\alpha}_{q;l}(r)V^{\alpha}_{l''m''}(r)
!    u^{\alpha}_{q';l'}(r)r^2dr & l''>0 \end{cases}, $$
!   where $u^{\alpha}_{q;l}$ is the $q$th APW radial function for angular
!   momentum $l$; $H$ is the Hamiltonian of the radial Schr\"{o}dinger equation;
!   and $V^{\alpha}_{l''m''}$ is the muffin-tin Kohn-Sham potential. Similar
!   integrals are calculated for APW-local-orbital and
!   local-orbital-local-orbital contributions.
!
! !REVISION HISTORY:
!   Created December 2003 (JKD)
!   Updated for compressed muffin-tin functions, March 2016 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ias,nthd
integer nr,nri,nro
integer iro,ir,npi,i
integer l1,l2,l3,m2,lm2
integer io,jo,ilo,jlo
real(8) t1
! allocatable arrays
real(8), allocatable :: fr(:)
! external functions
real(8) fintgt
external fintgt
! begin loops over atoms and species
call omp_hold(natmtot,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(fr,is,nr,nri,nro,iro,npi) &
!$OMP PRIVATE(l1,l2,l3,io,jo,ir,t1) &
!$OMP PRIVATE(lm2,m2,i,ilo,jlo) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ias=1,natmtot
  allocate(fr(nrmtmax))
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  nro=nr-nri
  iro=nri+1
  npi=npmti(is)
!---------------------------!
!     APW-APW integrals     !
!---------------------------!
  do l1=0,lmaxapw
    do io=1,apword(l1,is)
      do l3=0,lmaxapw
        do jo=1,apword(l3,is)
          if (l1.eq.l3) then
            do ir=1,nr
              fr(ir)=apwfr(ir,1,io,l1,ias)*apwfr(ir,2,jo,l3,ias)*r2sp(ir,is)
            end do
            t1=fintgt(-1,nr,rsp(:,is),fr)
            haa(1,jo,l3,io,l1,ias)=t1/y00
          else
            haa(1,jo,l3,io,l1,ias)=0.d0
          end if
          if (l1.ge.l3) then
            lm2=1
            do l2=1,lmaxi
              do m2=-l2,l2
                lm2=lm2+1
                i=lm2
                do ir=1,nri
                  t1=apwfr(ir,1,io,l1,ias)*apwfr(ir,1,jo,l3,ias)*r2sp(ir,is)
                  fr(ir)=t1*vsmt(i,ias)
                  i=i+lmmaxi
                end do
                do ir=iro,nr
                  t1=apwfr(ir,1,io,l1,ias)*apwfr(ir,1,jo,l3,ias)*r2sp(ir,is)
                  fr(ir)=t1*vsmt(i,ias)
                  i=i+lmmaxo
                end do
                t1=fintgt(-1,nr,rsp(:,is),fr)
                haa(lm2,jo,l3,io,l1,ias)=t1
                haa(lm2,io,l1,jo,l3,ias)=t1
              end do
            end do
            do l2=lmaxi+1,lmaxo
              do m2=-l2,l2
                lm2=lm2+1
                i=npi+lm2
                do ir=iro,nr
                  t1=apwfr(ir,1,io,l1,ias)*apwfr(ir,1,jo,l3,ias)*r2sp(ir,is)
                  fr(ir)=t1*vsmt(i,ias)
                  i=i+lmmaxo
                end do
                t1=fintgt(-1,nro,rsp(iro,is),fr(iro))
                haa(lm2,jo,l3,io,l1,ias)=t1
                haa(lm2,io,l1,jo,l3,ias)=t1
              end do
            end do
          end if
        end do
      end do
    end do
  end do
!-------------------------------------!
!     local-orbital-APW integrals     !
!-------------------------------------!
  do ilo=1,nlorb(is)
    l1=lorbl(ilo,is)
    do l3=0,lmaxapw
      do io=1,apword(l3,is)
        if (l1.eq.l3) then
          do ir=1,nr
            fr(ir)=lofr(ir,1,ilo,ias)*apwfr(ir,2,io,l3,ias)*r2sp(ir,is)
          end do
          t1=fintgt(-1,nr,rsp(:,is),fr)
          hloa(1,io,l3,ilo,ias)=t1/y00
        else
          hloa(1,io,l3,ilo,ias)=0.d0
        end if
        lm2=1
        do l2=1,lmaxi
          do m2=-l2,l2
            lm2=lm2+1
            i=lm2
            do ir=1,nri
              t1=lofr(ir,1,ilo,ias)*apwfr(ir,1,io,l3,ias)*r2sp(ir,is)
              fr(ir)=t1*vsmt(i,ias)
              i=i+lmmaxi
            end do
            do ir=nri+1,nr
              t1=lofr(ir,1,ilo,ias)*apwfr(ir,1,io,l3,ias)*r2sp(ir,is)
              fr(ir)=t1*vsmt(i,ias)
              i=i+lmmaxo
            end do
            t1=fintgt(-1,nr,rsp(:,is),fr)
            hloa(lm2,io,l3,ilo,ias)=t1
          end do
        end do
        do l2=lmaxi+1,lmaxo
          do m2=-l2,l2
            lm2=lm2+1
            i=npi+lm2
            do ir=iro,nr
              t1=lofr(ir,1,ilo,ias)*apwfr(ir,1,io,l3,ias)*r2sp(ir,is)
              fr(ir)=t1*vsmt(i,ias)
              i=i+lmmaxo
            end do
            t1=fintgt(-1,nro,rsp(iro,is),fr(iro))
            hloa(lm2,io,l3,ilo,ias)=t1
          end do
        end do
      end do
    end do
  end do
!-----------------------------------------------!
!     local-orbital-local-orbital integrals     !
!-----------------------------------------------!
  do ilo=1,nlorb(is)
    l1=lorbl(ilo,is)
    do jlo=1,nlorb(is)
      l3=lorbl(jlo,is)
      if (l1.eq.l3) then
        do ir=1,nr
          fr(ir)=lofr(ir,1,ilo,ias)*lofr(ir,2,jlo,ias)*r2sp(ir,is)
        end do
        t1=fintgt(-1,nr,rsp(:,is),fr)
        hlolo(1,jlo,ilo,ias)=t1/y00
      else
        hlolo(1,jlo,ilo,ias)=0.d0
      end if
      lm2=1
      do l2=1,lmaxi
        do m2=-l2,l2
          lm2=lm2+1
          i=lm2
          do ir=1,nri
            t1=lofr(ir,1,ilo,ias)*lofr(ir,1,jlo,ias)*r2sp(ir,is)
            fr(ir)=t1*vsmt(i,ias)
            i=i+lmmaxi
          end do
          do ir=iro,nr
            t1=lofr(ir,1,ilo,ias)*lofr(ir,1,jlo,ias)*r2sp(ir,is)
            fr(ir)=t1*vsmt(i,ias)
            i=i+lmmaxo
          end do
          t1=fintgt(-1,nr,rsp(:,is),fr)
          hlolo(lm2,jlo,ilo,ias)=t1
        end do
      end do
      do l2=lmaxi+1,lmaxo
        do m2=-l2,l2
          lm2=lm2+1
          i=npi+lm2
          do ir=iro,nr
            t1=lofr(ir,1,ilo,ias)*lofr(ir,1,jlo,ias)*r2sp(ir,is)
            fr(ir)=t1*vsmt(i,ias)
            i=i+lmmaxo
          end do
          t1=fintgt(-1,nro,rsp(iro,is),fr(iro))
          hlolo(lm2,jlo,ilo,ias)=t1
        end do
      end do
    end do
  end do
  deallocate(fr)
! end loops over atoms and species
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
return
end subroutine
!EOC

