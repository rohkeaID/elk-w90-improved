
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genapwfr
! !INTERFACE:
subroutine genapwfr
! !USES:
use modmain
! !DESCRIPTION:
!   Generates the APW radial functions. This is done by integrating the scalar
!   relativistic Schr\"{o}dinger equation (or its energy deriatives) at the
!   current linearisation energies using the spherical part of the Kohn-Sham
!   potential. The number of radial functions at each $l$-value is given by the
!   variable {\tt apword} (at the muffin-tin boundary, the APW functions have
!   continuous derivatives up to order ${\tt apword}-1$). Within each $l$, these
!   functions are orthonormalised with the Gram-Schmidt method. The radial
!   Hamiltonian is applied to the orthonormalised functions and the results are
!   stored in the global array {\tt apwfr}.
!
! !REVISION HISTORY:
!   Created March 2003 (JKD)
!   Copied to equivalent atoms, February 2010 (A. Kozhevnikov and JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ja,ias,jas
integer nr,ir,nn,l,io,jo
real(8) e,t1
! automatic arrays
logical done(natmmax)
real(8) vr(nrmtmax),fr(nrmtmax)
real(8) p0(nrmtmax,apwordmax),p1(nrmtmax),p1s(apwordmax)
real(8) q0(nrmtmax),q1(nrmtmax),ep0(nrmtmax,apwordmax)
! external functions
real(8) fintgt
external fintgt
do is=1,nspecies
  nr=nrmt(is)
  done(:)=.false.
  do ia=1,natoms(is)
    if (done(ia)) cycle
    ias=idxas(ia,is)
! use spherical part of potential
    vr(1:nr)=vsmt(1,1:nr,ias)*y00
    do l=0,lmaxapw
      do io=1,apword(l,is)
! linearisation energy accounting for energy derivative
        e=apwe(io,l,ias)+dble(apwdm(io,l,is))*deapwlo
! integrate the radial Schrodinger equation
        call rschrodint(solsc,l,e,nr,rsp(:,is),vr,nn,p0(:,io),p1,q0,q1)
        ep0(1:nr,io)=e*p0(1:nr,io)
! normalise radial functions
        fr(1:nr)=p0(1:nr,io)**2
        t1=fintgt(-1,nr,rsp(:,is),fr)
        t1=1.d0/sqrt(abs(t1))
        call dscal(nr,t1,p0(:,io),1)
        p1s(io)=t1*p1(nr)
        call dscal(nr,t1,ep0(:,io),1)
! subtract linear combination of previous vectors
        do jo=1,io-1
          fr(1:nr)=p0(1:nr,io)*p0(1:nr,jo)
          t1=-fintgt(-1,nr,rsp(:,is),fr)
          call daxpy(nr,t1,p0(:,jo),1,p0(:,io),1)
          p1s(io)=p1s(io)+t1*p1s(jo)
          call daxpy(nr,t1,ep0(:,jo),1,ep0(:,io),1)
        end do
! normalise radial functions again
        fr(1:nr)=p0(1:nr,io)**2
        t1=fintgt(-1,nr,rsp(:,is),fr)
        t1=abs(t1)
        if (t1.lt.1.d-25) then
          write(*,*)
          write(*,'("Error(genapwfr): degenerate APW radial functions")')
          write(*,'(" for species ",I4)') is
          write(*,'(" atom ",I4)') ia
          write(*,'(" angular momentum ",I4)') l
          write(*,'(" and order ",I4)') io
          write(*,*)
          stop
        end if
        t1=1.d0/sqrt(t1)
        call dscal(nr,t1,p0(:,io),1)
        p1s(io)=t1*p1s(io)
        call dscal(nr,t1,ep0(:,io),1)
! divide by r and store in global array
        do ir=1,nr
          t1=1.d0/rsp(ir,is)
          apwfr(ir,1,io,l,ias)=t1*p0(ir,io)
          apwfr(ir,2,io,l,ias)=t1*ep0(ir,io)
        end do
! derivative at the muffin-tin surface
        apwdfr(io,l,ias)=(p1s(io)-p0(nr,io)*t1)*t1
      end do
    end do
    done(ia)=.true.
! copy to equivalent atoms
    do ja=1,natoms(is)
      if ((.not.done(ja)).and.(eqatoms(ia,ja,is))) then
        jas=idxas(ja,is)
        do l=0,lmaxapw
          do io=1,apword(l,is)
            call dcopy(nr,apwfr(:,1,io,l,ias),1,apwfr(:,1,io,l,jas),1)
            call dcopy(nr,apwfr(:,2,io,l,ias),1,apwfr(:,2,io,l,jas),1)
            apwdfr(io,l,jas)=apwdfr(io,l,ias)
          end do
        end do
        done(ja)=.true.
      end if
    end do
! end loop over atoms and species
  end do
end do
return
end subroutine
!EOC
