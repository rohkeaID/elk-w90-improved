
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: linengy
! !INTERFACE:
subroutine linengy
! !USES:
use modmain
use modmpi
! !DESCRIPTION:
!   Calculates the new linearisation energies for both the APW and local-orbital
!   radial functions. See the routine {\tt findband}.
!
! !REVISION HISTORY:
!   Created May 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
logical fnd
integer is,ia,ja,ias,jas
integer l,ilo,io,jo,nnf
! automatic arrays
logical done(natmmax)
real(8) vr(nrmtmax)
nnf=0
! begin loops over atoms and species
do is=1,nspecies
  done(:)=.false.
  do ia=1,natoms(is)
    if (done(ia)) cycle
    ias=idxas(ia,is)
    vr(1:nrmt(is))=vsmt(1,1:nrmt(is),ias)*y00
!-----------------------!
!     APW functions     !
!-----------------------!
    do l=0,lmaxapw
      do io=1,apword(l,is)
        if (apwve(io,l,is)) then
! check if previous radial functions have same default energies
          do jo=1,io-1
            if (apwve(jo,l,is)) then
              if (abs(apwe0(io,l,is)-apwe0(jo,l,is)).lt.1.d-4) then
                apwe(io,l,ias)=apwe(jo,l,ias)
                goto 10
              end if
            end if
          end do
! find the band energy starting from default
          apwe(io,l,ias)=apwe0(io,l,is)
          call findband(solsc,l,nrmt(is),rsp(:,is),vr,epsband,demaxbnd, &
           apwe(io,l,ias),fnd)
          if (.not.fnd) nnf=nnf+1
        else
! set linearisation energy automatically
          if (autolinengy) apwe(io,l,ias)=efermi+dlefe
        end if
10 continue
      end do
    end do
!---------------------------------!
!     local-orbital functions     !
!---------------------------------!
    do ilo=1,nlorb(is)
      do io=1,lorbord(ilo,is)
        if (lorbve(io,ilo,is)) then
! check if previous radial functions have same default energies
          do jo=1,io-1
            if (lorbve(jo,ilo,is)) then
              if (abs(lorbe0(io,ilo,is)-lorbe0(jo,ilo,is)).lt.1.d-4) then
                lorbe(io,ilo,ias)=lorbe(jo,ilo,ias)
                goto 20
              end if
            end if
          end do
          l=lorbl(ilo,is)
! find the band energy starting from default
          lorbe(io,ilo,ias)=lorbe0(io,ilo,is)
          call findband(solsc,l,nrmt(is),rsp(:,is),vr,epsband,demaxbnd, &
           lorbe(io,ilo,ias),fnd)
          if (.not.fnd) nnf=nnf+1
        else
! set linearisation energy automatically
          if (autolinengy) lorbe(io,ilo,ias)=efermi+dlefe
        end if
20 continue
      end do
    end do
    done(ia)=.true.
! copy to equivalent atoms
    do ja=1,natoms(is)
      if ((.not.done(ja)).and.(eqatoms(ia,ja,is))) then
        jas=idxas(ja,is)
        do l=0,lmaxapw
          do io=1,apword(l,is)
            apwe(io,l,jas)=apwe(io,l,ias)
          end do
        end do
        do ilo=1,nlorb(is)
          do io=1,lorbord(ilo,is)
            lorbe(io,ilo,jas)=lorbe(io,ilo,ias)
          end do
        end do
        done(ja)=.true.
      end if
    end do
! end loops over atoms and species
  end do
end do
if ((nnf.gt.0).and.mp_mpi) then
  write(*,*)
  write(*,'("Warning(linengy): could not find ",I3," linearisation energies &
   &in s.c. loop ",I5)') nnf,iscl
end if
return
end subroutine
!EOC
