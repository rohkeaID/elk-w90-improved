
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: fsmfield
! !INTERFACE:
subroutine fsmfield
! !USES:
use modmain
! !DESCRIPTION:
!   Updates the effective magnetic field, ${\bf B}_{\rm FSM}$, required for
!   fixing the spin moment to a given value, ${\bf M}_{\rm FSM}$. This is done
!   by adding a vector to the field which is proportional to the difference
!   between the moment calculated in the $i$th self-consistent loop and the
!   required moment:
!   $$ {\bf B}_{\rm FSM}^{i+1}={\bf B}_{\rm FSM}^i+\lambda\left({\bf M}^i
!    -{\bf M}_{\rm FSM}\right), $$
!   where $\lambda$ is a scaling factor.
!
! !REVISION HISTORY:
!   Created March 2005 (JKD)
!EOP
!BOC
implicit none
! local variables
integer idm,is,ia,ias,ir
real(8) v1(3),v2(3),t1
if ((.not.spinpol).or.(fsmtype.eq.0)) return
! fixed spin direction not valid for collinear magnetism
if ((.not.ncmag).and.(fsmtype.lt.0)) return
! determine the global effective field
if ((abs(fsmtype).eq.1).or.(abs(fsmtype).eq.3)) then
  if (ncmag) then
    v1(:)=momtot(:)
  else
    v1(:)=0.d0
    v1(3)=momtot(1)
  end if
  v2(:)=v1(:)-momfix(:)
  if (ncmag) then
    bfsmc(:)=bfsmc(:)+taufsm*v2(:)
  else
    bfsmc(1)=bfsmc(1)+taufsm*v2(3)
  end if
! make sure that the constraining field is perpendicular to the fixed moment
! for fixed direction calculations (Y. Kvashnin and LN)
  if (fsmtype.lt.0) call r3vo(momfix,bfsmc)
! add to the Kohn-Sham field
  do idm=1,ndmag
    t1=bfsmc(idm)
    do ias=1,natmtot
      is=idxis(ias)
      bsmt(1:npcmt(is),ias,idm)=bsmt(1:npcmt(is),ias,idm)+t1
    end do
    do ir=1,ngtot
      bsir(ir,idm)=bsir(ir,idm)+t1*cfunir(ir)
    end do
  end do
end if
if ((abs(fsmtype).eq.2).or.(abs(fsmtype).eq.3)) then
! determine the muffin-tin fields for fixed local moments
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
! if any component is >= 1000 then do not fix the moment
      t1=sum(abs(mommtfix(:,ia,is)))
      if (t1.ge.1000.d0) cycle
      if (ncmag) then
        v1(:)=mommt(:,ias)
      else
        v1(:)=0.d0
        v1(3)=mommt(1,ias)
      end if
      v2(:)=v1(:)-mommtfix(:,ia,is)
      if (ncmag) then
        bfsmcmt(:,ias)=bfsmcmt(:,ias)+taufsm*v2(:)
      else
        bfsmcmt(1,ias)=bfsmcmt(1,ias)+taufsm*v2(3)
      end if
! fixed spin direction
      if (fsmtype.lt.0) call r3vo(mommtfix(:,ia,is),bfsmcmt(:,ias))
! add to the Kohn-Sham muffin-tin field
      do idm=1,ndmag
        t1=bfsmcmt(idm,ias)
        bsmt(1:npcmt(is),ias,idm)=bsmt(1:npcmt(is),ias,idm)+t1
      end do
    end do
  end do
end if
return
end subroutine
!EOC

