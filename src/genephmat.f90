
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genephmat(iq,ik,de,dvphmt,dvphir,ephmat)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: iq,ik
real(8), intent(in) :: de
complex(8), intent(in) :: dvphmt(lmmaxvr,nrcmtmax,natmtot,nbph)
complex(8), intent(in) :: dvphir(ngtot,nbph)
complex(8), intent(out) :: ephmat(nstsv,nstsv,nbph)
! local variables
integer jk,jkq,isym,i,j,l
integer nst,nstq,ist,jst
real(8) vpql(3)
! automatic arrays
integer idx(nstsv),idxq(nstsv)
! allocatable arrays
complex(8), allocatable :: wfmt(:,:,:,:,:),wfir(:,:,:)
complex(8), allocatable :: wfmtq(:,:,:,:,:),wfirq(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:),zrhoir(:)
! external functions
complex(8) zfinp
external zfinp
! equivalent reduced k-point
jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
! k+q-vector in lattice coordinates
vpql(:)=vkl(:,ik)+vql(:,iq)
! find reduced k-point index corresponding to k+q
call findkpt(vpql,isym,jkq)
! index to states in energy window around Fermi energy
nst=0
nstq=0
do ist=1,nstsv
  if (abs(evalsv(ist,jk)-efermi).lt.de) then
    nst=nst+1
    idx(nst)=ist
  end if
  if (abs(evalsv(ist,jkq)-efermi).lt.de) then
    nstq=nstq+1
    idxq(nstq)=ist
  end if
end do
! generate the wavefunctions for all states at k and k+q
allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nst))
allocate(wfir(ngtot,nspinor,nst))
call genwfsvp(.false.,.false.,nst,idx,vkl(:,ik),wfmt,ngtot,wfir)
allocate(wfmtq(lmmaxvr,nrcmtmax,natmtot,nspinor,nstq))
allocate(wfirq(ngtot,nspinor,nstq))
call genwfsvp(.false.,.false.,nstq,idxq,vpql,wfmtq,ngtot,wfirq)
allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot),zrhoir(ngtot))
ephmat(:,:,:)=0.d0
do i=1,nstq
  ist=idxq(i)
  do j=1,nst
    jst=idx(j)
! note that the complex conjugate of the density is found because zfinp
! conjugates the first function
    call genzrho(.true.,.true.,wfmt(:,:,:,:,j),wfir(:,:,j),wfmtq(:,:,:,:,i), &
     wfirq(:,:,i),zrhomt,zrhoir)
    do l=1,nbph
      ephmat(ist,jst,l)=zfinp(zrhomt,zrhoir,dvphmt(:,:,:,l),dvphir(:,l))
    end do
  end do
end do
deallocate(wfmt,wfmtq,wfir,wfirq,zrhomt,zrhoir)
return
end subroutine

