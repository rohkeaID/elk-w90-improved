
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genephmat(iq,ik,de,a,dvmt,dvir,ephmat)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: iq,ik
real(8), intent(in) :: de
complex(8), intent(in) :: a(nbph,nbph)
complex(8), intent(in) :: dvmt(npcmtmax,natmtot,nbph),dvir(ngtot,nbph)
complex(8), intent(out) :: ephmat(nstsv,nstsv,nbph)
! local variables
integer jk,jkq,isym,ld
integer nst,nstq,ist,jst
integer ispn,jspn,i,j
integer is,ia,ias,l
integer nrc,nrci,npc
integer igp,igpq,ifg
real(8) vpql(3)
! automatic arrays
integer idx(nstsv),idxq(nstsv)
integer ngp(nspnfv),ngpq(nspnfv)
complex(8) x(nbph)
! allocatable arrays
integer, allocatable :: igpig(:,:),igpqig(:,:)
complex(8), allocatable :: wfmt(:,:,:,:),wfgp(:,:,:)
complex(8), allocatable :: wfmtq(:,:,:,:),wfgpq(:,:,:)
complex(8), allocatable :: zfmt1(:),zfmt2(:)
complex(8), allocatable :: wfir1(:),wfir2(:),z(:)
! external functions
complex(8) zfmtinp,zdotc
external zfmtinp,zdotc
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
allocate(igpig(ngkmax,nspnfv))
allocate(wfmt(npcmtmax,natmtot,nspinor,nst),wfgp(ngkmax,nspinor,nst))
call genwfsvp(.false.,.true.,nst,idx,ngridg,igfft,vkl(:,ik),ngp,igpig,wfmt, &
 ngkmax,wfgp)
allocate(igpqig(ngkmax,nspnfv))
allocate(wfmtq(npcmtmax,natmtot,nspinor,nstq),wfgpq(ngkmax,nspinor,nstq))
call genwfsvp(.false.,.true.,nstq,idxq,ngridg,igfft,vpql,ngpq,igpqig,wfmtq, &
 ngkmax,wfgpq)
! zero the electron-phonon coupling matrix elements
ephmat(:,:,:)=0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
allocate(zfmt1(npcmtmax),zfmt2(npcmtmax))
do i=1,nstq
  ist=idxq(i)
  do j=1,nst
    jst=idx(j)
    x(:)=0.d0
    do is=1,nspecies
      nrc=nrcmt(is)
      nrci=nrcmti(is)
      npc=npcmt(is)
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        if (spinpol) then
          zfmt1(1:npc)=wfmtq(1:npc,ias,1,i)*conjg(wfmt(1:npc,ias,1,j)) &
                      +wfmtq(1:npc,ias,2,i)*conjg(wfmt(1:npc,ias,2,j))
        else
          zfmt1(1:npc)=wfmtq(1:npc,ias,1,i)*conjg(wfmt(1:npc,ias,1,j))
        end if
        call zfsht(nrc,nrci,zfmt1,zfmt2)
        do l=1,nbph
          x(l)=x(l)+zfmtinp(nrc,nrci,rcmt(:,is),r2cmt(:,is),zfmt2,dvmt(:,ias,l))
        end do
      end do
    end do
    ephmat(ist,jst,:)=x(:)
  end do
end do
deallocate(zfmt1,zfmt2)
!---------------------------!
!     interstitial part     !
!---------------------------!
allocate(wfir1(ngtot),wfir2(ngtot),z(ngkmax))
do j=1,nst
  jst=idx(j)
  do ispn=1,nspinor
    jspn=jspnfv(ispn)
! Fourier transform wavefunction to real-space
    wfir1(:)=0.d0
    do igp=1,ngp(jspn)
      ifg=igfft(igpig(igp,jspn))
      wfir1(ifg)=wfgp(igp,ispn,j)
    end do
    call zfftifc(3,ngridg,1,wfir1)
    do l=1,nbph
! apply potential derivative to wavefunction
      wfir2(:)=dvir(:,l)*wfir1(:)
! Fourier transform to G+p+q-space
      call zfftifc(3,ngridg,-1,wfir2)
      do igpq=1,ngpq(jspn)
        ifg=igfft(igpqig(igpq,jspn))
        z(igpq)=wfir2(ifg)
      end do
      do i=1,nstq
        ist=idxq(i)
! compute inner product
        ephmat(ist,jst,l)=ephmat(ist,jst,l)+zdotc(ngpq(jspn),wfgpq(:,ispn,i), &
         1,z,1)
      end do
    end do
  end do
end do
deallocate(wfir1,wfir2,z)
! convert to phonon coordinates
ld=nstsv**2
do i=1,nstq
  ist=idxq(i)
  do j=1,nst
    jst=idx(j)
    x(:)=ephmat(ist,jst,:)
    call zgemv('T',nbph,nbph,zone,a,nbph,x,1,zzero,ephmat(ist,jst,1),ld)
  end do
end do
deallocate(igpig,igpqig)
deallocate(wfmt,wfgp,wfmtq,wfgpq)
return
end subroutine

