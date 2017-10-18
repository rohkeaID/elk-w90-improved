
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hmlxbsek(ik2)
use modmain
implicit none
! arguments
integer, intent(in) :: ik2
! local variables
integer ik1,ist1,ist2,jst1,jst2
integer i1,i2,j1,j2,a1,a2,b1,b2
integer is,ias,l
real(8) t0
complex(8) zgq0,z1
! automatic arrays
integer idx(nstsv)
! allocatable arrays
complex(8), allocatable :: wfmt1(:,:,:,:),wfir1(:,:,:)
complex(8), allocatable :: wfmt2(:,:,:,:),wfir2(:,:,:)
complex(8), allocatable :: zrhomt(:,:),zrhoir(:)
complex(8), allocatable :: zvclmt(:,:,:),zvclir(:,:)
complex(8), allocatable :: zfmt(:)
! external functions
complex(8) zfinp
external zfinp
! allocate local arrays
allocate(wfmt1(npcmtmax,natmtot,nspinor,nstsv),wfir1(ngtot,nspinor,nstsv))
allocate(wfmt2(npcmtmax,natmtot,nspinor,nstsv),wfir2(ngtot,nspinor,nstsv))
allocate(zrhomt(npcmtmax,natmtot),zrhoir(ngtot))
allocate(zvclmt(npcmtmax,natmtot,nvcbse),zvclir(ngtot,nvcbse))
allocate(zfmt(npcmtmax))
! index to all states
do ist1=1,nstsv
  idx(ist1)=ist1
end do
! calculate the wavefunctions for all states of k-point ik2
call genwfsvp(.false.,.false.,nstsv,idx,vkl(:,ik2),wfmt2,ngtot,wfir2)
l=0
do i2=1,nvbse
  ist2=istbse(i2,ik2)
  do j2=1,ncbse
    jst2=jstbse(j2,ik2)
    a2=ijkbse(i2,j2,ik2)
    l=l+1
! calculate the complex overlap density
    call genzrho(.true.,.true.,wfmt2(:,:,:,ist2),wfir2(:,:,ist2), &
     wfmt2(:,:,:,jst2),wfir2(:,:,jst2),zrhomt,zrhoir)
! compute the Coulomb potential
    call genzvclmt(nrcmt,nrcmti,nrcmtmax,rcmt,npcmtmax,zrhomt,zvclmt(:,:,l))
    call zpotcoul(nrcmt,nrcmti,npcmt,npcmti,nrcmtmax,rcmt,ng2gk,1,gc,ngvec, &
     jlgrmt,ylmg,sfacg,zrhoir,npcmtmax,zvclmt(:,:,l),zvclir(:,l),zgq0)
  end do
end do
t0=occmax*wkptnr
! start loop over ik1
do ik1=1,nkptnr
  if (ik1.eq.ik2) then
    wfmt1(:,:,:,:)=wfmt2(:,:,:,:)
    wfir1(:,:,:)=wfir2(:,:,:)
  else
    call genwfsvp(.false.,.false.,nstsv,idx,vkl(:,ik1),wfmt1,ngtot,wfir1)
  end if
  do i1=1,nvbse
    ist1=istbse(i1,ik1)
    do j1=1,ncbse
      jst1=jstbse(j1,ik1)
      a1=ijkbse(i1,j1,ik1)
! calculate the complex overlap density
      call genzrho(.true.,.true.,wfmt1(:,:,:,ist1),wfir1(:,:,ist1), &
       wfmt1(:,:,:,jst1),wfir1(:,:,jst1),zrhomt,zrhoir)
      l=0
      do i2=1,nvbse
        ist2=istbse(i2,ik2)
        do j2=1,ncbse
          jst2=jstbse(j2,ik2)
          a2=ijkbse(i2,j2,ik2)
          l=l+1
! compute the matrix element
          z1=t0*zfinp(zrhomt,zrhoir,zvclmt(:,:,l),zvclir(:,l))
          hmlbse(a1,a2)=hmlbse(a1,a2)+z1
! compute off-diagonal blocks if required
          if (bsefull) then
            b1=a1+nbbse
            b2=a2+nbbse
            hmlbse(b1,b2)=hmlbse(b1,b2)-conjg(z1)
! conjugate the potential
            do ias=1,natmtot
              is=idxis(ias)
              call zfmtconj(nrcmt(is),nrcmti(is),npcmt(is),zvclmt(:,ias,l))
            end do
            zvclir(:,l)=conjg(zvclir(:,l))
            z1=t0*zfinp(zrhomt,zrhoir,zvclmt(:,:,l),zvclir(:,l))
            hmlbse(a1,b2)=hmlbse(a1,b2)+z1
            hmlbse(b1,a2)=hmlbse(b1,a2)-conjg(z1)
          end if
        end do
      end do
    end do
  end do
end do
deallocate(wfmt1,wfmt2,wfir1,wfir2)
deallocate(zrhomt,zrhoir,zvclmt,zvclir,zfmt)
return
end subroutine
