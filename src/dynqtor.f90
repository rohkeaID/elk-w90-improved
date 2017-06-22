
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dynqtor(dynq,dynr)
use modmain
use modphonon
implicit none
! arguments
complex(8), intent(in) :: dynq(nbph,nbph,nqpt)
complex(8), intent(out) :: dynr(nbph,nbph,nqptnr)
! local variables
integer ir,iq,i,j,n
integer isym,lspl
integer i1,i2,i3,j1,j2,j3
real(8) v1(3),v2(3),v3(3)
real(8) s(3,3),t1
complex(8) z1
! allocatable arrays
complex(8), allocatable :: dyns(:,:)
allocate(dyns(nbph,nbph))
dynr(:,:,:)=0.d0
! loop over q-vectors
do j1=0,ngridq(1)-1
  v1(1)=dble(j1)/dble(ngridq(1))
  do j2=0,ngridq(2)-1
    v1(2)=dble(j2)/dble(ngridq(2))
    do j3=0,ngridq(3)-1
      v1(3)=dble(j3)/dble(ngridq(3))
      iq=iqmap(j1,j2,j3)
! map v1 to the first Brillouin zone
      v2(:)=v1(:)
      call vecfbz(epslat,bvec,v2)
! rotate and add the dynamical matrix of the reduced q-point with all symmetries
      n=0
      dyns(:,:)=0.d0
      do isym=1,nsymcrys
        lspl=lsplsymc(isym)
        s(:,:)=dble(symlat(:,:,lspl))
        call r3mtv(s,vql(:,iq),v3)
        call vecfbz(epslat,bvec,v3)
        t1=abs(v2(1)-v3(1))+abs(v2(2)-v3(2))+abs(v2(3)-v3(3))
        if (t1.lt.epslat) then
          call dynsymapp(isym,vql(:,iq),dynq(:,:,iq),dyns)
          n=n+1
        end if
      end do
      if (n.eq.0) then
        write(*,*)
        write(*,'("Error(dynqtor): vector ",3G18.10)') v1
        write(*,'(" cannot be mapped to reduced q-point set")')
        write(*,*)
        stop
      end if
      t1=1.d0/dble(n)
      dyns(:,:)=t1*dyns(:,:)
! loop over R-vectors
      ir=0
      do i3=ngridq(3)/2-ngridq(3)+1,ngridq(3)/2
        do i2=ngridq(2)/2-ngridq(2)+1,ngridq(2)/2
          do i1=ngridq(1)/2-ngridq(1)+1,ngridq(1)/2
            ir=ir+1
            t1=twopi*(v1(1)*dble(i1)+v1(2)*dble(i2)+v1(3)*dble(i3))
            z1=cmplx(cos(t1),sin(t1),8)
            do i=1,nbph
              do j=1,nbph
                dynr(i,j,ir)=dynr(i,j,ir)+z1*dyns(i,j)
              end do
            end do
          end do
        end do
      end do
    end do
  end do
end do
t1=1.d0/dble(nqptnr)
dynr(:,:,:)=t1*dynr(:,:,:)
deallocate(dyns)
return
end subroutine

