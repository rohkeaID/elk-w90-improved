
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rotdmat(rspl,rspn,lmax,nspinor,ld,dmat1,dmat2)
implicit none
! arguments
real(8), intent(in) :: rspl(3,3),rspn(3,3)
integer, intent(in) :: lmax,nspinor
integer, intent(in) :: ld
complex(8), intent(in) :: dmat1(ld,nspinor,ld,nspinor)
complex(8), intent(inout) :: dmat2(ld,nspinor,ld,nspinor)
! local variables
integer lds,ispn,jspn,p
integer lmmax,l,lm1,lm2,nm
real(8), parameter :: eps=1.d-8
real(8) ang(3),angi(3)
real(8) rot(3,3),det,v(3),th
complex(8), parameter :: zzero=(0.d0,0.d0),zone=(1.d0,0.d0)
complex(8) su2(2,2),a(2,2),b(2,2)
! allocatable arrays
complex(8), allocatable :: dm(:,:,:,:),c(:,:),d(:,:)
! external functions
real(8) r3mdet
external r3mdet
lmmax=(lmax+1)**2
allocate(dm(ld,nspinor,ld,nspinor))
allocate(c(lmmax,lmmax),d(lmmax,lmmax))
! find the determinant of the spatial rotation matrix
det=r3mdet(rspl)
if (det.gt.0.d0) then
  p=1
else
  p=-1
end if
! make the rotation matrix proper
rot(:,:)=dble(p)*rspl(:,:)
! compute the Euler angles of the spatial rotation
call roteuler(rot,ang)
! inverse rotation: the matrix is to be rotated, not the spherical harmonics
angi(1)=-ang(3)
angi(2)=-ang(2)
angi(3)=-ang(1)
! determine the rotation matrix for complex spherical harmonics
call ylmrot(p,angi,lmax,lmmax,d)
! apply (l,m) rotation matrix as U*D*conjg(U')
lds=ld*nspinor
do ispn=1,nspinor
  do jspn=1,nspinor
    lm1=1
    do l=0,lmax
      nm=2*l+1
      call zgemm('N','N',nm,lmmax,nm,zone,d(lm1,lm1),lmmax, &
       dmat1(lm1,ispn,1,jspn),lds,zzero,c(lm1,1),lmmax)
      lm1=lm1+nm
    end do
    lm1=1
    do l=0,lmax
      nm=2*l+1
      call zgemm('N','C',lmmax,nm,nm,zone,c(1,lm1),lmmax,d(lm1,lm1),lmmax, &
       zzero,dm(1,ispn,lm1,jspn),lds)
      lm1=lm1+nm
    end do
  end do
end do
! spin rotation if required
if (nspinor.eq.2) then
! convert spin rotation matrix to axis-angle form
  call rotaxang(eps,rspn,det,v,th)
! find the SU(2) representation of the rotation matrix
  call axangsu2(v,th,su2)
! apply SU(2) symmetry matrix as U*D*conjg(U*) andd add to dmat2
  do lm1=1,lmmax
    do lm2=1,lmmax
      a(:,:)=dm(lm1,:,lm2,:)
      call z2mm(su2,a,b)
      call z2mmct(b,su2,a)
      dmat2(lm1,:,lm2,:)=dmat2(lm1,:,lm2,:)+a(:,:)
    end do
  end do
else
  dmat2(1:lmmax,1,1:lmmax,1)=dmat2(1:lmmax,1,1:lmmax,1)+dm(1:lmmax,1,1:lmmax,1)
end if
deallocate(dm,c,d)
return
end subroutine

