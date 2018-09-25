
! Copyright (C) 2016 A. Davydov, A. Sanna, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine dysonr(ik,wr,swfm,sf)
use modmain
use modgw
use modomp
implicit none
! arguments
integer, intent(in) :: ik
real(8), intent(in) :: wr(nwplot)
complex(8), intent(in) :: swfm(nstsv,nstsv,0:nwfm)
real(8), intent(out) :: sf(nwplot)
! local variables
integer ist,jst,iw
integer info,nthd
real(8) w,e,sum,t1
complex(8) z1
! allocatable arrays
integer, allocatable :: ipiv(:)
complex(8), allocatable :: swr(:,:,:),gs(:),g(:,:),work(:)
allocate(swr(nstsv,nstsv,nwplot))
swr(:,:,:)=0.d0
call omp_hold(nstsv,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ist) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do jst=1,nstsv
  do ist=1,nstsv
    if ((gwdiag.gt.0).and.(ist.ne.jst)) cycle
! perform analytic continuation from the imaginary to the real axis
    call acgwse(ist,jst,swfm,wr,swr)
  end do
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
! solve the Dyson equation for each frequency
call omp_hold(nwplot,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(gs,g,ipiv,work) &
!$OMP PRIVATE(w,ist,jst,e,t1) &
!$OMP PRIVATE(z1,info,sum) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do iw=1,nwplot
  allocate(gs(nstsv),g(nstsv,nstsv))
  allocate(ipiv(nstsv),work(nstsv))
  w=wr(iw)
! compute the diagonal matrix G_s
  do ist=1,nstsv
    e=evalsv(ist,ik)-efermi
    t1=sign(swidth,e)
    gs(ist)=1.d0/cmplx(w-e,t1,8)
  end do
! compute 1 - G_s Sigma
  do ist=1,nstsv
    z1=-gs(ist)
    g(ist,:)=z1*swr(ist,:,iw)
    g(ist,ist)=g(ist,ist)+1.d0
  end do
! invert this matrix
  call zgetrf(nstsv,nstsv,g,nstsv,ipiv,info)
  if (info.eq.0) call zgetri(nstsv,g,nstsv,ipiv,work,nstsv,info)
  if (info.ne.0) then
    write(*,*)
    write(*,'("Error(dysonr): unable to solve the Dyson equation")')
    write(*,'(" for frequency ",I6)') iw
    write(*,*)
    stop
  end if
! compute G = (1 - G_s Sigma)^(-1) G_s
  do jst=1,nstsv
    z1=gs(jst)
    g(:,jst)=z1*g(:,jst)
  end do
! determine the spectral function
  sum=0.d0
  do ist=1,nstsv
    sum=sum+abs(aimag(g(ist,ist)))
  end do
  sf(iw)=sum*occmax/pi
  deallocate(gs,g,ipiv,work)
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
deallocate(swr)
return
end subroutine

