
! Copyright (C) 2016 A. Davydov, A. Sanna, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine dysonr(ik,swfm,sf)
use modmain
use modgw
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: swfm(nstsv,nstsv,0:nwfm)
real(8), intent(out) :: sf(nwplot)
! local variables
integer ist,jst,iw,jw,info
real(8) dw,w,e,sum,t1
complex(8) z1
! allocatable arrays
integer, allocatable :: ipiv(:)
complex(8), allocatable :: wfm(:),ufm(:),wr(:),ur(:)
complex(8), allocatable :: swr(:,:,:),gs(:),g(:,:),work(:)
allocate(wfm(0:nwfm),ufm(0:nwfm),wr(nwplot),ur(nwplot))
allocate(swr(nstsv,nstsv,nwplot))
! complex fermionic frequencies
do iw=-nwfm,nwfm,2
  jw=(iw+nwfm)/2
  wfm(jw)=cmplx(0.d0,wgw(iw),8)
end do
! complex array of real axis frequencies
dw=(wplot(2)-wplot(1))/dble(nwplot)
do iw=1,nwplot
  w=dw*dble(iw-1)+wplot(1)
  wr(iw)=cmplx(w,0.d0,8)
end do
do ist=1,nstsv
  do jst=1,nstsv
    ufm(:)=swfm(ist,jst,:)
    call pade(nwfm+1,wfm,ufm,nwplot,wr,ur)
! store the real axis Sigma
    swr(ist,jst,:)=ur(:)
  end do
end do
deallocate(wfm,ufm,ur)
! solve the Dyson equation for each frequency
allocate(gs(nstsv),g(nstsv,nstsv))
allocate(ipiv(nstsv),work(nstsv))
do iw=1,nwplot
  w=dble(wr(iw))
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
end do
deallocate(ipiv,work)
deallocate(wr,swr,gs,g)
return
end subroutine

