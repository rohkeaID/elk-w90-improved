
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine grad2rfmt(nr,nri,r,rfmt,g2rfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: r(nr)
real(8), intent(in) :: rfmt(*)
real(8), intent(out) :: g2rfmt(*)
! local variables
integer nro,iro,ir
integer l,m,lm,npi,i
real(8) t1
! automatic arrays
real(8) ri(nr),ri2(nr),fr(nr),cf(3,nr)
! tabulate 1/r and 1/r^2
do ir=1,nr
  ri(ir)=1.d0/r(ir)
  ri2(ir)=ri(ir)**2
end do
nro=nr-nri
iro=nri+1
npi=lmmaxi*nri
lm=0
do l=0,lmaxi
  t1=-dble(l*(l+1))
  do m=-l,l
    lm=lm+1
! use a cubic spline to compute radial derivatives
    i=lm
    do ir=1,nri
      fr(ir)=rfmt(i)
      i=i+lmmaxi
    end do
    do ir=iro,nr
      fr(ir)=rfmt(i)
      i=i+lmmaxo
    end do
    call spline(nr,r,fr,cf)
! apply Laplacian
    i=lm
    do ir=1,nri
      g2rfmt(i)=2.d0*(ri(ir)*cf(1,ir)+cf(2,ir))+ri2(ir)*t1*rfmt(i)
      i=i+lmmaxi
    end do
    do ir=iro,nr
      g2rfmt(i)=2.d0*(ri(ir)*cf(1,ir)+cf(2,ir))+ri2(ir)*t1*rfmt(i)
      i=i+lmmaxo
    end do
  end do
end do
do l=lmaxi+1,lmaxo
  t1=-dble(l*(l+1))
  do m=-l,l
    lm=lm+1
    i=npi+lm
    do ir=iro,nr
      fr(ir)=rfmt(i)
      i=i+lmmaxo
    end do
    call spline(nro,r(iro),fr(iro),cf(1,iro))
    i=npi+lm
    do ir=iro,nr
      g2rfmt(i)=2.d0*(ri(ir)*cf(1,ir)+cf(2,ir))+ri2(ir)*t1*rfmt(i)
      i=i+lmmaxo
    end do
  end do
end do
! apply smoothing if required
call rfmtsm(msmooth,nr,nri,g2rfmt)
return
end subroutine

