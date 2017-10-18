
! Copyright (C) 2017 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genjlgqr(gqc,jlgqr)
use modmain
implicit none
! arguments
real(8), intent(in) :: gqc(ngrf)
real(8), intent(out) :: jlgqr(njcmax,nspecies,ngrf)
! local variables
integer ig,is,lmax
integer nrci,irc,i
real(8) t1,t2
! generate spherical Bessel functions on the coarse radial mesh over all species
do ig=1,ngrf
  t1=gqc(ig)
  do is=1,nspecies
    nrci=nrcmti(is)
    lmax=lmaxi
    i=1
    do irc=1,nrcmt(is)
      t2=t1*rcmt(irc,is)
      call sbessel(lmax,t2,jlgqr(i,is,ig))
      i=i+lmax+1
      if (irc.eq.nrci) lmax=lmaxo
    end do
  end do
end do
return
end subroutine

