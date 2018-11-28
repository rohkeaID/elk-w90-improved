
! Copyright (C) 2017-18 Arsenii Gerasimov, Yaroslav Kvashnin and Lars Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: zfplotUNK
! !INTERFACE:
subroutine zfplotUNK(np,vpl,zfmt,zfir,fp)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   fnum : plot file number (in,integer)
!   nf   : number of functions (in,integer)
!   zfmt : complex muffin-tin function (in,complex(npmtmax,natmtot,nf))
!   zfir : complex intersitial function (in,complex(ngtot,nf))
! !DESCRIPTION:
!
!
! !REVISION HISTORY:
!   Created September 2018 (Arsenii Gerasimov)
!EOP
!BOC
implicit none
! arguments
integer,    intent(in   ) :: np
real(8),    intent(in   ) :: vpl(3,np)
complex(8), intent(in   ) :: zfmt(npcmtmax,natmtot),zfir(ngtot)
complex(8), intent(  out) :: fp(np)
! local variables
integer ias,is,ip
!complex(8) fp_
! allocatable arrays
complex(8), allocatable :: zfmt1(:,:,:)
complex(8), allocatable :: zfft(:)
!-------------------------------------------------------------------------------
! Unpack the muffin-tin function
allocate(zfmt1(lmmaxo,nrcmtmax,natmtot))

do ias=1,natmtot
  is=idxis(ias)
  call zfmtpack(.false.,nrcmt(is),nrcmti(is),zfmt(:,ias),zfmt1(:,:,ias))
end do

! Fourier transform zfir to G-space
allocate(zfft(ngtot))
zfft(:)=zfir(:)
call zfftifc(3,ngridg,-1,zfft)
! Begin loop over all points
do ip=1,np
  call zfip(ip)
end do

deallocate(zfmt1,zfft)
return

contains

subroutine zfip(ip)
implicit none
! arguments
integer,    intent(in   ) :: ip
! local variables
integer is,ia,ias,nrc,nrci
integer ir0,ir,lmax,l,m,lm
integer ig,ifg,i1,i2,i3,i,j
real(8) r,t1
complex(8) t2
complex(8) sum
! automatic arrays
real(8) v1(3),v2(3),v3(3),v4(3),v5(3),tp(2)
complex(8) ylm(lmmaxo),ya(4)
!-------------------------------------------------------------------------------
v2(:) = vpl(:,ip)
call r3frac(epslat,v2)
! Convert point to Cartesian coordinates
call r3mv(avec,v2,v1)
! Check if point is in a muffin-tin
do is = 1,nspecies
  nrc = nrcmt(is)
  nrci = nrcmti(is)
  do ia = 1,natoms(is)
    ias = idxas(ia,is)
    v2(:) = v1(:) - atposc(:,ia,is)
    do i1 = -1,1
      v3(:) = v2(:) + dble(i1)*avec(:,1)
      do i2 = -1,1
        v4(:) = v3(:) + dble(i2)*avec(:,2)
        do i3 = -1,1
          v5(:) = v4(:) + dble(i3)*avec(:,3)
          t1 = v5(1)**2 + v5(2)**2 + v5(3)**2
          if ( t1 .lt. r2cmt(nrci,is) ) then
            call sphcrd(v5,r,tp)
            call genylm(lmaxo,tp,ylm)
            do ir = 1,nrc
              if ( rsp(ir,is) .ge. r ) then
                if ( ir .le. 2 ) then
                  ir0 = 1
                else if ( ir .gt. nrcmt(is)-2 ) then
                  ir0 = nrcmt(is) - 3
                else
                  ir0 = ir - 2
                end if
                r = max(r,rsp(1,is))
                if ( ir0 .le. nrci ) then
                  lmax = lmaxi
                else
                  lmax = lmaxo
                end if
                sum = dcmplx(0.d0,0.d0)
                do l = 0,lmax
                  do m = -l,l
                    lm = idxlm(l,m)
                    do j = 1,4
                      i = ir0 + j - 1
                      ya(j) = zfmt1(lm,i,ias)
                    end do
                    t2 = poly4(rsp(ir0,is),ya,r)
                    sum = sum + t2*ylm(lm)
                  end do
                end do

                goto 10
              end if
            end do
          end if
        end do
      end do
    end do
  end do
end do
! Otherwise use direct Fourier transform of interstitial function
sum = dcmplx(0.d0,0.d0)
do ig = 1,ngvec
  ifg = igfft(ig)
  t1 = vgc(1,ig)*v1(1) + vgc(2,ig)*v1(2) + vgc(3,ig)*v1(3)
  sum = sum + zfft(ifg)*cmplx( cos(t1),sin(t1),8 )
end do
10 continue
fp(ip) = sum
return
end subroutine

complex(8) function poly4(xa,ya,x)
implicit none
! arguments
real(8),    intent(in) :: xa(4),x
complex(8), intent(in) :: ya(4)
! local variables
real(8)    x0,x1,x2,x3
complex(8) y0,y1,y2,y3
complex(8) c1,c2,c3,t0,t1,t2,t3,t4,t5,t6
!-------------------------------------------------------------------------------
! evaluate the polynomial coefficient
x0 = xa(1)
x1 = xa(2) - x0; x2 = xa(3) - x0; x3 = xa(4) - x0
y0 = ya(1)
y1 = ya(2) - y0; y2 = ya(3) - y0; y3 = ya(4) - y0
t1 = x1*x2*y3;   t2 = x2*x3*y1;   t3 = x3*x1*y2
t4 = x1 - x2;    t5 = x1 - x3;    t6 = x2 - x3
t0 = 1.d0/( x1*x2*x3*t4*t5*t6 )
c3 = t1*t4 + t2*t6 - t3*t5
t4 = x1**2; t5 = x2**2; t6 = x3**2
c2 = t1*( t5 - t4 ) + t2*( t6 - t5 ) + t3*( t4 - t6 )
c1 = t1*( x2*t4 - x1*t5 ) + t2*( x3*t5 - x2*t6 ) + t3*( x1*t6 - x3*t4 )
t1 = x - x0
! evaluate the polynomial
poly4 = y0 + t0*t1*( c1 + t1*( c2 + c3*t1 ) )
return
end function

end subroutine
