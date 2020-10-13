
! Copyright (C) 2015 Manh Duc Le, 2017-18 Arsenii Gerasimov, Yaroslav Kvashnin
! and Lars Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genw90twf
! !INTERFACE:
subroutine genw90twf(twfmt,twfir)
! !USES:
use modmain
use modw90
! !INPUT/OUTPUT PARAMETERS:
!   twfmt : muffin-tin trial wavefunction in spherical coordinates
!           (inout,complex(npcmtmax,natmtot,nspinor,wann_nproj))
!   twfir : interstitial trial wavefunction in real-space
!           (inout,complex(ngtot,nspinor,wann_nproj))
! !DESCRIPTION:
!   Generates the trial orbitals and trial wavefunctions.
!
! !REVISION HISTORY:
!   Created January 2015 (Manh Duc Le),
!   Modified July 2018 (Yaroslav Kvashnin, Arsenii Gerasimov)
!EOP
!BOC
implicit none
complex(8), intent(inout) :: twfmt(npcmtmax,natmtot,nspinor,wann_nproj)
complex(8), intent(inout) :: twfir(ngtot,nspinor,wann_nproj)
! local variables
integer    i,is,ia,ias,lmaxproj,lmmaxproj,lm,l,m,n,nr,ir,irc
integer    irco,nrc,nrci
real(8)    s2,s3,s6,s12
complex(8) z1
! allocatable arrays
real(8),    allocatable :: rlm(:)
real(8),    allocatable :: wann_proj_yaxis(:) !AG
real(8),    allocatable :: wann_proj_rot(:,:) !AG
complex(8), allocatable :: wann_projclm_rot(:) !AG
complex(8), allocatable :: twfmt1(:)
! automatic arrays
logical noweight(wann_nproj)
integer irlm(16)
real(8) v1(3)
! fixed values
data irlm / 1, &            ! s
            2,3,1, &        ! pz,px,py -> py,pz,px
            3,4,2,5,1, &    ! z2,xz,yz,xy2y,xy -> xy,yz,z2,zx,x2y2
            4,5,3,6,2,7,1 / !   z3,xz2,yz2,zx2y2,xyz,xx2y2,yx2y2 ->
                            !-> yx2y2,xyz,yz2,z3,xz2,zx2y2,xx2y2
s2  = 1./sqrt(2.)
s3  = 1./sqrt(3.)
s6  = 1./sqrt(6.)
s12 = 1./sqrt(12.)
!-------------------------------------------------------------------------------

! Determine the maximum lmax for projections (<=3)
lmaxproj=maxval(wann_proj_l)
if( any( wann_proj_l.lt. 0 ) .and. lmaxproj.lt.1 ) lmaxproj=1 ! sp,sp2,sp3 orbitals
if( any( wann_proj_l.lt.-3 ) .and. lmaxproj.lt.2 ) lmaxproj=2 ! sp3d,sp3d2 orbitals
if( lmaxproj .gt. 3 ) then
  write(*,*)
  write(*,'("error(genw90twf): projections with l>3 not supported")')
  write(*,*)
  stop
end if
lmmaxproj = lmaxproj*( lmaxproj + 2 ) + 1

allocate(rlm(lmmaxproj))
allocate(wann_projulr(nrmtmax,0:lmaxo,natmtot,wann_nproj))
allocate(wann_projclm(lmaxo*( lmaxo + 2 ) + 1,natmtot,wann_nproj))
allocate(wann_proj_haswt(wann_nproj,natmtot))

! At the moment, non-atom centred projections are treated as random (!)
noweight = .true.
wann_projulr = 0.d0
do is = 1,nspecies
  nr = nrcmt(is)
  do ia = 1,natoms(is)
    ias = idxas(ia,is)
    ! Code taken from writegeom.f90
    if (molecule) then
      ! Map lattice coordinates to [-0.5,0.5)
      v1(:) = atposl(:,ia,is)
      do i = 1,3
        if ( v1(i) .gt. 0.5d0 ) v1(i) = v1(i) - 1.d0
      end do
    else
      ! Otherwise write lattice coordinates
      v1(:) = atposl(:,ia,is)
    end if
    wann_proj_haswt(:,ias) = .false.
    do n = 1,wann_nproj
      ! Only treats atom centred projections in this version (!)
      if( abs( wann_proj_site(1,n) - v1(1) ) .lt. 0.0001 .and. &
         &abs( wann_proj_site(2,n) - v1(2) ) .lt. 0.0001 .and. &
         &abs( wann_proj_site(3,n) - v1(3) ) .lt. 0.0001) then
        !write(*,*) 'Atoms positions from genw90twf: v1 first and wann_proj_site next'
        !write(*,*) v1
        !write(*,*) wann_proj_site(:,n)
        wann_proj_haswt(n,ias) = .true.
        noweight(n) = .false.
        do l = 0,lmaxproj
          irc = 1
          do ir = 1,nr
            wann_projulr(ir,l,ias,n) = apwfr(irc,1,1,l,ias) ! Just use the APW radial fn
                                                            ! (on coarse mesh)
            irc = irc + lradstp
          end do
        end do
        rlm = 0
        wann_projclm(:,ias,n) = zzero
        if( wann_proj_l(n) .ge. 0 ) then
          ! Converts projection angular part from Wannier90 indexing (tab 3.1)
                                                                  ! to standard
          l = wann_proj_l(n)
          if( wann_proj_m(n) .gt. (2*l + 1) ) goto 100
          m = irlm(wann_proj_m(n) + l*l) - 1 - l
          lm = l*(l + 1) + m + 1
          rlm(lm) = 1
          !write(*,*)'n,l,m,lm,w90_m = ',n,l,m,lm,wann_proj_m(n)
        else
          select case( wann_proj_l(n) )
            ! Lin. comb. of real spherical harmonics from table 3.2 of Wannier90
                                                                    ! user guide
            case (-1)                            ! sp
              rlm(1) = s2                        ! 1/sqrt(2) s
              select case ( wann_proj_m(n) )
                case (1)
                  rlm(4) = s2                    ! 1/sqrt(2) px
                case (2)
                  rlm(4) = -s2                   ! -1/sqrt(2) px
                case default
                  goto 100
              end select
            case (-2)                            ! sp2
              rlm(1) = s3                        ! 1/sqrt(3) s
              select case ( wann_proj_m(n) )
                case (1)
                  rlm(4) = -s6; rlm(2) =  s2     ! -1/sqrt(6) px + 1/sqrt(2) py
                case (2)
                  rlm(4) = -s6; rlm(2) = -s2     ! -1/sqrt(6) px - 1/sqrt(2) py
                case (3)
                  rlm(4) = 2.*s6;                ! 2/sqrt(6) px
                case default
                  goto 100
              end select
            case (-3)                            ! sp3
              rlm(1) = 1/2.                      ! 1/sqrt(3) s
              select case ( wann_proj_m(n) )
                case (1)    !py,pz,px
                  rlm(2:4) = (/ 1., 1., 1./)/.2  ! (+px+py+pz)/2
                case (2)
                  rlm(2:4) = (/-1.,-1., 1./)/.2  ! (+px-py-pz)/2
                case (3)
                  rlm(2:4) = (/ 1.,-1.,-1./)/.2  ! (-px+py-pz)/2
                case (4)
                  rlm(2:4) = (/-1., 1.,-1./)/.2  ! (-px-py+pz)/2
                case default
                  goto 100
              end select
            case (-4)                            ! sp3d
              select case ( wann_proj_m(n) )
                case (1)    ! s, py,pz,px
                  rlm(1:4) = (/s3, s2,0d0,-s6/)  ! 1/sqrt(3) s - 1/sqrt(6) px +
                                                                ! + 1/sqrt(2) py
                case (2)
                  rlm(1:4) = (/s3,-s2,0d0,-s6/)  ! 1/sqrt(3) s - 1/sqrt(6) px -
                                                                ! - 1/sqrt(2) py
                case (3)
                  rlm(1:4) = (/s3,0d0,0d0,2*s6/) ! 1/sqrt(3) s + 2/sqrt(6) px
                case (4)
                  rlm(3)   =  s2; rlm(7) = s2    ! 1/sqrt(2) pz + 1/sqrt(2) dz2
                case (5)
                  rlm(3)   = -s2; rlm(7) = s2    ! -1/sqrt(2) pz + 1/sqrt(2) dz2
                case default
                  goto 100
              end select
            case (-5)                            ! sp3d2
              rlm(1) = s6                        ! 1/sqrt(6) s
              select case ( wann_proj_m(n) )
                case (1)
                  rlm(4) = -s2                   ! -1/sqrt(2) px
                  rlm(7) = -s12; rlm(9)=1./2.    ! -1/sqrt(12) dz2 + 1/2 dx2-y2
                case (2)
                  rlm(4) =  s2                   ! 1/sqrt(2) px
                  rlm(7) = -s12; rlm(9)=1./2.    ! -1/sqrt(12) dz2 + 1/2 dx2-y2
                case (3)
                  rlm(4) = -s2                   ! -1/sqrt(2) px
                  rlm(7) = -s12; rlm(9)=-1./2.   ! -1/sqrt(12) dz2 + 1/2 dx2-y2
                case (4)
                  rlm(4) =  s2                   ! 1/sqrt(2) px
                  rlm(7) = -s12; rlm(9)=-1./2.   ! -1/sqrt(12) dz2 + 1/2 dx2-y2
                case (5)
                  rlm(3) = -s2;  rlm(7)=s3       ! -1/sqrt(2) pz + 1/sqrt(3) dz2
                case (6)
                  rlm(3) =  s2;  rlm(7)=s3       ! 1/sqrt(2) pz + 1/sqrt(3) dz2
                case default
                  goto 100
              end select
          end select
        end if
        ! Converts from real to complex spherical harmonics
        call rtozflm(lmaxproj,rlm,wann_projclm(:,ias,n))
        !write(*,*) 'ias,n,wann_projclm: ',ias,n,wann_projclm(:,ias,n)
      end if
    end do
  end do ! Loop over atoms
end do ! Loop over species

! Create the trial wavefunction matrix
allocate(twfmt1(npcmtmax))
allocate(wann_proj_yaxis(3)) !AG
allocate(wann_proj_rot(3,3)) !AG
allocate(wann_projclm_rot(lmaxo*( lmaxo + 2 ) + 1)) !AG
twfmt = zzero
twfir = zzero ! Trial wavefunction is zero outside muffin-tin.
!do ispn=1,nspinor !yk
  do is = 1,nspecies
    nrc = nrcmt(is)
    nrci = nrcmti(is)
    irco = nrci + 1
    do ia = 1,natoms(is)
      ias = idxas(ia,is)
      do n = 1,wann_nproj
        if( .not.wann_proj_haswt(n,ias) ) cycle
        ! Zero the wavefunction
        twfmt1 = zzero
        call r3cross(wann_proj_zaxis(:,n),wann_proj_xaxis(:,n),wann_proj_yaxis(:)) !AG
        wann_proj_rot(:,1) = wann_proj_xaxis(:,n) !AG
        wann_proj_rot(:,2) = wann_proj_yaxis(:)   !AG
        wann_proj_rot(:,3) = wann_proj_zaxis(:,n) !AG
        call rotzflm(wann_proj_rot(:,:),&
                                        0,&
                                    lmaxo,&
                                        1,&
                        lmaxo*(lmaxo+2)+1,&
                    wann_projclm(:,ias,n),&
                    wann_projclm_rot(:)) !AG
        i = 1
        do ir = 1,nrci
          lm = 0
          do l = 0,lmaxi
            do m = -l,l
              lm = lm+1
              z1 = wann_projclm_rot(lm)
              if( abs(z1).gt.1.d-14 ) then
                twfmt1(i) = twfmt1(i) + z1*wann_projulr(ir,l,ias,n)
              end if
              i = i+1
            end do
          end do
        end do
        do ir = nrci+1,nrc
          lm = 0
          do l = 0,lmaxo
            do m = -l,l
              lm = lm+1
              z1 = wann_projclm_rot(lm)
              if( abs(z1) .gt. 1.d-14 ) then
                twfmt1(i) = twfmt1(i) + z1*wann_projulr(ir,l,ias,n)
              end if
              i = i + 1
            end do
          end do
        end do
        if ( nspinor .eq. 2 ) then
          call zbsht(nrc,nrci,twfmt1,twfmt(:,ias,&
                                          int(1.5d0-(wann_proj_spin(n)/2d0)),n))
        else
          call zbsht(nrc,nrci,twfmt1,twfmt(:,ias,1,n))
        end if
      end do !n; projections
    end do ! ia; atoms
  end do ! is; species
!end do
deallocate(wann_projclm_rot,wann_proj_rot,wann_proj_yaxis,twfmt1)

deallocate(rlm)

return

100 write(*,*)
    write(*,'("error(genw90twf): (l=",I2,",m=",I2,") combination not recognised.")')&
                                                        wann_proj_l(n),wann_proj_m(n)
    write(*,*)
    stop

end subroutine
!EOC
