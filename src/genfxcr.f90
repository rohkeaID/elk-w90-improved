
! Copyright (C) 2012 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genfxcr(tsh,fxcmt,fxcir)
use modmain
use modtddft
use modfxcifc
implicit none
! arguments
logical, intent(in) :: tsh
real(8), intent(out) :: fxcmt(lmmaxvr,nrmtmax,natmtot)
real(8), intent(out) :: fxcir(ngtot)
! local variables
integer idm,is,ia,ias
integer nr,nri,ir,i,n
real(8) t1
real(8), allocatable :: rho(:),rhoup(:),rhodn(:),mag(:,:)
real(8), allocatable :: fxc(:),fxcuu(:),fxcud(:),fxcdd(:)
! number of independent spin components
n=lmmaxvr*nrmtmax
allocate(rho(n),fxc(n))
if (spinpol) then
  allocate(mag(n,3))
  n=max(n,ngtot)
  allocate(rhoup(n),rhodn(n))
  allocate(fxcuu(n),fxcud(n),fxcdd(n))
end if
!---------------------------!
!     muffin-tin kernel     !
!---------------------------!
do is=1,nspecies
  nr=nrmt(is)
  nri=nrmtinr(is)
  n=lmmaxvr*nr
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! compute the density in spherical coordinates
    call rbsht(nr,nri,1,rhomt(:,:,ias),1,rho)
    if (spinpol) then
!------------------------!
!     spin-polarised     !
!------------------------!
! magnetisation in spherical coordinates
      do idm=1,ndmag
        call rbsht(nr,nri,1,magmt(:,:,ias,idm),1,mag(:,idm))
      end do
      if (ncmag) then
! non-collinear (use Kubler's trick)
        do i=1,n
! compute rhoup=(rho+|m|)/2 and rhodn=(rho-|m|)/2
          t1=sqrt(mag(i,1)**2+mag(i,2)**2+mag(i,3)**2)
          rhoup(i)=0.5d0*(rho(i)+t1)
          rhodn(i)=0.5d0*(rho(i)-t1)
        end do
      else
! collinear
        do i=1,n
! compute rhoup=(rho+m_z)/2 and rhodn=(rho-m_z)/2
          rhoup(i)=0.5d0*(rho(i)+mag(i,1))
          rhodn(i)=0.5d0*(rho(i)-mag(i,1))
        end do
      end if
! compute fxc
      call fxcifc(fxctype,n=n,rhoup=rhoup,rhodn=rhodn,fxcuu=fxcuu,fxcud=fxcud, &
       fxcdd=fxcdd)
! form the scalar quantity dv/drho
      do i=1,n
        fxc(i)=0.25d0*(fxcuu(i)+2.d0*fxcud(i)+fxcdd(i))
      end do
    else
!--------------------------!
!     spin-unpolarised     !
!--------------------------!
      call fxcifc(fxctype,n=n,rho=rho,fxc=fxc)
    end if
    if (tsh) then
! convert fxc to spherical harmonics if required
      call rfsht(nr,nri,1,fxc,1,fxcmt(:,:,ias))
    else
      call dcopy(n,fxc,1,fxcmt(:,:,ias),1)
    end if
  end do
end do
!-----------------------------!
!     interstitial kernel     !
!-----------------------------!
if (spinpol) then
!------------------------!
!     spin-polarised     !
!------------------------!
  if (ncmag) then
! non-collinear
    do ir=1,ngtot
      t1=sqrt(magir(ir,1)**2+magir(ir,2)**2+magir(ir,3)**2)
      rhoup(ir)=0.5d0*(rhoir(ir)+t1)
      rhodn(ir)=0.5d0*(rhoir(ir)-t1)
    end do
  else
! collinear
    do ir=1,ngtot
      rhoup(ir)=0.5d0*(rhoir(ir)+magir(ir,1))
      rhodn(ir)=0.5d0*(rhoir(ir)-magir(ir,1))
    end do
  end if
! compute fxc
  call fxcifc(fxctype,n=ngtot,rhoup=rhoup,rhodn=rhodn,fxcuu=fxcuu,fxcud=fxcud, &
   fxcdd=fxcdd)
  do ir=1,ngtot
    fxcir(ir)=0.25d0*(fxcuu(ir)+2.d0*fxcud(ir)+fxcdd(ir))
  end do
else
!--------------------------!
!     spin-unpolarised     !
!--------------------------!
  call fxcifc(fxctype,n=ngtot,rho=rhoir,fxc=fxcir)
end if
deallocate(rho,fxc)
if (spinpol) then
  deallocate(mag,rhoup,rhodn)
  deallocate(fxcuu,fxcud,fxcdd)
end if
return
end subroutine

