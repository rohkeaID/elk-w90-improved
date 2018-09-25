
! Copyright (C) 2018 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rhomagq
use modmain
use modulr
use modomp
implicit none
! local variables
integer idm,is,ias,ir
integer npc,i,nthd
! allocatable arrays
complex(8), allocatable :: zfft(:)
! partial Fourier transform of muffin-tin density
do ias=1,natmtot
  is=idxis(ias)
  npc=npcmt(is)
  call omp_hold(npc,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zfft) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
  do i=1,npc
    allocate(zfft(nqpt))
    zfft(:)=rhormt(i,ias,:)
    call zfftifc(3,ngridq,-1,zfft)
    rhoqmt(i,ias,:)=zfft(:)
    deallocate(zfft)
  end do
!$OMP END DO
!$OMP END PARALLEL
  call omp_free(nthd)
end do
! partial Fourier transform of interstitial density
call omp_hold(ngtot,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zfft) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ir=1,ngtot
  allocate(zfft(nqpt))
  zfft(:)=rhorir(ir,:)
  call zfftifc(3,ngridq,-1,zfft)
  rhoqir(ir,:)=zfft(:)
  deallocate(zfft)
end do
!$OMP END DO
!$OMP END PARALLEL
call omp_free(nthd)
! partial Fourier transform of the muffin-tin magnetisation
do idm=1,ndmag
  do ias=1,natmtot
    is=idxis(ias)
    npc=npcmt(is)
    call omp_hold(npc,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zfft) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
    do i=1,npc
      allocate(zfft(nqpt))
      zfft(:)=magrmt(i,ias,idm,:)
      call zfftifc(3,ngridq,-1,zfft)
      magqmt(i,ias,idm,:)=dble(zfft(:))
      deallocate(zfft)
    end do
!$OMP END DO
!$OMP END PARALLEL
    call omp_free(nthd)
  end do
end do
! partial Fourier transform of the interstitial magnetisation
do idm=1,ndmag
  call omp_hold(ngtot,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zfft) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
  do ir=1,ngtot
    allocate(zfft(nqpt))
    zfft(:)=magrir(ir,idm,:)
    call zfftifc(3,ngridq,-1,zfft)
    magqir(ir,idm,:)=zfft(:)
    deallocate(zfft)
  end do
!$OMP END DO
!$OMP END PARALLEL
  call omp_free(nthd)
end do
return
end subroutine

