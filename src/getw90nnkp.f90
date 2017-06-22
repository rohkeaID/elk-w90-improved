! Copyright (C) 2015 Jon Lafuente and Manh Duc Le
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: getw90nnkp
! !INTERFACE:
subroutine getw90nnkp
! !USES:
use modmain
use modw90
! !DESCRIPTION:
!   Determines the list of neighbouring k-points on which to calculate the
!   Amn and Mmn overlap matrices required by wannier90.
!
! !REVISION HISTORY:
!   Created July 2015 (Jon Lafuente and Manh Duc Le)
!EOP
!BOC
implicit none
! local variables
integer i,j,ik,jk,ig,inkp,id,ib,jb,jd,nshellused,info
real(8) dist0,dist1,delta
logical lpar,lb1
! allocatable arrays
integer, allocatable :: nn(:),shell_list(:)
real(8), allocatable :: dist(:),bvectors(:,:,:)
real(8), allocatable :: amat(:,:),umat(:,:),vmat(:,:),smat(:,:),sing(:)
! automatic arrays
real(8) bvector(3),bweight(6),work(60)
real(8), parameter  :: targ(6)=(/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/)
! Allocates arrays
allocate(nn(wann_maxshell))
allocate(dist(wann_maxshell))
! First cycle to find distances and number of neighbours in each k-point shells,
dist = 9.d9
nn = 0
bweight = 0.d0
dist1=0.d0
do id=1,wann_maxshell
  do ig=1,ngvec
    if (any(abs(ivg(:,ig)).gt.1)) cycle
    do jk=1,nkpt
      bvector = vkc(:,1)-(vkc(:,jk)+vgc(:,ig))
      dist0 = sqrt( dot_product(bvector,bvector) )
      if(dist0.gt.(dist1+wann_tol).and.dist0.le.(dist(id)+wann_tol)) dist(id)=dist0
    end do
  end do
  do ig=1,ngvec
    if (any(abs(ivg(:,ig)).gt.2)) cycle
    do jk=1,nkpt
      bvector = vkc(:,1)-(vkc(:,jk)+vgc(:,ig))
      dist0 = sqrt( dot_product(bvector,bvector) )
      if (abs(dist0-dist(id)).le.wann_tol) nn(id) = nn(id) + 1
    end do
  end do
  dist1 = dist(id)
  if(dist1.eq.9.d9) then
    wann_maxshell=id-1
    exit
  end if
end do
! Loops over the shells, to find all bvectors
allocate(bvectors(3,wann_maxshell,maxval(nn)))
ok: do id=1,wann_maxshell
  ! Finds the bvectors belonging to this shell
  ib = 1
  do ik=1,nkpt
    do ig=1,ngvec
      if (any(abs(ivg(:,ig)).gt.2)) cycle
      do jk=1,nkpt
        bvector = vkc(:,1)-(vkc(:,jk)+vgc(:,ig))
        dist0 = sqrt( dot_product(bvector,bvector) )
        if (abs(dist0-dist(id)).le.wann_tol) then
          bvectors(:,id,ib) = bvector
          ib=ib+1
          if(ib.gt.nn(id)) cycle ok
        end if
      end do
    end do
  end do
end do ok
! Loops over the shells, to calculate if the condition specified by equation B1 of
!   Marzari and Vanderbilt (Phys. Rev. B 56, 12847, 1997) is satisfied.
!   Code here follows that in kmesh.F90 from Wannier90
allocate(shell_list(wann_maxshell))
shell_list = -1
nshellused = 0
do id=1,wann_maxshell
! Check if shell is not parallel to existing shell (cosine=1)
  if(id.gt.1) then
    lpar = .false.
    do ib=1,nn(id)
      do jd=1,id-1
        do jb=1,nn(jd)
          delta=dot_product(bvectors(:,id,ib),bvectors(:,jd,jb)) / &
           sqrt(dot_product(bvectors(:,id,ib),bvectors(:,id,ib)) * &
                dot_product(bvectors(:,jd,jb),bvectors(:,jd,jb)))
          if(abs(abs(delta)-1.d0)<wann_tol) then
            lpar=.true.
            exit
          end if
        end do
        if(lpar) exit
      end do
      if(lpar) exit
    end do
    if(lpar) cycle
  end if
! Calculates the singular value decomposition of the |b^2| values found so far.
  nshellused = nshellused + 1
  shell_list(nshellused) = id
  allocate(amat(6,nshellused));          amat = 0.d0
  allocate(umat(6,6));                   umat = 0.d0
  allocate(vmat(nshellused,nshellused)); vmat = 0.d0
  allocate(smat(nshellused,6));          smat = 0.d0
  allocate(sing(nshellused));            sing = 0.d0
  do jd=1,nshellused
    do jb=1,nn(shell_list(jd))
      amat(1,jd)=amat(1,jd)+bvectors(1,shell_list(jd),jb)*bvectors(1,shell_list(jd),jb)
      amat(2,jd)=amat(2,jd)+bvectors(2,shell_list(jd),jb)*bvectors(2,shell_list(jd),jb)
      amat(3,jd)=amat(3,jd)+bvectors(3,shell_list(jd),jb)*bvectors(3,shell_list(jd),jb)
      amat(4,jd)=amat(4,jd)+bvectors(1,shell_list(jd),jb)*bvectors(2,shell_list(jd),jb)
      amat(5,jd)=amat(5,jd)+bvectors(2,shell_list(jd),jb)*bvectors(3,shell_list(jd),jb)
      amat(6,jd)=amat(6,jd)+bvectors(3,shell_list(jd),jb)*bvectors(1,shell_list(jd),jb)
    end do
  end do
  info=0
  call dgesvd('A','A',6,nshellused,amat,6,sing,umat,6,vmat,nshellused,work,60,info)
  if(info.ne.0) then
    write(*,*)
    write(*,'("Error(getw90nnkp): error computing singular value decomposition")')
    write(*,*)
    stop
  end if
  if(any(abs(sing).lt.wann_tol)) then
    shell_list(nshellused) = -1
    nshellused = nshellused - 1
    deallocate(amat)
    deallocate(umat)
    deallocate(vmat)
    deallocate(smat)
    deallocate(sing)
    cycle
  end if
! Calculates the weights from the SVD
  do jd=1,nshellused
    smat(jd,jd) = 1/sing(jd)
  end do
  bweight(1:nshellused)=matmul(transpose(vmat),matmul(smat,matmul(transpose(umat),targ)))
  lb1=.true.
  do i=1,3
    do j=1,3
      delta = 0.d0
      do jd=1,nshellused
        do jb=1,nn(shell_list(jd))
          delta=delta+bweight(jd)*bvectors(i,shell_list(jd),jb)*bvectors(j,shell_list(jd),jb)
        end do
      end do
      if(i.eq.j.and.abs(delta-1.d0).gt.wann_tol) lb1=.false.
      if(i.ne.j.and.abs(delta).gt.wann_tol) lb1=.false.
      if(.not.lb1) exit
    end do
    if(.not.lb1) exit
  end do 
  deallocate(amat)
  deallocate(umat)
  deallocate(vmat)
  deallocate(smat)
  deallocate(sing)
  if(lb1) exit
end do
if(.not.lb1) then
  write(*,*)
  write(*,'("Error(getw90nnkp): Could not satisfy the B1 condition with any shells.")')
  write(*,*)
  stop
end if
wann_nntot = 0
do id=1,nshellused
  wann_nntot = wann_nntot + nn(shell_list(id))
end do
if(allocated(wann_nnkp)) deallocate(wann_nnkp)
allocate(wann_nnkp(5,wann_nntot*nkpt))
inkp = 1
do ik=1,nkpt
  do ig=1,ngvec
    if (any(abs(ivg(:,ig)).gt.1)) cycle
    do jk=1,nkpt
      do jd=1,nshellused
        dist0 = sqrt( (vkc(1,ik)-(vkc(1,jk)+vgc(1,ig)))**2 + &
                      (vkc(2,ik)-(vkc(2,jk)+vgc(2,ig)))**2 + &
                      (vkc(3,ik)-(vkc(3,jk)+vgc(3,ig)))**2 )
        if (abs(dist0-dist(shell_list(jd))).le.wann_tol) then
          wann_nnkp(1,inkp) = ik
          wann_nnkp(2,inkp) = jk
          wann_nnkp(3:5,inkp) = ivg(:,ig)
          inkp = inkp + 1
        end if
      end do
    end do
  end do
end do
deallocate(nn)
deallocate(dist)
deallocate(shell_list)
deallocate(bvectors)
end subroutine getw90nnkp
!EOC
