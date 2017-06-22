
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine testcheck
implicit none
! local variables
logical exist
integer i,j,k,n
integer nv_,nv,vt_,vt,iv_,iv
real(8) rv_,rv,a,b
real(8) tol,t1,t2
complex(8) zv_,zv
character(256) fname_,fname,descr
n=0
do i=0,999
  write(fname_,'("TEST",I3.3,".OUT_")') i
  inquire(file=trim(fname_),exist=exist)
  if (exist) then
    write(fname,'("TEST",I3.3,".OUT")') i
    inquire(file=trim(fname),exist=exist)
    if (.not.exist) then
      write(*,*)
      write(*,'("Error(testcheck): file ",A," does not exist")') trim(fname)
      write(*,*)
      stop
    end if
    open(90,file=trim(fname_),action='READ',form='FORMATTED')
    open(91,file=trim(fname),action='READ',form='FORMATTED')
    read(90,*,err=10) descr
    read(91,*,err=20) descr
    read(90,*,err=10) vt_,nv_
    read(91,*,err=20) vt,nv
    if (vt_.ne.vt) then
      write(*,*)
      write(*,'("Error(testcheck): differing variable type")')
      write(*,'(" for quantity ''",A,"''")') trim(descr)
      write(*,'(" ",A," : ",I8)') trim(fname_),vt_
      write(*,'(" ",A,"  : ",I8)') trim(fname),vt
      write(*,*)
      stop
    end if
    if (nv_.ne.nv) then
      write(*,*)
      write(*,'("Error(testcheck): differing number of variables")')
      write(*,'(" for quantity ''",A,"''")') trim(descr)
      write(*,'(" ",A," : ",I8)') trim(fname_),nv_
      write(*,'(" ",A,"  : ",I8)') trim(fname),nv
      write(*,*)
      stop
    end if
    if (nv.le.0) then
      write(*,*)
      write(*,'("Error(testcheck): nv <= 0 : ",I8)') nv
      write(*,*)
      stop
    end if
    if (vt.eq.1) then
! integer variables
      do j=1,nv
        read(90,*,err=10) k,iv_
        if (j.ne.k) goto 10
        read(91,*,err=20) k,iv
        if (j.ne.k) goto 20
        if (iv.ne.iv_) then
          write(*,*)
          write(*,'("Error(testcheck): variable ",I8," is different")') j
          write(*,'(" for quantity ''",A,"''")') trim(descr)
          write(*,'(" ",A," : ",I8)') trim(fname_),iv_
          write(*,'(" ",A,"  : ",I8)') trim(fname),iv
          write(*,*)
          stop
        end if
      end do
    else if (vt.eq.2) then
! real variables
      read(90,*,err=10) tol
      read(91,*,err=20) tol
      do j=1,nv
        read(90,*,err=10) k,rv_
        if (j.ne.k) goto 10
        read(91,*,err=20) k,rv
        if (j.ne.k) goto 20
        t1=abs(rv_-rv)
        t2=abs(rv_)*tol
        if ((t1.gt.t2).and.(abs(rv_).gt.1.d-4)) then
          write(*,*)
          write(*,'("Error(testcheck): variable ",I8," outside tolerance")') j
          write(*,'(" for quantity ''",A,"''")') trim(descr)
          write(*,'(" ",A," (correct value)",T40," : ",G22.12)') trim(fname_), &
           rv_
          write(*,'(" ",A,T40," : ",G22.12)') trim(fname),rv
          write(*,'(" absolute difference",T40," : ",G22.12)') t1
          write(*,'(" required relative tolerance",T40," : ",G22.12)') tol
          write(*,'(" required absolute tolerance",T40," : ",G22.12)') t2
          write(*,*)
          stop
        end if
      end do
    else if (vt.eq.3) then
! complex variables
      read(90,*,err=10) tol
      read(91,*,err=20) tol
      do j=1,nv
        read(90,*,err=10) k,a,b
        zv_=cmplx(a,b,8)
        if (j.ne.k) goto 10
        read(91,*,err=20) k,a,b
        zv=cmplx(a,b,8)
        if (j.ne.k) goto 20
        t1=abs(zv_-zv)
        t2=abs(zv_)*tol
        if ((t1.gt.t2).and.(abs(rv_).gt.1.d-4)) then
          write(*,*)
          write(*,'("Error(testcheck): variable ",I8," outside tolerance")') j
          write(*,'(" for quantity ''",A,"''")') trim(descr)
          write(*,'(" ",A," (correct value)",T40," : ",2G22.12)') &
           trim(fname_),zv_
          write(*,'(" ",A,T40," : ",2G22.12)') trim(fname),zv
          write(*,'(" difference",T40," : ",G22.12)') t1
          write(*,'(" required relative tolerance",T40," : ",G22.12)') tol
          write(*,'(" required absolute tolerance",T40," : ",G22.12)') t2
          write(*,*)
          stop
        end if
      end do
    else
      write(*,*)
      write(*,'("Error(testcheck): variable type not defined : ",I8)') vt
      write(*,*)
      stop
    end if
    close(90)
    close(91)
    n=n+1
  end if
end do
if (n.eq.0) then
  write(*,*)
  write(*,'("Warning(testcheck): no tests found")')
else
  write(*,*)
  write(*,'("Info(testcheck): passed all tests")')
end if
return
10 continue
write(*,*)
write(*,'("Error(testcheck): error reading from ",A)') trim(fname_)
write(*,*)
stop
20 continue
write(*,*)
write(*,'("Error(testcheck): error reading from ",A)') trim(fname)
write(*,*)
stop
end subroutine

