
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writespecies(symb,name,zn,mass,rmin,rm,rmax,nrm,nst,n,l,k,occ,eval)
use modmain
use modmpi
implicit none
! arguments
character(*), intent(in) :: symb,name
real(8), intent(in) :: zn,mass
real(8), intent(in) :: rmin,rm,rmax
integer, intent(in) :: nrm,nst
integer, intent(in) :: n(nst),l(nst),k(nst)
real(8), intent(in) :: occ(nst)
real(8), intent(in) :: eval(nst)
! local variables
integer lmax,nlo
integer ist,jst,i
logical core(maxstsp),lorb(maxstsp)
! default APW band energy
real(8), parameter :: e0=0.15d0
! find which states belong to core
do ist=1,nst
  if (eval(ist).lt.ecvcut) then
    core(ist)=.true.
  else
    core(ist)=.false.
  end if
end do
! check that the state for same n and l but different k is also core
do ist=1,nst
  if (core(ist)) then
    do jst=1,nst
      if ((n(ist).eq.n(jst)).and.(l(ist).eq.l(jst))) core(jst)=.true.
    end do
  end if
end do
lmax=1
do ist=1,nst
  if (.not.core(ist)) lmax=max(lmax,l(ist))
end do
! determine the local orbitals
nlo=lmax+1
lorb(:)=.false.
do ist=1,nst
  if (.not.core(ist)) then
    if ((l(ist).eq.0).or.(l(ist).lt.k(ist))) then
      if ((eval(ist).lt.esccut).or.(l(ist).ge.2)) then
        lorb(ist)=.true.
        nlo=nlo+1
      end if
    end if
  end if
end do
if (mp_mpi) then
  open(55,file=trim(symb)//'.in',form='FORMATTED')
  write(55,'(" ''",A,"''",T45,": spsymb")') trim(symb)
  write(55,'(" ''",A,"''",T45,": spname")') trim(name)
  write(55,'(G14.6,T45,": spzn")') zn
  write(55,'(G18.10,T45,": spmass")') mass
  write(55,'(G14.6,2F10.4,I6,T45,": rminsp, rmt, rmaxsp, nrmt")') rmin,rm, &
   rmax,nrm
  write(55,'(I4,T45,": nstsp")') nst
  write(55,'(3I4,G14.6,L1,T45,": nsp, lsp, ksp, occsp, spcore")') n(1),l(1), &
   k(1),occ(1),core(1)
  do ist=2,nst
    write(55,'(3I4,G14.6,L1)') n(ist),l(ist),k(ist),occ(ist),core(ist)
  end do
  write(55,'(I4,T45,": apword")') 1
  write(55,'(F10.4,I4,"  ",L1,T45,": apwe0, apwdm, apwve")') e0,0,.false.
  write(55,'(I4,T45,": nlx")') 0
  write(55,'(I4,T45,": nlorb")') nlo
  do i=0,lmax
    write(55,'(2I4,T45,": lorbl, lorbord")') i,2
    write(55,'(F10.4,I4,"  ",L1,T45,": lorbe0, lorbdm, lorbve")') e0,0,.false.
    write(55,'(F10.4,I4,"  ",L1)') e0,1,.false.
  end do
  do ist=1,nst
    if (lorb(ist)) then
      write(55,'(2I4,T45,": lorbl, lorbord")') l(ist),3
      write(55,'(F10.4,I4,"  ",L1,T45,": lorbe0, lorbdm, lorbve")') e0,0,.false.
      write(55,'(F10.4,I4,"  ",L1)') e0,1,.false.
      write(55,'(F10.4,I4,"  ",L1)') eval(ist),0,.true.
    end if
  end do
  close(55)
end if
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
return
end subroutine

