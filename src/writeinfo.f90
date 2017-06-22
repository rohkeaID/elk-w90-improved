
! Copyright (C) 2002-2009 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writeinfo
! !INTERFACE:
subroutine writeinfo(fnum)
! !USES:
use modmain
use modmpi
use moddftu
use modrdm
! !INPUT/OUTPUT PARAMETERS:
!   fnum : unit specifier for INFO.OUT file (in,integer)
! !DESCRIPTION:
!   Outputs basic information about the run to the file {\tt INFO.OUT}. Does not
!   close the file afterwards.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!   Updated with DFT+U quantities July 2009 (FC)
!EOP
!BOC
implicit none
! arguments
integer fnum
! local variables
integer i,is,ia,k,l
character(10) dat,tim
write(fnum,'("+----------------------------+")')
write(fnum,'("| Elk version ",I1.1,".",I1.1,".",I2.2," started |")') version
write(fnum,'("+----------------------------+")')
call date_and_time(date=dat,time=tim)
write(fnum,*)
write(fnum,'("Date (YYYY-MM-DD) : ",A4,"-",A2,"-",A2)') dat(1:4),dat(5:6), &
 dat(7:8)
write(fnum,'("Time (hh:mm:ss)   : ",A2,":",A2,":",A2)') tim(1:2),tim(3:4), &
 tim(5:6)
if (np_mpi.gt.1) then
  write(fnum,*)
  write(fnum,'("Using MPI, number of processes : ",I8)') np_mpi
end if
if (notelns.gt.0) then
  write(fnum,*)
  write(fnum,'("Notes :")')
  do i=1,notelns
    write(fnum,'(A)') notes(i)
  end do
end if
write(fnum,*)
write(fnum,'("All units are atomic (Hartree, Bohr, etc.)")')
write(fnum,*)
select case(task)
case(0,1,28,29,200,201,350,351,440)
  if (trdstate) then
    write(fnum,'("+------------------------------------------+")')
    write(fnum,'("| Ground-state run resuming from STATE.OUT |")')
    write(fnum,'("+------------------------------------------+")')
  else
    write(fnum,'("+-------------------------------------------------+")')
    write(fnum,'("| Ground-state run starting from atomic densities |")')
    write(fnum,'("+-------------------------------------------------+")')
  end if
case(2,3)
  if (trdstate) then
    write(fnum,'("+---------------------------------------------------+")')
    write(fnum,'("| Geometry optimisation run resuming from STATE.OUT |")')
    write(fnum,'("+---------------------------------------------------+")')
  else
    write(fnum,'("+------------------------------------------------------+")')
    write(fnum,'("| Geometry optimisation starting from atomic densities |")')
    write(fnum,'("+------------------------------------------------------+")')
  end if
case(5)
  write(fnum,'("+-------------------------------+")')
  write(fnum,'("| Ground-state Hartree-Fock run |")')
  write(fnum,'("+-------------------------------+")')
case(300)
  write(fnum,'("+----------------------------------------------+")')
  write(fnum,'("| Reduced density matrix functional theory run |")')
  write(fnum,'("+----------------------------------------------+")')
case default
  write(*,*)
  write(*,'("Error(writeinfo): task not defined : ",I8)') task
  write(*,*)
  stop
end select
write(fnum,*)
write(fnum,'("Lattice vectors :")')
write(fnum,'(3G18.10)') avec(1,1),avec(2,1),avec(3,1)
write(fnum,'(3G18.10)') avec(1,2),avec(2,2),avec(3,2)
write(fnum,'(3G18.10)') avec(1,3),avec(2,3),avec(3,3)
write(fnum,*)
write(fnum,'("Reciprocal lattice vectors :")')
write(fnum,'(3G18.10)') bvec(1,1),bvec(2,1),bvec(3,1)
write(fnum,'(3G18.10)') bvec(1,2),bvec(2,2),bvec(3,2)
write(fnum,'(3G18.10)') bvec(1,3),bvec(2,3),bvec(3,3)
write(fnum,*)
write(fnum,'("Unit cell volume      : ",G18.10)') omega
write(fnum,'("Brillouin zone volume : ",G18.10)') omegabz
do is=1,nspecies
  write(fnum,*)
  write(fnum,'("Species : ",I4," (",A,")")') is,trim(spsymb(is))
  write(fnum,'(" parameters loaded from : ",A)') trim(spfname(is))
  write(fnum,'(" name : ",A)') trim(spname(is))
  write(fnum,'(" nuclear charge    : ",G18.10)') spzn(is)
  write(fnum,'(" electronic charge : ",G18.10)') spze(is)
  write(fnum,'(" atomic mass : ",G18.10)') spmass(is)
  write(fnum,'(" muffin-tin radius : ",G18.10)') rmt(is)
  write(fnum,'(" number of radial points in muffin-tin : ",I6)') nrmt(is)
  write(fnum,'(" number on inner part of muffin-tin    : ",I6)') nrmtinr(is)
  write(fnum,'(" atomic positions (lattice), magnetic fields (Cartesian) :")')
  do ia=1,natoms(is)
    write(fnum,'(I4," : ",3F12.8,"  ",3F12.8)') ia,atposl(:,ia,is), &
     bfcmt(:,ia,is)
  end do
end do
write(fnum,*)
write(fnum,'("Total number of atoms per unit cell : ",I4)') natmtot
write(fnum,*)
write(fnum,'("Spin treatment :")')
if (spinpol) then
  write(fnum,'(" spin-polarised")')
else
  write(fnum,'(" spin-unpolarised")')
end if
if (spinorb) then
  write(fnum,'(" spin-orbit coupling")')
end if
if (spincore) then
  write(fnum,'(" spin-polarised core states")')
end if
if (spinpol) then
  write(fnum,'(" global magnetic field (Cartesian) : ",3G18.10)') bfieldc
  if (ncmag) then
    write(fnum,'(" non-collinear magnetisation")')
  else
    write(fnum,'(" collinear magnetisation in z-direction")')
  end if
end if
if (spinsprl) then
  write(fnum,'(" spin-spiral state assumed")')
  write(fnum,'("  q-vector (lattice)   : ",3G18.10)') vqlss
  write(fnum,'("  q-vector (Cartesian) : ",3G18.10)') vqcss
  write(fnum,'("  q-vector length      : ",G18.10)') sqrt(vqcss(1)**2 &
   +vqcss(2)**2+vqcss(3)**2)
end if
if (fsmtype.ne.0) then
  write(fnum,'(" fixed spin moment (FSM) calculation, type : ",I4)') fsmtype
  if (fsmtype.lt.0) then
    write(fnum,'("  only moment direction is fixed")')
  end if
end if
if ((abs(fsmtype).eq.1).or.(abs(fsmtype).eq.3)) then
  write(fnum,'("  fixing total moment to (Cartesian) :")')
  write(fnum,'("  ",3G18.10)') momfix
end if
if ((abs(fsmtype).eq.2).or.(abs(fsmtype).eq.3)) then
  write(fnum,'("  fixing local muffin-tin moments to (Cartesian) :")')
  do is=1,nspecies
    write(fnum,'("  species : ",I4," (",A,")")') is,trim(spsymb(is))
    do ia=1,natoms(is)
      write(fnum,'("   ",I4,3G18.10)') ia,mommtfix(:,ia,is)
    end do
  end do
end if
if (ftmtype.ne.0) then
  write(fnum,*)
  write(fnum,'(" fixed tensor moment (FTM) calculation, type : ",I4)') ftmtype
end if
if (efieldpol) then
  write(fnum,*)
  write(fnum,'("Constant electric field applied across unit cell")')
  write(fnum,'(" field strength : ",3G18.10)') efieldc
end if
write(fnum,*)
write(fnum,'("Number of Bravais lattice symmetries : ",I4)') nsymlat
write(fnum,'("Number of crystal symmetries         : ",I4)') nsymcrys
if (tsyminv) then
  write(fnum,'("Crystal has inversion symmetry")')
else
  write(fnum,'("Crystal has no inversion symmetry")')
end if
if (tefvr) then
  write(fnum,'("Real symmetric eigensolver will be used")')
else
  write(fnum,'("Complex Hermitian eigensolver will be used")')
end if
write(fnum,*)
if (autokpt) then
  write(fnum,'("Radius of sphere used to determine k-point grid density : ",&
   &G18.10)') radkpt
end if
write(fnum,'("k-point grid : ",3I6)') ngridk
write(fnum,'("k-point offset : ",3G18.10)') vkloff
if (reducek.eq.0) then
  write(fnum,'("k-point set is not reduced")')
else if (reducek.eq.1) then
  write(fnum,'("k-point set is reduced with full crystal symmetry group")')
else if (reducek.eq.2) then
  write(fnum,'("k-point set is reduced with symmorphic symmetries only")')
else
  write(*,*)
  write(*,'("Error(writeinfo): undefined k-point reduction type : ",I8)') &
   reducek
  write(*,*)
  stop
end if
write(fnum,'("Total number of k-points : ",I8)') nkpt
write(fnum,*)
write(fnum,'("Muffin-tin radius times maximum |G+k| : ",G18.10)') rgkmax
if (isgkmax.eq.-2) then
  write(fnum,'(" using gkmax = rgkmax / 2")')
else
  if ((isgkmax.ge.1).and.(isgkmax.le.nspecies)) then
    write(fnum,'(" using radius of species ",I4," (",A,")")') isgkmax, &
     trim(spsymb(isgkmax))
  else if (isgkmax.eq.-1) then
    write(fnum,'(" using average radius")')
  else
    write(fnum,'(" using smallest radius")')
  end if
end if
write(fnum,'("Maximum |G+k| for APW functions       : ",G18.10)') gkmax
write(fnum,'("Maximum (1/2)|G+k|^2                  : ",G18.10)') 0.5d0*gkmax**2
write(fnum,'("Maximum |G| for potential and density : ",G18.10)') gmaxvr
write(fnum,'("Constant for pseudocharge density : ",I4)') npsd
write(fnum,'("Radial integration step length : ",I4)') lradstp
write(fnum,*)
write(fnum,'("G-vector grid sizes : ",3I6)') ngridg(:)
write(fnum,'("Number of G-vectors : ",I8)') ngvec
write(fnum,*)
write(fnum,'("Maximum angular momentum used for")')
write(fnum,'(" APW functions                      : ",I4)') lmaxapw
write(fnum,'(" H and O matrix elements outer loop : ",I4)') lmaxmat
write(fnum,'(" potential and density              : ",I4)') lmaxvr
write(fnum,'(" inner part of muffin-tin           : ",I4)') lmaxinr
write(fnum,*)
write(fnum,'("Total nuclear charge    : ",G18.10)') chgzn
write(fnum,'("Total core charge       : ",G18.10)') chgcrtot
write(fnum,'("Total valence charge    : ",G18.10)') chgval
write(fnum,'("Total excess charge     : ",G18.10)') chgexs
write(fnum,'("Total electronic charge : ",G18.10)') chgtot
write(fnum,*)
write(fnum,'("Effective Wigner radius, r_s : ",G18.10)') rwigner
write(fnum,*)
write(fnum,'("Number of empty states         : ",I4)') nempty
write(fnum,'("Total number of valence states : ",I4)') nstsv
write(fnum,'("Total number of core states    : ",I4)') nstcr
write(fnum,*)
if (lorbcnd) then
  write(fnum,'("Conduction state local-orbitals added automatically")')
end if
write(fnum,'("Total number of local-orbitals : ",I4)') nlotot
if (tefvit) then
  write(fnum,*)
  write(fnum,'("Using iterative diagonalisation for the first-variational &
   &eigenvalue equation")')
end if
write(fnum,*)
if (task.eq.5) then
  write(fnum,'("Hartree-Fock calculation using Kohn-Sham states")')
  if (hybrid) then
    write(fnum,'(" hybrid functional, coefficient : ",G18.10)') hybridc
  end if
end if
if (xctype(1).lt.0) then
  write(fnum,'("Optimised effective potential (OEP) and exact exchange (EXX)")')
  write(fnum,'(" Phys. Rev. B 53, 7024 (1996)")')
  write(fnum,'("Correlation functional : ",3I6)') abs(xctype(1)),xctype(2:3)
  write(fnum,'(" ",A)') trim(xcdescr)
else
  write(fnum,'("Exchange-correlation functional : ",3I6)') xctype(:)
  write(fnum,'(" ",A)') trim(xcdescr)
end if
if (xcgrad.eq.0) then
  write(fnum,'(" Local density approximation (LDA)")')
else if ((xcgrad.eq.1).or.(xcgrad.eq.2)) then
  write(fnum,'(" Generalised gradient approximation (GGA)")')
else if (xcgrad.eq.3) then
  write(fnum,'(" meta-GGA; using kinetic energy density")')
end if
if (dftu.ne.0) then
  write(fnum,*)
  write(fnum,'("DFT+U calculation")')
  if (dftu.eq.1) then
    write(fnum,'(" fully localised limit (FLL)")')
    write(fnum,'(" see Phys. Rev. B 52, R5467 (1995)")')
  else if (dftu.eq.2) then
    write(fnum,'(" around mean field (AMF)")')
    write(fnum,'(" see Phys. Rev. B 49, 14211 (1994)")')
  else if (dftu.eq.3) then
    write(fnum,'(" interpolation between FLL and AMF")')
    write(fnum,'(" see Phys. Rev. B 67, 153106 (2003)")')
  else
    write(*,*)
    write(*,'("Error(writeinfo): dftu not defined : ",I8)') dftu
    write(*,*)
    stop
  end if
  do i=1,ndftu
    is=idftu(1,i)
    l=idftu(2,i)
    if (inpdftu.eq.1) then
      write(fnum,'(" species : ",I4," (",A,")",", l = ",I2,", U = ",F12.8, &
       &", J = ",F12.8)') is,trim(spsymb(is)),l,ujdu(1,i),ujdu(2,i)
    else if (inpdftu.eq.2) then
      write(fnum,'(" species : ",I4," (",A,")",", l = ",I2)') is, &
       trim(spsymb(is)),l
      write(fnum,'(" Slater integrals are provided as input")')
      do k=0,2*l,2
        write(fnum,'(" F^(",I1,") = ",F12.8)') k,fdu(k,i)
      end do
    else if (inpdftu.eq.3) then
      write(fnum,'(" species : ",I4," (",A,")",", l = ",I2)') is, &
       trim(spsymb(is)),l
      write(fnum,'(" Racah parameters are provided as input")')
      do k=0,l
        write(fnum,'(" E^(",I1,") = ",F12.8)') k,edu(k,i)
      end do
    else if (inpdftu.eq.4) then
      write(fnum,'(" species : ",I4," (",A,")",", l = ",I2)') is, &
       trim(spsymb(is)),l
      write(fnum,'(" Slater integrals are calculated by means of Yukawa &
       &potential")')
      write(fnum,'(" Yukawa potential screening length (a.u^-1) : ",F12.8)') &
       lambdadu(i)
    else if(inpdftu.eq.5) then
      write(fnum,'(" species : ",I4," (",A,")",", l = ",I2)') is, &
       trim(spsymb(is)),l
      write(fnum,'(" Slater integrals are calculated by means of Yukawa &
       &potential")')
      write(fnum,'(" Yukawa potential screening length corresponds to U = ",&
       &F12.8)') udufix(i)
    end if
  end do
end if
if (task.eq.300) then
  write(fnum,*)
  write(fnum,'("RDMFT calculation")')
  write(fnum,'(" see arXiv:0801.3787v1 [cond-mat.mtrl-sci]")')
  write(fnum,'(" RDMFT exchange-correlation type : ",I4)') rdmxctype
  if (rdmxctype.eq.1) then
    write(fnum,'("  Hartree-Fock functional")')
  else if (rdmxctype.eq.2) then
    write(fnum,'("  Power functional, exponent : ",G18.10)') rdmalpha
  end if
end if
write(fnum,*)
write(fnum,'("Smearing type : ",I4)') stype
write(fnum,'(" ",A)') trim(sdescr)
if (autoswidth) then
  write(fnum,'("Automatic determination of smearing width")')
else
  write(fnum,'("Smearing width : ",G18.10)') swidth
  write(fnum,'("Effective electronic temperature (K) : ",G18.10)') swidth/kboltz
end if
write(fnum,*)
write(fnum,'("Mixing type : ",I4)') mixtype
write(fnum,'(" ",A)') trim(mixdescr)
call flushifc(fnum)
return
end subroutine
!EOC

