
! Copyright (C) 2002-2009 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modmain

!----------------------------!
!     lattice parameters     !
!----------------------------!
! lattice vectors stored column-wise
real(8) avec(3,3)
! magnitude of random displacements added to lattice vectors
real(8) rndavec
! inverse of lattice vector matrix
real(8) ainv(3,3)
! reciprocal lattice vectors
real(8) bvec(3,3)
! inverse of reciprocal lattice vector matrix
real(8) binv(3,3)
! unit cell volume
real(8) omega
! Brillouin zone volume
real(8) omegabz
! any vector with length less than epslat is considered zero
real(8) epslat

!--------------------------!
!     atomic variables     !
!--------------------------!
! maximum allowed species
integer, parameter :: maxspecies=8
! maximum allowed atoms per species
integer, parameter :: maxatoms=200
! number of species
integer nspecies
! number of atoms for each species
integer natoms(maxspecies)
! maximum number of atoms over all the species
integer natmmax
! total number of atoms
integer natmtot
! index to atoms and species
integer idxas(maxatoms,maxspecies)
! inverse atoms and species indices
integer idxis(maxatoms*maxspecies)
integer idxia(maxatoms*maxspecies)
! molecule is .true. is the system is an isolated molecule
logical molecule
! primcell is .true. if primitive unit cell is to be found automatically
logical primcell
! atomic positions in lattice coordinates
real(8) atposl(3,maxatoms,maxspecies)
! atomic positions in Cartesian coordinates
real(8) atposc(3,maxatoms,maxspecies)
! magnitude of random displacements added to the atomic positions
real(8) rndatposc

!----------------------------------!
!     atomic species variables     !
!----------------------------------!
! species files path
character(256) sppath
! species filenames
character(256) spfname(maxspecies)
! species name
character(256) spname(maxspecies)
! species symbol
character(256) spsymb(maxspecies)
! species nuclear charge
real(8) spzn(maxspecies)
! ptnucl is .true. if the nuclei are to be treated as point charges, if .false.
! the nuclei have a finite spherical distribution
logical ptnucl
! nuclear radius
real(8) rnucl(maxspecies)
! nuclear volume
real(8) vnucl(maxspecies)
! number of radial mesh points to nuclear radius
integer nrnucl(maxspecies)
! number of coarse radial mesh points to nuclear radius
integer nrcnucl(maxspecies)
! species electronic charge
real(8) spze(maxspecies)
! species mass
real(8) spmass(maxspecies)
! smallest radial point for each species
real(8) rminsp(maxspecies)
! effective infinity for species
real(8) rmaxsp(maxspecies)
! number of radial points to effective infinity for each species
integer nrsp(maxspecies)
! maximum nrsp over all the species
integer nrspmax
! maximum allowed states for each species
integer, parameter :: maxstsp=40
! number of states for each species
integer nstsp(maxspecies)
! maximum nstsp over all the species
integer nstspmax
! core-valence cut-off energy for species file generation
real(8) ecvcut
! semi-core-valence cut-off energy for species file generation
real(8) esccut
! state principle quantum number for each species
integer nsp(maxstsp,maxspecies)
! state l value for each species
integer lsp(maxstsp,maxspecies)
! state k value for each species
integer ksp(maxstsp,maxspecies)
! spcore is .true. if species state is core
logical spcore(maxstsp,maxspecies)
! total number of core states
integer nstcr
! state eigenvalue for each species
real(8) evalsp(maxstsp,maxspecies)
! state occupancy for each species
real(8) occsp(maxstsp,maxspecies)
! species radial mesh
real(8), allocatable :: rsp(:,:)
! r^2 on radial mesh
real(8), allocatable :: r2sp(:,:)
! species charge density
real(8), allocatable :: rhosp(:,:)
! species self-consistent potential
real(8), allocatable :: vrsp(:,:)
! exchange-correlation type for atomic species (the converged ground-state of
! the crystal does not depend on this choice)
integer xctsp(3)

!---------------------------------------------------------------!
!     muffin-tin radial mesh and angular momentum variables     !
!---------------------------------------------------------------!
! scale factor for number of muffin-tin points
real(8) nrmtscf
! number of muffin-tin radial points for each species
integer nrmt(maxspecies)
! maximum nrmt over all the species
integer nrmtmax
! optional default muffin-tin radius for all atoms
real(8) rmtall
! minimum allowed distance between muffin-tin surfaces
real(8) rmtdelta
! muffin-tin radii
real(8) rmt(maxspecies)
! total muffin-tin volume
real(8) omegamt
! radial step length for coarse mesh
integer lradstp
! number of coarse radial mesh points
integer nrcmt(maxspecies)
! maximum nrcmt over all the species
integer nrcmtmax
! coarse muffin-tin radial mesh
real(8), allocatable :: rcmt(:,:)
! r^2 on coarse radial mesh
real(8), allocatable :: r2cmt(:,:)
! maximum allowable angular momentum for augmented plane waves
integer, parameter :: maxlapw=50
! maximum angular momentum for augmented plane waves
integer lmaxapw
! (lmaxapw+1)^2
integer lmmaxapw
! maximum angular momentum on the outer part of the muffin-tin
integer lmaxo
! (lmaxo+1)^2
integer lmmaxo
! maximum angular momentum on the inner part of the muffin-tin
integer lmaxi
! (lmaxi+1)^2
integer lmmaxi
! fraction of muffin-tin radius which constitutes the inner part
real(8) fracinr
! number of fine/coarse radial points on the inner part of the muffin-tin
integer nrmti(maxspecies),nrcmti(maxspecies)
! index to (l,m) pairs
integer, allocatable :: idxlm(:,:)
! inverse index to (l,m) pairs
integer, allocatable :: idxil(:),idxim(:)
! number of fine/coarse points in packed muffin-tins
integer npmti(maxspecies),npmt(maxspecies)
integer npcmti(maxspecies),npcmt(maxspecies)
! maximum number of points over all packed muffin-tins
integer npmtmax,npcmtmax

!--------------------------------!
!     spin related variables     !
!--------------------------------!
! spinpol is .true. for spin-polarised calculations
logical spinpol
! spinorb is .true. for spin-orbit coupling
logical spinorb
! scale factor of spin-orbit coupling term in Hamiltonian
real(8) socscf
! dimension of magnetisation and magnetic vector fields (1 or 3)
integer ndmag
! ncmag is .true. if the magnetisation is non-collinear, i.e. when ndmag = 3
logical ncmag
! if cmagz is .true. then collinear magnetism along the z-axis is enforced
logical cmagz
! spcpl is .true. if the up and down spins are coupled
logical spcpl
! fixed spin moment type
!  0      : none
!  1 (-1) : total moment (direction)
!  2 (-2) : individual muffin-tin moments (direction)
!  3 (-3) : total and muffin-tin moments (direction)
integer fsmtype
! fixed total spin magnetic moment
real(8) momfix(3)
! fixed spin moment global effective field in Cartesian coordinates
real(8) bfsmc(3)
! muffin-tin fixed spin moments
real(8) mommtfix(3,maxatoms,maxspecies)
! muffin-tin fixed spin moment effective fields in Cartesian coordinates
real(8), allocatable :: bfsmcmt(:,:)
! fixed spin moment field step size
real(8) taufsm
! second-variational spinor dimension (1 or 2)
integer nspinor
! global external magnetic field in Cartesian coordinates
real(8) bfieldc(3)
! initial field
real(8) bfieldc0(3)
! external magnetic field in each muffin-tin in Cartesian coordinates
real(8) bfcmt(3,maxatoms,maxspecies)
! initial field as read in from input file
real(8) bfcmt0(3,maxatoms,maxspecies)
! magnitude of random vectors added to muffin-tin fields
real(8) rndbfcmt
! external magnetic fields are multiplied by reducebf after each s.c. loop
real(8) reducebf
! spinsprl is .true. if a spin-spiral is to be calculated
logical spinsprl
! ssdph is .true. if the muffin-tin spin-spiral magnetisation is de-phased
logical ssdph
! number of spin-dependent first-variational functions per state
integer nspnfv
! map from second- to first-variational spin index
integer jspnfv(2)
! spin-spiral q-vector in lattice coordinates
real(8) vqlss(3)
! spin-spiral q-vector in Cartesian coordinates
real(8) vqcss(3)
! current q-point in spin-spiral supercell calculation
integer iqss
! number of primitive unit cells in spin-spiral supercell
integer nscss
! number of fixed spin direction points on the sphere for finding the magnetic
! anisotropy energy (MAE)
integer npmae0,npmae
! (theta,phi) coordinates for each MAE direction
real(8), allocatable :: tpmae(:,:)

!---------------------------------------------!
!     electric field and vector potential     !
!---------------------------------------------!
! tefield is .true. if a polarising constant electric field is applied
logical tefield
! electric field vector in Cartesian coordinates
real(8) efieldc(3)
! electric field vector in lattice coordinates
real(8) efieldl(3)
! tafield is .true. if a constant vector potential is applied
logical tafield
! vector potential A-field which couples to paramagnetic current
real(8) afieldc(3)
! A-field in lattice coordinates
real(8) afieldl(3)

!----------------------------!
!     symmetry variables     !
!----------------------------!
! type of symmetry allowed for the crystal
!  0 : only the identity element is used
!  1 : full symmetry group is used
!  2 : only symmorphic symmetries are allowed
integer symtype
! number of Bravais lattice point group symmetries
integer nsymlat
! Bravais lattice point group symmetries
integer symlat(3,3,48)
! determinants of lattice symmetry matrices (1 or -1)
integer symlatd(48)
! index to inverses of the lattice symmetries
integer isymlat(48)
! lattice point group symmetries in Cartesian coordinates
real(8) symlatc(3,3,48)
! tshift is .true. if atomic basis is allowed to be shifted
logical tshift
! tsyminv is .true. if the crystal has inversion symmetry
logical tsyminv
! maximum of symmetries allowed
integer, parameter :: maxsymcrys=1024
! number of crystal symmetries
integer nsymcrys
! crystal symmetry translation vector in lattice coordinates
real(8) vtlsymc(3,maxsymcrys)
! tv0symc is .true. if the translation vector is zero
logical tv0symc(maxsymcrys)
! spatial rotation element in lattice point group for each crystal symmetry
integer lsplsymc(maxsymcrys)
! global spin rotation element in lattice point group for each crystal symmetry
integer lspnsymc(maxsymcrys)
! equivalent atom index for each crystal symmetry
integer, allocatable :: ieqatom(:,:,:)
! eqatoms(ia,ja,is) is .true. if atoms ia and ja are equivalent
logical, allocatable :: eqatoms(:,:,:)
! number of site symmetries
integer, allocatable :: nsymsite(:)
! site symmetry spatial rotation element in lattice point group
integer, allocatable :: lsplsyms(:,:)
! site symmetry global spin rotation element in lattice point group
integer, allocatable :: lspnsyms(:,:)

!----------------------------!
!     G-vector variables     !
!----------------------------!
! G-vector cut-off for interstitial potential and density
real(8) gmaxvr
! G-vector grid sizes
integer ngridg(3)
! G-vector grid sizes for coarse grid (G < 2*gkmax)
integer ngdc(3)
! total number of G-vectors
integer ngtot
! total number of G-vectors for coarse grid (G < 2*gkmax)
integer ngtc
! integer grid intervals for each direction
integer intgv(2,3)
! number of G-vectors with G < gmaxvr
integer ngvec
! number of G-vectors for coarse grid (G < 2*gkmax)
integer ngvc
! G-vector integer coordinates (i1,i2,i3)
integer, allocatable :: ivg(:,:)
! map from (i1,i2,i3) to G-vector index
integer, allocatable :: ivgig(:,:,:)
! map from G-vector index to FFT array
integer, allocatable :: igfft(:)
! map from G-vector index to FFT array for coarse grid (G < 2*gkmax)
integer, allocatable :: igfc(:)
! G-vectors in Cartesian coordinates
real(8), allocatable :: vgc(:,:)
! length of G-vectors
real(8), allocatable :: gc(:)
! spherical Bessel functions j_l(|G|R_mt)
real(8), allocatable :: jlgrmt(:,:,:)
! spherical harmonics of the G-vectors
complex(8), allocatable :: ylmg(:,:)
! structure factors for the G-vectors
complex(8), allocatable :: sfacg(:,:)
! smooth step function form factors for all species and G-vectors
real(8), allocatable :: ffacg(:,:)
! characteristic function in G-space: 0 inside the muffin-tins and 1 outside
complex(8), allocatable :: cfunig(:)
! characteristic function in real-space: 0 inside the muffin-tins and 1 outside
real(8), allocatable :: cfunir(:)
! characteristic function in real-space for coarse grid (G < 2*gkmax)
real(8), allocatable :: cfrc(:)
! Coulomb Green's function in G-space = 4 pi / G^2
real(8), allocatable :: gclg(:)

!---------------------------!
!     k-point variables     !
!---------------------------!
! autokpt is .true. if the k-point set is determined automatically
logical autokpt
! radius of sphere used to determine k-point density when autokpt is .true.
real(8) radkpt
! k-point grid sizes
integer ngridk(3)
! k-point offset
real(8) vkloff(3)
! corners of box in lattice coordinates containing the k-points, the zeroth
! vector is the origin
real(8) kptboxl(3,0:3)
! type of reduction to perform on k-point set
!  0 : no reduction
!  1 : reduce with full crystal symmetry group
!  2 : reduce with symmorphic symmetries only
integer reducek
! number of point group symmetries used for k-point reduction
integer nsymkpt
! point group symmetry matrices used for k-point reduction
integer symkpt(3,3,48)
! total number of reduced k-points
integer nkpt
! total number of non-reduced k-points
integer nkptnr
! locations of k-points on integer grid
integer, allocatable :: ivk(:,:)
! map from (i1,i2,i3) to reduced k-point index
integer, allocatable :: ivkik(:,:,:)
! map from (i1,i2,i3) to non-reduced k-point index
integer, allocatable :: ivkiknr(:,:,:)
! k-points in lattice coordinates
real(8), allocatable :: vkl(:,:)
! k-points in Cartesian coordinates
real(8), allocatable :: vkc(:,:)
! reduced k-point weights
real(8), allocatable :: wkpt(:)
! weight of each non-reduced k-point
real(8) wkptnr
! k-point at which to determine effective mass tensor
real(8) vklem(3)
! displacement size for computing the effective mass tensor
real(8) deltaem
! number of displacements in each direction
integer ndspem

!------------------------------!
!     G+k-vector variables     !
!------------------------------!
! species for which the muffin-tin radius will be used for calculating gkmax
integer isgkmax
! smallest muffin-tin radius times gkmax
real(8) rgkmax
! maximum |G+k| cut-off for APW functions
real(8) gkmax
! number of G+k-vectors for augmented plane waves
integer, allocatable :: ngk(:,:)
! maximum number of G+k-vectors over all k-points
integer ngkmax
! index from G+k-vectors to G-vectors
integer, allocatable :: igkig(:,:,:)
! G+k-vectors in lattice coordinates
real(8), allocatable :: vgkl(:,:,:,:)
! G+k-vectors in Cartesian coordinates
real(8), allocatable :: vgkc(:,:,:,:)
! length of G+k-vectors
real(8), allocatable :: gkc(:,:,:)
! (theta, phi) coordinates of G+k-vectors
real(8), allocatable :: tpgkc(:,:,:,:)
! structure factors for the G+k-vectors
complex(8), allocatable :: sfacgk(:,:,:,:)

!---------------------------!
!     q-point variables     !
!---------------------------!
! q-point grid sizes
integer ngridq(3)
! integer grid intervals for the q-points
integer intq(2,3)
! type of reduction to perform on q-point set (see reducek)
integer reduceq
! number of point group symmetries used for q-point reduction
integer nsymqpt
! point group symmetry matrices used for q-point reduction
integer symqpt(3,3,48)
! total number of reduced q-points
integer nqpt
! total number of non-reduced q-points
integer nqptnr
! map from non-reduced grid to reduced index
integer, allocatable :: iqmap(:,:,:)
! map from non-reduced grid to non-reduced index
integer, allocatable :: iqmapnr(:,:,:)
! locations of q-points on integer grid
integer, allocatable :: ivq(:,:)
! map from (i1,i2,i3) to q-vector index
integer, allocatable :: ivqiq(:,:,:)
! map from q-vector index to FFT array
integer, allocatable :: iqfft(:)
! q-points in lattice coordinates
real(8), allocatable :: vql(:,:)
! q-points in Cartesian coordinates
real(8), allocatable :: vqc(:,:)
! length of q-vectors
real(8), allocatable :: qc(:)
! q-point weights
real(8), allocatable :: wqpt(:)
! weight for each non-reduced q-point
real(8) wqptnr
! index of q = 0 point
integer iq0
! regularised Coulomb Green's function in q-space
real(8), allocatable :: gclq(:)

!-----------------------------------------------------!
!     spherical harmonic transform (SHT) matrices     !
!-----------------------------------------------------!
! trotsht is .true. if the spherical cover used for the SHT is to be rotated
logical trotsht
data trotsht / .false. /
! spherical cover rotation matrix
real(8) rotsht(3,3)
! real backward SHT matrix for lmaxi
real(8), allocatable :: rbshti(:,:)
! real forward SHT matrix for lmaxi
real(8), allocatable :: rfshti(:,:)
! real backward SHT matrix for lmaxo
real(8), allocatable :: rbshto(:,:)
! real forward SHT matrix for lmaxo
real(8), allocatable :: rfshto(:,:)
! complex backward SHT matrix for lmaxi
complex(8), allocatable :: zbshti(:,:)
! complex forward SHT matrix for lmaxi
complex(8), allocatable :: zfshti(:,:)
! complex backward SHT matrix for lmaxo
complex(8), allocatable :: zbshto(:,:)
! complex forward SHT matrix for lmaxo
complex(8), allocatable :: zfshto(:,:)

!---------------------------------------------------------------!
!     density, potential and exchange-correlation variables     !
!---------------------------------------------------------------!
! exchange-correlation functional type
integer xctype(3)
! exchange-correlation functional description
character(512) xcdescr
! exchange-correlation functional spin requirement
integer xcspin
! exchange-correlation functional density gradient requirement
integer xcgrad
! small constant used to stabilise non-collinear GGA
real(8) dncgga
! muffin-tin and interstitial charge density
real(8), allocatable :: rhomt(:,:),rhoir(:)
! trhonorm is .true. if the density is to be normalised after every iteration
logical trhonorm
! muffin-tin and interstitial magnetisation vector field
real(8), allocatable :: magmt(:,:,:),magir(:,:)
! tcden is .true. if the current density is to be calculated
logical tcden
! muffin-tin and interstitial current density vector field
real(8), allocatable :: cdmt(:,:,:),cdir(:,:)
! amount of smoothing to be applied to the exchange-correlation potentials and
! magnetic field
integer msmooth
! muffin-tin and interstitial Coulomb potential
real(8), allocatable :: vclmt(:,:),vclir(:)
! Poisson solver pseudocharge density constant
integer npsd
! lmaxo+npsd+1
integer lnpsd
! muffin-tin and interstitial exchange energy density
real(8), allocatable :: exmt(:,:),exir(:)
! muffin-tin and interstitial correlation energy density
real(8), allocatable :: ecmt(:,:),ecir(:)
! muffin-tin and interstitial exchange-correlation potential
real(8), allocatable :: vxcmt(:,:),vxcir(:)
! muffin-tin and interstitial Kohn-Sham effective potential
real(8), allocatable :: vsmt(:,:),vsir(:)
! G-space interstitial Kohn-Sham effective potential
complex(8), allocatable :: vsig(:)
! muffin-tin and interstitial exchange-correlation magnetic field
real(8), allocatable :: bxcmt(:,:,:),bxcir(:,:)
! muffin-tin Kohn-Sham effective magnetic field in spherical coordinates
real(8), allocatable :: bsmt(:,:,:)
! interstitial Kohn-Sham effective magnetic field
real(8), allocatable :: bsir(:,:)
! nosource is .true. if the field is to be made source-free
logical nosource
! tssxc is .true. if scaled spin exchange-correlation (SSXC) is to be used
logical tssxc
! SSXC scaling factor
real(8) ssxc
! spin-orbit coupling radial function
real(8), allocatable :: socfr(:,:)
! kinetic energy density
real(8), allocatable :: taumt(:,:,:),tauir(:,:)
! core kinetic energy density
real(8), allocatable :: taucr(:,:,:)
! taudft is .true. if meta-GGA is to be treated as a tau-DFT functional
logical taudft
! tau-DFT exchange-correlation potential
real(8), allocatable :: wxcmt(:,:),wxcir(:)
! tau-DFT Kohn-Sham potential
real(8), allocatable :: wsmt(:,:),wsir(:)
! Tran-Blaha '09 constant c [Phys. Rev. Lett. 102, 226401 (2009)]
real(8) c_tb09
! tc_tb09 is .true. if the Tran-Blaha constant has been read in
logical tc_tb09
! if trdstate is .true. the density and potential can be read from STATE.OUT
logical trdstate
! temperature in degrees kelvin
real(8) tempk

!--------------------------!
!     mixing variables     !
!--------------------------!
! type of mixing to use for the potential
integer mixtype
! mixing type description
character(256) mixdescr
! adaptive mixing parameter
real(8) beta0
real(8) betamax
! subspace dimension for Broyden mixing
integer mixsdb
! Broyden mixing parameters alpha and w0
real(8) broydpm(2)

!----------------------------------------------!
!     charge, moment and current variables     !
!----------------------------------------------!
! tolerance for error in total charge
real(8) epschg
! total nuclear charge
real(8) chgzn
! core charges
real(8) chgcr(maxspecies)
! total core charge
real(8) chgcrtot
! core leakage charge
real(8), allocatable :: chgcrlk(:)
! total valence charge
real(8) chgval
! excess charge
real(8) chgexs
! total charge
real(8) chgtot
! calculated total charge
real(8) chgcalc
! interstitial region charge
real(8) chgir
! muffin-tin charges
real(8), allocatable :: chgmt(:)
! total muffin-tin charge
real(8) chgmttot
! effective Wigner radius
real(8) rwigner
! total moment
real(8) momtot(3)
! total moment magnitude
real(8) momtotm
! interstitial region moment
real(8) momir(3)
! muffin-tin moments
real(8), allocatable :: mommt(:,:)
! total muffin-tin moment
real(8) mommttot(3)
! total current
real(8) curtot(3)
! total current magnitude
real(8) curtotm

!-----------------------------------------!
!     APW and local-orbital variables     !
!-----------------------------------------!
! energy step used for numerical calculation of energy derivatives
real(8) deapwlo
! maximum allowable APW order
integer, parameter :: maxapword=4
! APW order
integer apword(0:maxlapw,maxspecies)
! maximum of apword over all angular momenta and species
integer apwordmax
! APW initial linearisation energies
real(8) apwe0(maxapword,0:maxlapw,maxspecies)
! APW linearisation energies
real(8), allocatable :: apwe(:,:,:)
! APW derivative order
integer apwdm(maxapword,0:maxlapw,maxspecies)
! apwve is .true. if the linearisation energies are allowed to vary
logical apwve(maxapword,0:maxlapw,maxspecies)
! APW radial functions
real(8), allocatable :: apwfr(:,:,:,:,:)
! derivate of radial functions at the muffin-tin surface
real(8), allocatable :: apwdfr(:,:,:)
! maximum number of local-orbitals
integer, parameter :: maxlorb=200
! maximum allowable local-orbital order
integer, parameter :: maxlorbord=5
! number of local-orbitals
integer nlorb(maxspecies)
! maximum nlorb over all species
integer nlomax
! total number of local-orbitals
integer nlotot
! local-orbital order
integer lorbord(maxlorb,maxspecies)
! local-orbital angular momentum
integer lorbl(maxlorb,maxspecies)
! maximum lorbl over all species
integer lolmax
! (lolmax+1)^2
integer lolmmax
! local-orbital initial energies
real(8) lorbe0(maxlorbord,maxlorb,maxspecies)
! local-orbital energies
real(8), allocatable :: lorbe(:,:,:)
! local-orbital derivative order
integer lorbdm(maxlorbord,maxlorb,maxspecies)
! lorbve is .true. if the linearisation energies are allowed to vary
logical lorbve(maxlorbord,maxlorb,maxspecies)
! local-orbital radial functions
real(8), allocatable :: lofr(:,:,:,:)
! band energy search tolerance
real(8) epsband
! maximum allowed change in energy during band energy search; enforced only if
! default energy is less than zero
real(8) demaxbnd
! minimum default linearisation energy over all APWs and local-orbitals
real(8) e0min
! if autolinengy is .true. then the fixed linearisation energies are set to the
! Fermi energy minus dlefe
logical autolinengy
! difference between linearisation and Fermi energies when autolinengy is .true.
real(8) dlefe
! lorbcnd is .true. if conduction state local-orbitals should be added
logical lorbcnd
! conduction state local-orbital order
integer lorbordc
! excess order of the APW and local-orbital functions
integer nxoapwlo
! excess local orbitals
integer nxlo

!-------------------------------------------!
!     overlap and Hamiltonian variables     !
!-------------------------------------------!
! overlap and Hamiltonian matrices sizes at each k-point
integer, allocatable :: nmat(:,:)
! maximum nmat over all k-points
integer nmatmax
! index to the position of the local-orbitals in the H and O matrices
integer, allocatable :: idxlo(:,:,:)
! APW-local-orbital overlap integrals
real(8), allocatable :: oalo(:,:,:)
! local-orbital-local-orbital overlap integrals
real(8), allocatable :: ololo(:,:,:)
! APW-APW Hamiltonian integrals
real(8), allocatable :: haa(:,:,:,:,:,:)
! local-orbital-APW Hamiltonian integrals
real(8), allocatable :: hloa(:,:,:,:,:)
! local-orbital-local-orbital Hamiltonian integrals
real(8), allocatable :: hlolo(:,:,:,:)
! complex Gaunt coefficient array
complex(8), allocatable :: gntyry(:,:,:)
! tefvr is .true. if the first-variational eigenvalue equation is to be solved
! as a real symmetric problem
logical tefvr
! tefvit is .true. if the first-variational eigenvalue equation is to be solved
! iteratively
logical tefvit
! minimum and maximum allowed number of eigenvalue equation iterations
integer minitefv,maxitefv
! eigenvalue mixing parameter for iterative solver
real(8) befvit
! iterative solver convergence tolerance
real(8) epsefvit
! type of eigenvalue solver to be used
integer evtype

!--------------------------------------------!
!     eigenvalue and occupancy variables     !
!--------------------------------------------!
! number of empty states per atom and spin
real(8) nempty0
! number of empty states
integer nempty
! number of first-variational states
integer nstfv
! number of second-variational states
integer nstsv
! smearing type
integer stype
! smearing function description
character(256) sdescr
! smearing width
real(8) swidth
! autoswidth is .true. if the smearing width is to be determined automatically
logical autoswidth
! effective mass used in smearing width formula
real(8) mstar
! maximum allowed occupancy (1 or 2)
real(8) occmax
! convergence tolerance for occupancies
real(8) epsocc
! second-variational occupation numbers
real(8), allocatable :: occsv(:,:)
! Fermi energy for second-variational states
real(8) efermi
! scissor correction applied when computing response functions
real(8) scissor
! density of states at the Fermi energy
real(8) fermidos
! estimated indirect and direct band gaps
real(8) bandgap(2)
! k-points of indirect and direct gaps
integer ikgap(3)
! error tolerance for the first-variational eigenvalues
real(8) evaltol
! second-variational eigenvalues
real(8), allocatable :: evalsv(:,:)
! tevecsv is .true. if second-variational eigenvectors are calculated
logical tevecsv
! maximum number of k-point and states indices in user-defined list
integer, parameter :: maxkst=20
! number of k-point and states indices in user-defined list
integer nkstlist
! user-defined list of k-point and state indices
integer kstlist(2,maxkst)

!------------------------------!
!     core state variables     !
!------------------------------!
! occupancies for core states
real(8), allocatable :: occcr(:,:)
! eigenvalues for core states
real(8), allocatable :: evalcr(:,:)
! radial wavefunctions for core states
real(8), allocatable :: rwfcr(:,:,:,:)
! radial charge density for core states
real(8), allocatable :: rhocr(:,:,:)
! spincore is .true. if the core is to be treated as spin-polarised
logical spincore
! number of core spin-channels
integer nspncr

!--------------------------!
!     energy variables     !
!--------------------------!
! eigenvalue sum
real(8) evalsum
! electron kinetic energy
real(8) engykn
! core electron kinetic energy
real(8) engykncr
! nuclear-nuclear energy
real(8) engynn
! electron-nuclear energy
real(8) engyen
! Hartree energy
real(8) engyhar
! Coulomb energy (E_nn + E_en + E_H)
real(8) engycl
! electronic Coulomb potential energy
real(8) engyvcl
! Madelung term
real(8) engymad
! exchange-correlation potential energy
real(8) engyvxc
! exchange-correlation effective field energy
real(8) engybxc
! energy of external global magnetic field
real(8) engybext
! exchange energy
real(8) engyx
! correlation energy
real(8) engyc
! electronic entropy
real(8) entrpy
! entropic contribution to free energy
real(8) engyts
! total energy
real(8) engytot

!------------------------------------!
!     force and stress variables     !
!------------------------------------!
! tforce is .true. if force should be calculated
logical tforce
! Hellmann-Feynman force on each atom
real(8), allocatable :: forcehf(:,:)
! incomplete basis set (IBS) force on each atom
real(8), allocatable :: forceibs(:,:)
! total force on each atom
real(8), allocatable :: forcetot(:,:)
! previous total force on each atom
real(8), allocatable :: forcetotp(:,:)
! maximum force magnitude over all atoms
real(8) forcemax
! tfav0 is .true. if the average force should be zero in order to prevent
! translation of the atomic basis
logical tfav0
! average force
real(8) forceav(3)
! atomic position optimisation type
!  0 : no optimisation
!  1 : unconstrained optimisation
integer atpopt
! maximum number of atomic position optimisation steps
integer maxatpstp
! default step size parameter for atomic position optimisation
real(8) tau0atp
! step size parameters for each atom
real(8), allocatable :: tauatp(:)
! number of strain tensors
integer nstrain
! current strain tensor
integer istrain
data istrain / 0 /
! strain tensors
real(8) strain(3,3,9)
! infinitesimal displacement parameter multiplied by the strain tensor for
! computing the stress tensor
real(8) deltast
! symmetry reduced stress tensor components
real(8) stress(9)
! previous stress tensor
real(8) stressp(9)
! stress tensor component magnitude maximum
real(8) stressmax
! lattice vector optimisation type
!  0 : no optimisation
!  1 : unconstrained optimisation
!  2 : iso-volumetric optimisation
integer latvopt
! maximum number of lattice vector optimisation steps
integer maxlatvstp
! default step size parameter for lattice vector optimisation
real(8) tau0latv
! step size for each stress tensor component acting on the lattice vectors
real(8) taulatv(9)

!-------------------------------!
!     convergence variables     !
!-------------------------------!
! maximum number of self-consistent loops
integer maxscl
! current self-consistent loop number
integer iscl
! Kohn-Sham potential convergence tolerance
real(8) epspot
! energy convergence tolerance
real(8) epsengy
! force convergence tolerance
real(8) epsforce
! stress tensor convergence tolerance
real(8) epsstress

!----------------------------------------------------------!
!     density of states, optics and response variables     !
!----------------------------------------------------------!
! number of energy intervals in the DOS/optics function plot
integer nwplot
! fine k-point grid size for integration of functions in the Brillouin zone
integer ngrkf
! smoothing level for DOS/optics function plot
integer nswplot
! energy interval for DOS/optics function plot
real(8) wplot(2)
! maximum angular momentum for the partial DOS plot
integer lmaxdos
! dosocc is .true. if the DOS is to be weighted by the occupancy
logical dosocc
! dosmsum is .true. if the partial DOS is to be summed over m
logical dosmsum
! dosssum is .true. if the partial DOS is to be summed over spin
logical dosssum
! number of optical matrix components required
integer noptcomp
! required optical matrix components
integer optcomp(3,27)
! intraband is .true. if the intraband term is to be added to the optical matrix
logical intraband
! lmirep is .true. if the (l,m) band characters should correspond to the
! irreducible representations of the site symmetries
logical lmirep
! spin-quantisation axis in Cartesian coordinates used when plotting the
! spin-resolved DOS (z-axis by default)
real(8) sqados(3)
! q-vector in lattice and Cartesian coordinates for calculating the matrix
! elements < i,k+q | exp(iq.r) | j,k >
real(8) vecql(3),vecqc(3)
! maximum initial-state energy allowed in ELNES transitions
real(8) emaxelnes
! structure factor energy window
real(8) wsfac(2)

!-------------------------------------!
!     1D/2D/3D plotting variables     !
!-------------------------------------!
! number of vertices in 1D plot
integer nvp1d
! total number of points in 1D plot
integer npp1d
! vertices in lattice coordinates for 1D plot
real(8), allocatable :: vvlp1d(:,:)
! distance to vertices in 1D plot
real(8), allocatable :: dvp1d(:)
! plot vectors in lattice coordinates for 1D plot
real(8), allocatable :: vplp1d(:,:)
! distance to points in 1D plot
real(8), allocatable :: dpp1d(:)
! corner vectors of 2D plot in lattice coordinates
real(8) vclp2d(3,0:2)
! grid sizes of 2D plot
integer np2d(2)
! corner vectors of 3D plot in lattice coordinates
real(8) vclp3d(3,0:3)
! grid sizes of 3D plot
integer np3d(3)

!----------------------------------------!
!     OEP and Hartree-Fock variables     !
!----------------------------------------!
! maximum number of core states over all species
integer ncrmax
! maximum number of OEP iterations
integer maxitoep
! initial value and scaling factors for OEP step size
real(8) tauoep(3)
! magnitude of the OEP residual
real(8) resoep
! exchange potential and magnetic field
real(8), allocatable :: vxmt(:,:),vxir(:)
real(8), allocatable :: bxmt(:,:,:),bxir(:,:)
! hybrid is .true. if a hybrid functional is to be used
logical hybrid
! hybrid functional mixing coefficient
real(8) hybridc

!-------------------------------------------------------------!
!     response function and perturbation theory variables     !
!-------------------------------------------------------------!
! |G| cut-off for response functions
real(8) gmaxrf
! energy cut-off for response functions
real(8) emaxrf
! number of G-vectors for response functions
integer ngrf
! number of response function frequencies
integer nwrf
! complex response function frequencies
complex(8), allocatable :: wrf(:)
! maximum number of spherical Bessel functions on the coarse radial mesh over
! all species
integer njcmax

!-------------------------------------------------!
!     Bethe-Salpeter equation (BSE) variables     !
!-------------------------------------------------!
! number of valence and conduction states for transitions
integer nvbse,ncbse
! default number of valence and conduction states
integer nvbse0,ncbse0
! maximum number of extra valence and conduction states
integer, parameter :: maxxbse=20
! number of extra valence and conduction states
integer nvxbse,ncxbse
! extra valence and conduction states
integer istxbse(maxxbse),jstxbse(maxxbse)
! total number of transitions
integer nvcbse
! size of blocks in BSE Hamiltonian matrix
integer nbbse
! size of BSE matrix (= 2*nbbse)
integer nmbse
! index from BSE valence states to second-variational states
integer, allocatable :: istbse(:,:)
! index from BSE conduction states to second-variational states
integer, allocatable :: jstbse(:,:)
! index from BSE valence-conduction pair and k-point to location in BSE matrix
integer, allocatable :: ijkbse(:,:,:)
! BSE Hamiltonian
complex(8), allocatable :: hmlbse(:,:)
! BSE Hamiltonian eigenvalues
real(8), allocatable :: evalbse(:)
! if bsefull is .true. then the full BSE Hamiltonian is calculated, otherwise
! only the Hermitian block
logical bsefull
! if hxbse/hdbse is .true. then the exchange/direct term is included in the BSE
! Hamiltonian
logical hxbse,hdbse

!--------------------------!
!     timing variables     !
!--------------------------!
! initialisation
real(8) timeinit
! Hamiltonian and overlap matrix set up
real(8) timemat
! first-variational calculation
real(8) timefv
! second-variational calculation
real(8) timesv
! charge density calculation
real(8) timerho
! potential calculation
real(8) timepot
! force calculation
real(8) timefor

!-----------------------------!
!     numerical constants     !
!-----------------------------!
real(8), parameter :: pi=3.1415926535897932385d0
real(8), parameter :: twopi=6.2831853071795864769d0
real(8), parameter :: fourpi=12.566370614359172954d0
! spherical harmonic for l=m=0
real(8), parameter :: y00=0.28209479177387814347d0
! complex constants
complex(8), parameter :: zzero=(0.d0,0.d0)
complex(8), parameter :: zone=(1.d0,0.d0)
complex(8), parameter :: zi=(0.d0,1.d0)
! array of i^l and (-i)^l values
complex(8), allocatable :: zil(:),zilc(:)
! Pauli spin matrices:
! sigma_x = ( 0  1 )   sigma_y = ( 0 -i )   sigma_z = ( 1  0 )
!           ( 1  0 )             ( i  0 )             ( 0 -1 )
complex(8) sigmat(2,2,3)
data sigmat / (0.d0,0.d0), (1.d0,0.d0), (1.d0,0.d0), (0.d0,0.d0), &
              (0.d0,0.d0), (0.d0,1.d0),(0.d0,-1.d0), (0.d0,0.d0), &
              (1.d0,0.d0), (0.d0,0.d0), (0.d0,0.d0),(-1.d0,0.d0) /
! Hartree in eV (CODATA 2010)
real(8), parameter :: ha_ev=27.21138505d0
! Boltzmann constant in eV/kelvin (CODATA 2010)
real(8), parameter :: kb_ev=8.6173324d-5
! Boltzmann constant in Hartree/kelvin
real(8), parameter :: kboltz=kb_ev/ha_ev
! speed of light in atomic units (=1/alpha) (CODATA 2010)
real(8), parameter :: sol=137.035999074d0
! scaled speed of light
real(8) solsc
! electron g-factor (CODATA 2010)
real(8), parameter :: gfacte=2.00231930436153d0
! Hartree in SI units (CODATA 2010)
real(8), parameter :: ha_si=4.35974434d-18
! Hartree in inverse meters (CODATA 2010)
real(8), parameter :: ha_im=2.194746313708d7
! Bohr radius in SI units (CODATA 2010)
real(8), parameter :: br_si=0.52917721092d-10
! Planck constant in SI units (CODATA 2010)
real(8), parameter :: hbar_si=1.054571726d-34
! electron charge in SI units (CODATA 2010)
real(8), parameter :: e_si=1.602176565d-19
! atomic unit of magnetic flux density in SI
real(8), parameter :: b_si=hbar_si/(e_si*br_si**2)
! atomic unit of electric field in SI
real(8), parameter :: ef_si=ha_si/(e_si*br_si)
! atomic unit of time in SI
real(8), parameter :: t_si=hbar_si/ha_si
! electron mass in SI (CODATA 2010)
real(8), parameter :: em_si=9.10938291d-31
! atomic mass unit in SI (CODATA 2010)
real(8), parameter :: amu_si=1.660538921d-27
! atomic mass unit in electron masses
real(8), parameter :: amu=amu_si/em_si

!---------------------------------!
!     miscellaneous variables     !
!---------------------------------!
! code version
integer version(3)
data version / 5,2,14 /
! maximum number of tasks
integer, parameter :: maxtasks=40
! number of tasks
integer ntasks
! task array
integer tasks(maxtasks)
! current task
integer task
! tlast is .true. if the calculation is on the last self-consistent loop
logical tlast
! tstop is .true. if the STOP file exists
logical tstop
! number of self-consistent loops after which STATE.OUT is written
integer nwrite
! if wrtvars is .true. then variables are written to VARIABLES.OUT
logical wrtvars
! filename extension for files generated by gndstate
character(256) filext
! default file extension
data filext / '.OUT' /
! scratch space path
character(256) scrpath
! maximum number of note lines
integer, parameter :: maxnlns=30
! number of note lines
integer notelns
! notes to include in INFO.OUT
character(80) notes(maxnlns)

end module

