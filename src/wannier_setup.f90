
subroutine wannier_setup(seed__name,mp_grid_loc,num_kpts_loc,&
     real_lattice_loc,recip_lattice_loc,kpt_latt_loc,num_bands_tot, &
     num_atoms_loc,atom_symbols_loc,atoms_cart_loc, gamma_only_loc,spinors_loc, &
     nntot_loc,nnlist_loc,nncell_loc,num_bands_loc,num_wann_loc, &
     proj_site_loc,proj_l_loc,proj_m_loc,proj_radial_loc,proj_z_loc, &
     proj_x_loc,proj_zona_loc,exclude_bands_loc,proj_s_loc,proj_s_qaxis_loc)

  implicit none

  character(len=*),  intent(in   ) :: seed__name
  integer,           intent(in   ) :: mp_grid_loc(3)
  integer,           intent(in   ) :: num_kpts_loc
  real(8),           intent(in   ) :: real_lattice_loc(3,3)
  real(8),           intent(in   ) :: recip_lattice_loc(3,3)
  real(8),           intent(in   ) :: kpt_latt_loc(3,num_kpts_loc)
  integer,           intent(in   ) :: num_bands_tot
  integer,           intent(in   ) :: num_atoms_loc
  character(len=*),  intent(in   ) :: atom_symbols_loc(num_atoms_loc)
  real(8),           intent(in   ) :: atoms_cart_loc(3,num_atoms_loc)
  logical,           intent(in   ) :: gamma_only_loc
  logical,           intent(in   ) :: spinors_loc
  integer,           intent(  out) :: nntot_loc
  integer,           intent(  out) :: nnlist_loc(num_kpts_loc,12)  ! AG: num_nnmax=12
  integer,           intent(  out) :: nncell_loc(3,num_kpts_loc,12)! AG: num_nnmax=12
  integer,           intent(  out) :: num_bands_loc
  integer,           intent(  out) :: num_wann_loc
  real(8),           intent(  out) :: proj_site_loc(3,num_bands_tot)
  integer,           intent(  out) :: proj_l_loc(num_bands_tot)
  integer,           intent(  out) :: proj_m_loc(num_bands_tot)
  integer,           intent(  out) :: proj_radial_loc(num_bands_tot)
  real(8),           intent(  out) :: proj_z_loc(3,num_bands_tot)
  real(8),           intent(  out) :: proj_x_loc(3,num_bands_tot)
  real(8),           intent(  out) :: proj_zona_loc(num_bands_tot)
  integer,           intent(  out) :: exclude_bands_loc(num_bands_tot)
  integer, optional, intent(  out) :: proj_s_loc(num_bands_tot)
  real(8), optional, intent(  out) :: proj_s_qaxis_loc(3,num_bands_tot)
!-------------------------------------------------------------------------------
  write(*,*)
  write(*,'("Error(libwannier): libwannier not or improperly installed")')
  write(*,'("                   You cannot use tasks 602-605 without Wannier90 library")')
  write(*,*)

  stop

end subroutine wannier_setup
