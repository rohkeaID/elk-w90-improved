
subroutine wannier_run(seed__name,mp_grid_loc,num_kpts_loc, &
     real_lattice_loc,recip_lattice_loc,kpt_latt_loc,num_bands_loc, &
     num_wann_loc,nntot_loc,num_atoms_loc,atom_symbols_loc, &
     atoms_cart_loc,gamma_only_loc,M_matrix_loc,A_matrix_loc,eigenvalues_loc, &
     U_matrix_loc,U_matrix_opt_loc,lwindow_loc,wann_centres_loc, &
     wann_spreads_loc,spread_loc)

  implicit none

  character(len=*),     intent(in   ) :: seed__name
  integer,              intent(in   ) :: mp_grid_loc(3)
  integer,              intent(in   ) :: num_kpts_loc
  real(8),              intent(in   ) :: real_lattice_loc(3,3)
  real(8),              intent(in   ) :: recip_lattice_loc(3,3)
  real(8),              intent(in   ) :: kpt_latt_loc(3,num_kpts_loc)
  integer,              intent(in   ) :: num_bands_loc
  integer,              intent(in   ) :: num_wann_loc
  integer,              intent(in   ) :: nntot_loc
  integer,              intent(in   ) :: num_atoms_loc
  character(len=*),     intent(in   ) :: atom_symbols_loc(num_atoms_loc)
  real(8),              intent(in   ) :: atoms_cart_loc(3,num_atoms_loc)
  logical,              intent(in   ) :: gamma_only_loc
  complex(8),           intent(in   ) :: M_matrix_loc(num_bands_loc,num_bands_loc,nntot_loc,num_kpts_loc)
  complex(8),           intent(in   ) :: A_matrix_loc(num_bands_loc,num_wann_loc,num_kpts_loc)
  real(8),              intent(in   ) :: eigenvalues_loc(num_bands_loc,num_kpts_loc)
  complex(8),           intent(  out) :: U_matrix_loc(num_wann_loc,num_wann_loc,num_kpts_loc)
  complex(8), optional, intent(  out) :: U_matrix_opt_loc(num_bands_loc,num_wann_loc,num_kpts_loc)
  logical,    optional, intent(  out) :: lwindow_loc(num_bands_loc,num_kpts_loc)
  real(8),    optional, intent(  out) :: wann_centres_loc(3,num_wann_loc)
  real(8),    optional, intent(  out) :: wann_spreads_loc(num_wann_loc)
  real(8),    optional, intent(  out) :: spread_loc(3)
!-------------------------------------------------------------------------------
  write(*,*)
  write(*,'("Error(libwannier): libwannier not or improperly installed")')
  write(*,'("                   You cannot use tasks 602-605 without Wannier90 library")')
  write(*,*)

  stop

end subroutine wannier_run
