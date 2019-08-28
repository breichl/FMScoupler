!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS) Coupler.
!*
!* FMS Coupler is free software: you can redistribute it and/or modify
!* it under the terms of the GNU Lesser General Public License as
!* published by the Free Software Foundation, either version 3 of the
!* License, or (at your option) any later version.
!*
!* FMS Coupler is distributed in the hope that it will be useful, but
!* WITHOUT ANY WARRANTY; without even the implied warranty of
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!* General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS Coupler.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

module atm_wave_exchange_mod

  use mpp_mod,             only: mpp_pe, mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE
  use mpp_domains_mod,     only: mpp_get_compute_domain
  use fms_mod,             only: clock_flag_default
  use constants_mod,       only: RADIUS
  use xgrid_mod,           only: xmap_type, setup_xmap, xgrid_count, stock_move, &
                                 put_to_xgrid, get_from_xgrid
  use time_manager_mod,    only: time_type
  use data_override_mod,   only: data_override
  use atmos_model_mod,     only: atmos_data_type
  use land_model_mod,     only: land_data_type
  use ocean_model_mod,     only: ocean_public_type
  use wave_model_mod,      only: wave_data_type, atmos_wave_boundary_type
  use stock_constants_mod, only: Lnd_stock, Ice_stock, ISTOCK_WATER, ISTOCK_SIDE

  implicit none
  private


  !---- exchange grid maps -----

  type(xmap_type), save :: xmap_atm_wav
  integer, save         :: n_xgrid_atm_wav

  ! Exchange grid indices
  integer :: X2_GRID_ATM, X2_GRID_WAV

  public :: atm_wave_exchange_init, atm_to_wave

  integer :: cplClock, fluxLandIceClock
  logical :: do_runoff
  real    :: Dt_cpl
contains

  subroutine atm_wave_exchange_init(Atm, Wav, Atmos_wave_boundary)
    type(atmos_data_type),          intent(in)    :: Atm !< A derived data type to specify atmospheric boundary data
    type(wave_data_type),           intent(inout) :: Wav !< A derived data type to specify wave boundary data
    type(atmos_wave_boundary_type), intent(inout) :: atmos_wave_boundary !< A derived data type to specify properties
                                                                     !! passed from atmos to waves

    integer :: is, ie, js, je

    call setup_xmap(xmap_atm_wav, (/ 'ATM', 'WAV' /),       &
         (/ atm%Domain, Wav%Domain /),                    &
         "INPUT/grid_spec.nc", Atm%grid )
    ! exchange grid indices
    X2_GRID_ATM = 1; X2_GRID_WAV = 2;
    n_xgrid_atm_wav = max(xgrid_count(xmap_atm_wav),1)
    call mpp_get_compute_domain( Wav%domain, is, ie, js, je )

    !allocate land_ice_boundary
    allocate( atmos_wave_boundary%u_10_mpp(is:ie,js:je,1) )
    allocate( atmos_wave_boundary%v_10_mpp(is:ie,js:je,1) )

    atmos_wave_boundary%u_10_mpp(:,:,:) = 0.0
    atmos_wave_boundary%v_10_mpp(:,:,:) = 0.0

  end subroutine atm_wave_exchange_init

  subroutine atm_to_wave( Time, Atm, Wave, Atmos_wave_Boundary )
    type(time_type),                intent(in) :: Time !< Current time
    type(atmos_data_type),           intent(in) :: Atm
    type(wave_data_type),            intent(in) :: Wave
    type(atmos_wave_boundary_type), intent(inout):: Atmos_wave_Boundary

    real, dimension(n_xgrid_atm_wav) :: &
         ex_u_atm, &
         ex_v_atm


    integer :: remap_method

    remap_method = 1

    call put_to_xgrid (Atm%u_bot , 'ATM', ex_u_atm , xmap_atm_wav, remap_method=remap_method, complete=.false.)
    call put_to_xgrid (Atm%v_bot , 'ATM', ex_v_atm , xmap_atm_wav, remap_method=remap_method, complete=.true.)
    if (Wave%pe) then
       call get_from_xgrid(Atmos_Wave_Boundary%U_10_mpp, 'WAV', ex_u_atm, xmap_atm_wav)
       call get_from_xgrid(Atmos_Wave_Boundary%V_10_mpp, 'WAV', ex_v_atm, xmap_atm_wav)
    endif

  end subroutine atm_to_wave

end module atm_wave_exchange_mod
