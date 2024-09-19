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

module atm_ice_wave_exchange_mod

  !! FMS
  use FMS

  use atmos_model_mod,     only: atmos_data_type
  use land_model_mod,      only: land_data_type
  use ice_model_mod,       only: ice_data_type
  use ocean_model_mod,     only: ocean_public_type
  use wave_model_mod,      only: wave_data_type, atmos_wave_boundary_type, ice_wave_boundary_type

  implicit none
  private


  !---- exchange grid maps -----

  type(FmsXgridXmap_type), save :: xmap_atm_wav
  integer, save         :: n_xgrid_atm_wav
  type(FmsXgridXmap_type), save :: xmap_ice_wav
  integer, save         :: n_xgrid_ice_wav

  ! Exchange grid indices
  integer :: X2_GRID_ATM, X2_GRID_WAV

  public :: atm_wave_exchange_init, atm_to_wave, ice_wave_exchange_init, ice_to_wave

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

    call fms_xgrid_setup_xmap(xmap_atm_wav, (/ 'ATM', 'WAV' /),       &
         (/ atm%Domain, Wav%Domain /),                    &
         "INPUT/grid_spec.nc", Atm%grid )
    ! exchange grid indices
    X2_GRID_ATM = 1; X2_GRID_WAV = 2;
    n_xgrid_atm_wav = max(fms_xgrid_count(xmap_atm_wav),1)
    call fms_mpp_domains_get_compute_domain( Wav%domain, is, ie, js, je )

    !allocate land_ice_boundary
    allocate( atmos_wave_boundary%wavgrd_u10_mpp(is:ie,js:je,1) )
    allocate( atmos_wave_boundary%wavgrd_v10_mpp(is:ie,js:je,1) )

    atmos_wave_boundary%wavgrd_u10_mpp(:,:,:) = 0.0
    atmos_wave_boundary%wavgrd_v10_mpp(:,:,:) = 0.0

  end subroutine atm_wave_exchange_init

  subroutine ice_wave_exchange_init(Ice, Wav, Ice_wave_boundary)
    type(ice_data_type),          intent(in)      :: Ice !< A derived data type to specify ocean/ice boundary data
    type(wave_data_type),           intent(inout) :: Wav !< A derived data type to specify wave boundary data
    type(ice_wave_boundary_type), intent(inout)   :: Ice_wave_boundary !< A derived data type to specify properties
                                                                     !! passed from atmos to waves

    integer :: is, ie, js, je

    call fms_xgrid_setup_xmap(xmap_ice_wav, (/ 'WAV', 'OCN' /),       &
         (/ Wav%Domain, Ice%Domain /),                    &
         "INPUT/grid_spec.nc" )
    ! exchange grid indices
    n_xgrid_ice_wav = max(fms_xgrid_count(xmap_ice_wav),1)
    call fms_mpp_domains_get_compute_domain( Wav%domain, is, ie, js, je )

    !allocate land_ice_boundary
    allocate( ice_wave_boundary%wavgrd_ucurr_mpp(is:ie,js:je,1) )
    ice_wave_boundary%wavgrd_ucurr_mpp(:,:,:) = 0.0
    allocate( ice_wave_boundary%wavgrd_vcurr_mpp(is:ie,js:je,1) )
    ice_wave_boundary%wavgrd_vcurr_mpp(:,:,:) = 0.0

    allocate( wav%ustkb_mpp(is:ie,js:je,wav%num_stk_bands) )
    wav%ustkb_mpp(:,:,:) = 0.0
    allocate( wav%vstkb_mpp(is:ie,js:je,wav%num_stk_bands) )
    wav%vstkb_mpp(:,:,:) = 0.0

    ! This are a temporary and costly trick to make MPI work
    allocate( wav%glob_loc_X(is:ie,js:je) )
    wav%glob_loc_X(:,:) = 0
    allocate( wav%glob_loc_Y(is:ie,js:je) )
    wav%glob_loc_Y(:,:) = 0

    call fms_mpp_domains_get_compute_domain( Ice%domain, is, ie, js, je )
    allocate( ice_wave_boundary%icegrd_ustkb_mpp(is:ie,js:je,1,wav%num_stk_bands) )
    ice_wave_boundary%icegrd_ustkb_mpp(:,:,:,:) = 0.0
    allocate( ice_wave_boundary%icegrd_vstkb_mpp(is:ie,js:je,1,wav%num_stk_bands) )
    ice_wave_boundary%icegrd_vstkb_mpp(:,:,:,:) = 0.0

    return
  end subroutine ice_wave_exchange_init

  !> Does atmosphere TO wave operations (could do wave TO atmosphere operations too).
  subroutine atm_to_wave( Time, Atm, Wav, Atmos_wave_Boundary )
    type(FmsTime_type),                intent(in) :: Time !< Current time
    type(atmos_data_type),           intent(in) :: Atm
    type(wave_data_type),            intent(in) :: Wav
    type(atmos_wave_boundary_type), intent(inout):: Atmos_wave_Boundary

    real, dimension(n_xgrid_atm_wav) :: &
         ex_u_atm, &
         ex_v_atm

    integer :: remap_method

    remap_method = 1

    call fms_xgrid_put_to_xgrid (Atm%u_bot , 'ATM', ex_u_atm , xmap_atm_wav, remap_method=remap_method, complete=.false.)
    call fms_xgrid_put_to_xgrid (Atm%v_bot , 'ATM', ex_v_atm , xmap_atm_wav, remap_method=remap_method, complete=.true.)
    if (Wav%pe) then
       call fms_xgrid_get_from_xgrid(Atmos_Wave_Boundary%wavgrd_u10_mpp, 'WAV', ex_u_atm, xmap_atm_wav)
       call fms_xgrid_get_from_xgrid(Atmos_Wave_Boundary%wavgrd_v10_mpp, 'WAV', ex_v_atm, xmap_atm_wav)
    endif

  end subroutine atm_to_wave

  !> Does both ice TO wave and wave TO ice exchange grid operations.
  subroutine ice_to_wave( Time, Ice, Wav, Ice_wave_Boundary )
    type(FmsTime_type),                intent(in) :: Time !< Current time
    type(ice_data_type),            intent(in) :: Ice !< The ice module container
    type(wave_data_type),           intent(in) :: Wav !< The wave module container
    type(ice_wave_boundary_type), intent(inout):: Ice_wave_Boundary !< The ice-wave boundary container

    real, dimension(n_xgrid_ice_wav) :: &
         ex_ucurr,  & ! Exchange grid x-current
         ex_vcurr,  & ! Exchange grid y-current
         ex_ustokes,& ! Exchange grid x-Stokes drift
         ex_vstokes   ! Exchange grid y-Stokes drift

    integer :: remap_method ! Interpolation method (todo: list options)
    integer :: i_stk

    remap_method = 1

    ! -> Put Ocean (ice) parameters onto exchange grid
    call fms_xgrid_put_to_xgrid (Ice%u_surf(:,:,:) , 'OCN', ex_ucurr , xmap_ice_wav)
    call fms_xgrid_put_to_xgrid (Ice%v_surf(:,:,:) , 'OCN', ex_vcurr , xmap_ice_wav)

    ! -> Only on wave-PEs, bring wave information off exchange grid
    if (Wav%pe) then
       call fms_xgrid_get_from_xgrid(Ice_Wave_Boundary%wavgrd_ucurr_mpp(:,:,1), 'WAV', ex_ucurr, xmap_ice_wav)
       call fms_xgrid_get_from_xgrid(Ice_Wave_Boundary%wavgrd_vcurr_mpp(:,:,1), 'WAV', ex_vcurr, xmap_ice_wav)
    endif

      ! -> Put Wave parameters onto exchange grid
    do i_stk = 1,wav%num_Stk_bands
      call fms_xgrid_put_to_xgrid (Wav%ustkb_mpp(:,:,i_stk) , 'WAV', ex_ustokes , xmap_ice_wav)
      call fms_xgrid_put_to_xgrid (Wav%vstkb_mpp(:,:,i_stk) , 'WAV', ex_vstokes , xmap_ice_wav)
      ! -> Only on ice-PEs, bring ice information off exchange grid
      if (Ice%pe) then
        call fms_xgrid_get_from_xgrid(Ice_Wave_Boundary%icegrd_ustkb_mpp(:,:,:,i_stk), 'OCN', ex_ustokes, xmap_ice_wav)
        call fms_xgrid_get_from_xgrid(Ice_Wave_Boundary%icegrd_vstkb_mpp(:,:,:,i_stk), 'OCN', ex_vstokes, xmap_ice_wav)
      endif

    enddo

    return
  end subroutine ice_to_wave

end module atm_ice_wave_exchange_mod
