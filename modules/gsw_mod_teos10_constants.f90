!==========================================================================
module gsw_mod_teos10_constants
!==========================================================================

implicit none

integer, parameter :: gtc_r14 = selected_real_kind(14,30)

real (gtc_r14), parameter :: db2pa = 1d4
real (gtc_r14), parameter :: rec_db2pa = 1d-4

real (gtc_r14), parameter :: pa2db = 1d-4
real (gtc_r14), parameter :: rec_pa2db = 1d4

real (gtc_r14), parameter :: pi = 3.141592653589793d0

!  cp0  =  The "specific heat" for use                         [ J/(kg K) ]
!          with Conservative Temperature   

real (gtc_r14), parameter :: gsw_cp0 = 3991.86795711963d0

!  T0  =  the Celcius zero point.                                     [ K ]

real (gtc_r14), parameter :: gsw_t0 = 273.15d0

!  P0  =  Absolute Pressure of one standard atmosphere.              [ Pa ]

real (gtc_r14), parameter :: gsw_p0 = 101325d0

!  SSO  =  Standard Ocean Reference Salinity.                      [ g/kg ]

real (gtc_r14), parameter :: gsw_sso = 35.16504d0
real (gtc_r14), parameter :: gsw_sqrtsso = 5.930011804372737d0

!  uPS  =  unit conversion factor for salinities                   [ g/kg ]

real (gtc_r14), parameter :: gsw_ups = gsw_sso/35.d0

!  sfac  =  1/(40*gsw_ups)

real (gtc_r14), parameter :: gsw_sfac = 0.0248826675584615d0

!  C3515  =  Conductivity at (SP=35, t_68=15, p=0)                [ mS/cm ]

real (gtc_r14), parameter :: gsw_c3515 = 42.9140d0

!  SonCl  =  SP to Chlorinity ratio                           [ (g/kg)^-1 ]

real (gtc_r14), parameter :: gsw_soncl = 1.80655d0

!  valence_factor  =  valence factor of sea salt of Reference Composition
!                                                              [ unitless ]

real (gtc_r14), parameter :: gsw_valence_factor = 1.2452898d0

!  atomic_weight = mole-weighted atomic weight of sea salt of Reference 
!                  Composition                                    [ g/mol ]

real (gtc_r14), parameter :: gsw_atomic_weight = 31.4038218d0

end module

!--------------------------------------------------------------------------
