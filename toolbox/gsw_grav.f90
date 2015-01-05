!==========================================================================
elemental function gsw_grav (lat, p)  
!==========================================================================
!
! Calculates acceleration due to gravity as a function of latitude and as
!  a function of pressure in the ocean.
!
! lat  =  latitude in decimal degress north                [ -90 ... +90 ]  
!  p  =  sea pressure                                              [ dbar ]
! 
! gsw_grav : grav  =  gravitational acceleration               [ m s^-2 ]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_z_from_p

use gsw_mod_teos10_constants, only : deg2rad, gamma

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: lat, p  

real (r14) :: gsw_grav

real (r14) :: sin2, gs, z

sin2 = sin(lat*deg2rad)**2
gs = 9.780327d0*(1d0 + (5.2792d-3 + (2.32d-5*sin2))*sin2) 

z = gsw_z_from_p(p,lat)

gsw_grav = gs*(1d0 - gamma*z)           ! z is the height corresponding to p. 
                                        ! Note. In the ocean z is negative.
return
end function

!--------------------------------------------------------------------------
