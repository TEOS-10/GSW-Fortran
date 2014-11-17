!==========================================================================
function gsw_grav(lat,p)  
!==========================================================================

! Calculates acceleration due to gravity as a function of latitude and as
!  a function of pressure in the ocean.
!
! lat  =  latitude in decimal degress north                [ -90 ... +90 ]  
!  p  =  sea pressure                                              [ dbar ]
! 
! gsw_grav : grav  =  gravitational acceleration               [ m s^-2 ]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), parameter :: pi = 3.141592653589793d0

real (r14) :: lat, p, gsw_grav, gamma, deg2rad, x, sin2, gs, z, gsw_z_from_p

gamma = 2.26d-7
deg2rad = pi/180d0
x = sin(lat*deg2rad)  ! convert to radians
sin2 = x*x
gs = 9.780327d0*(1.0d0 + (5.2792d-3 + (2.32d-5*sin2))*sin2) 

z = gsw_z_from_p(p,lat)

gsw_grav = gs*(1 - gamma*z)             ! z is the height corresponding to p. 
                                        ! Note. In the ocean z is negative.

return
end

!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! basic thermodynamic properties in terms of in-situ t, based on the exact Gibbs function
!--------------------------------------------------------------------------

