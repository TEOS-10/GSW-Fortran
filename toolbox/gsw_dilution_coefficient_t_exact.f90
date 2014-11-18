!==========================================================================
elemental function gsw_dilution_coefficient_t_exact (sa, t, p)
!==========================================================================
!
!  Calculates the dilution coefficient of seawater.  The dilution 
!  coefficient of seawater is defined as the Absolute Salinity times the 
!  second derivative of the Gibbs function with respect to Absolute 
!  Salinity, that is, SA.*g_SA_SA.
!
!  SA =  Absolute Salinity                                         [ g/kg ]
!  t  =  in-situ temperature (ITS-90)                             [ deg C ]
!  p  =  sea pressure                                              [ dbar ]
!        ( i.e. absolute pressure - 10.1325 dbar ) 
!
!  dilution_coefficient_t_exact  =  dilution coefficient   [ (J/kg)(kg/g) ]
!--------------------------------------------------------------------------

use gsw_mod_teos10_constants, only : gsw_sfac

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, t, p

real (r14) :: gsw_dilution_coefficient_t_exact

real (r14) :: g08, x, x2, y, z

x2 = gsw_sfac*sa
x = sqrt(x2)
y = t*0.025d0
z = p*1d-4      ! note.the input pressure (p) is sea pressure in units of dbar.

g08 = 2.0d0*(8103.20462414788d0 + &
          y*(2175.341332000392d0 + &
	      y*(-274.2290036817964d0 + &
                  y*(197.4670779425016d0 + &
		      y*(-68.5590309679152d0 + 9.98788038278032d0*y))) - &
          90.6734234051316d0*z) + &
	      1.5d0*x*(-5458.34205214835d0 - 980.14153344888d0*y + &
                  (4.0d0/3.0d0)*x*(2247.60742726704d0 - &
		  340.1237483177863d0*1.25d0*x + 220.542973797483d0*y) + &
              180.142097805543d0*z) + &
          z*(-219.1676534131548d0 + &
	      (-16.32775915649044d0 - 120.7020447884644d0*z)*z))
    
g08 = x2*g08 + & 
          x*(-7296.43987145382d0 + &
	      z*(598.378809221703d0 + &
                  z*(-156.8822727844005d0 + &
		      (204.1334828179377d0 - 10.23755797323846d0*z)*z)) + &
              y*(-1480.222530425046d0 + &
	          z*(-525.876123559641d0 + &
                      (249.57717834054571d0 - 88.449193048287d0*z)*z) + &
                  y*(-129.1994027934126d0 + &
	              z*(1149.174198007428d0 + &
                          z*(-162.5751787551336d0 + 76.9195462169742d0*z)) + &
                  y*(-30.0682112585625d0 - 1380.9597954037708d0*z + &
                      y*(2.626801985426835d0 + 703.695562834065d0*z))))) + &
      11625.62913253464d0 + 1702.453469893412d0*y
    
gsw_dilution_coefficient_t_exact = 0.25d0*gsw_sfac*g08

! Note that this function avoids the singularity that occurs at SA = 0 if
! the straightforward expression for the dilution coefficient of seawater,
! SA*g_SA_SA is simply evaluated as SA.*gsw_gibbs(2,0,0,SA,t,p). 

return
end function

!--------------------------------------------------------------------------

