!==========================================================================
elemental function gsw_melting_seaice_equilibrium_sa_ct_ratio (sa, p, &
                                                       saturation_fraction)
!==========================================================================
!
!  Calculates the ratio of SA to CT changes when sea ice melts into 
!  seawater with both the seawater and the sea ice temperatures being  
!  almost equal to the equilibrium freezing temperature.  It is assumed  
!  that a small mass of seaice melts into an infinite mass of seawater.  If 
!  indeed the temperature of the seawater and the sea ice were both equal  
!  to the freezing temperature, then no melting or freezing would occur; an  
!  imbalance between these three temperatures is needed for freezing or 
!  melting to occur (the three temperatures being (1) the seawater 
!  temperature, (2) the sea ice temperature, and (3) the freezing 
!  temperature.  
!
!  Note that the output of this function, dSA/dCT is independent of the 
!  sea ice salinity, SA_seaice.  That is, the output applies equally to  
!  pure ice Ih and to sea ice with seaice salinity, SA_seaice.  This result 
!  is proven in the manuscript, McDougall et al. (2013).  
!
!  The output, melting_seaice_equilibrium_SA_CT_ratio, is dSA/dCT rather  
!  than dCT/dSA.  This is done so that when SA = 0, the output, dSA/dCT is 
!  zero whereas dCT/dSA would be infinite. 
!
!  SA  =  Absolute Salinity of seawater                            [ g/kg ]
!  p   =  sea pressure at which the melting occurs                 [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar ) 
!  saturation_fraction = the saturation fraction of dissolved air in 
!               seawater.  The saturation_fraction must be between 0 and 1.
!
!  melting_seaice_equilibrium_SA_CT_ratio = the ratio dSA/dCT of SA to CT  
!                            changes when sea ice melts into seawater, with   
!                            the seawater and sea ice being close to the  
!                            freezing temperature.             [ g/(kg K) ] 
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_ct_freezing, gsw_enthalpy, gsw_t_freezing
use gsw_mod_toolbox, only : gsw_enthalpy_ice, gsw_enthalpy_first_derivatives

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, p, saturation_fraction

real (r14) :: gsw_melting_seaice_equilibrium_sa_ct_ratio

real (r14) :: ctf, denominator, h, h_ih, t_seaice, h_hat_sa, h_hat_ct

ctf = gsw_ct_freezing(sa,p,saturation_fraction)
t_seaice = gsw_t_freezing(sa,p,saturation_fraction)

h = gsw_enthalpy(sa,ctf,p)
h_ih = gsw_enthalpy_ice(t_seaice,p)
call gsw_enthalpy_first_derivatives(sa,ctf,p,h_hat_sa,h_hat_ct)
          ! note that h_hat_ct is equal to cp0*(273.15 + t)/(273.15 + pt0)

denominator = h - h_ih - sa*h_hat_sa
gsw_melting_seaice_equilibrium_sa_ct_ratio = sa*h_hat_ct/denominator
         
return
end function

!--------------------------------------------------------------------------
