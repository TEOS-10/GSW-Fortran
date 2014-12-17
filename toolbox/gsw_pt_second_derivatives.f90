!==========================================================================
elemental subroutine gsw_pt_second_derivatives (sa, ct, pt_sa_sa, &
                                                pt_sa_ct, pt_ct_ct)
! =========================================================================
!
!  Calculates the following three second-order derivatives of potential 
!  temperature (the regular potential temperature which has a reference 
!  sea pressure of 0 dbar), 
!   (1) pt_SA_SA, the second derivative with respect to Absolute Salinity 
!       at constant Conservative Temperature,
!   (2) pt_SA_CT, the derivative with respect to Conservative Temperature
!       and Absolute Salinity, and
!   (3) pt_CT_CT, the second derivative with respect to Conservative 
!       Temperature at constant Absolute Salinity. 
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!
!  pt_SA_SA  =  The second derivative of potential temperature (the 
!               regular potential temperature which has reference sea 
!               pressure of 0 dbar) with respect to Absolute Salinity 
!               at constant Conservative Temperature.  
!               pt_SA_SA has units of:                     [ K/((g/kg)^2) ]
!  pt_SA_CT  =  The derivative of potential temperature with respect 
!               to Absolute Salinity and Conservative Temperature.   
!               pt_SA_CT has units of:                         [ 1/(g/kg) ]
!  pt_CT_CT  =  The second derivative of potential temperature (the 
!               regular one with p_ref = 0 dbar) with respect to 
!               Conservative Temperature at constant SA.  
!               pt_CT_CT has units of:                              [ 1/K ]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_pt_first_derivatives

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, ct
real (r14), intent(out), optional :: pt_sa_sa, pt_sa_ct, pt_ct_ct

real (r14) :: ct_l, ct_u, dct, dsa, pt_ct_l, pt_ct_u, pt_sa_l, pt_sa_u, sa_l
real (r14) :: sa_u

if (present(pt_sa_sa)) then

    dsa  = 1d-3                 ! increment of absolute salinity is 0.001 g/kg
    sa_l = max(sa - dsa, 0d0)
    sa_u = sa + dsa

    call gsw_pt_first_derivatives(sa_l,ct,pt_sa=pt_sa_l)
    call gsw_pt_first_derivatives(sa_u,ct,pt_sa=pt_sa_u)

    pt_sa_sa = (pt_sa_u - pt_sa_l)/(sa_u - sa_l)

end if

if (present(pt_sa_ct) .or. present(pt_ct_ct)) then

    dct  = 1d-2     ! increment of conservative temperature is 0.01 degrees c
    ct_l = ct - dct
    ct_u = ct + dct

    if (present(pt_sa_ct) .and. present(pt_ct_ct)) then

        call gsw_pt_first_derivatives(sa,ct_l,pt_sa_l,pt_ct_l)
        call gsw_pt_first_derivatives(sa,ct_u,pt_sa_u,pt_ct_u)

    else if (present(pt_sa_ct)) then

        call gsw_pt_first_derivatives(sa,ct_l,pt_sa=pt_sa_l)
        call gsw_pt_first_derivatives(sa,ct_u,pt_sa=pt_sa_u)

    else if (present(pt_ct_ct)) then

        call gsw_pt_first_derivatives(sa,ct_l,pt_ct=pt_ct_l)
        call gsw_pt_first_derivatives(sa,ct_u,pt_ct=pt_ct_u)

    end if

    if (present(pt_sa_ct)) pt_sa_ct = (pt_sa_u - pt_sa_l)/(ct_u - ct_l)

    if (present(pt_ct_ct)) pt_ct_ct = (pt_ct_u - pt_ct_l)/(ct_u - ct_l)

end if

return
end subroutine

!--------------------------------------------------------------------------
