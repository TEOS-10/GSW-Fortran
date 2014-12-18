!==========================================================================
elemental subroutine gsw_ct_second_derivatives (sa, pt, ct_sa_sa, ct_sa_pt, &
                                                ct_pt_pt)
!==========================================================================
!
!  Calculates the following three, second-order derivatives of Conservative 
!  Temperature
!   (1) CT_SA_SA, the second derivative with respect to Absolute Salinity  
!       at constant potential temperature (with p_ref = 0 dbar),
!   (2) CT_SA_pt, the derivative with respect to potential temperature
!       (the regular potential temperature which is referenced to 0 dbar)
!       and Absolute Salinity, and
!   (3) CT_pt_pt, the second derivative with respect to potential 
!       temperature (the regular potential temperature which is referenced 
!       to 0 dbar) at constant Absolute Salinity. 
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  pt  =  potential temperature (ITS-90)                          [ deg C ]   
!         (whose reference pressure is 0 dbar)
!
!  CT_SA_SA  =  The second derivative of Conservative Temperature with 
!               respect to Absolute Salinity at constant potential 
!               temperature (the regular potential temperature which 
!               has reference sea pressure of 0 dbar).  
!               CT_SA_SA has units of:                     [ K/((g/kg)^2) ]
!  CT_SA_pt  =  The derivative of Conservative Temperature with 
!               respect to potential temperature (the regular one with 
!               p_ref = 0 dbar) and Absolute Salinity.   
!               CT_SA_pt has units of:                        [ 1/(g/kg) ]
!  CT_pt_pt  =  The second derivative of Conservative Temperature with 
!               respect to potential temperature (the regular one with 
!               p_ref = 0 dbar) at constant SA.   
!               CT_pt_pt has units of:                              [ 1/K ]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_ct_first_derivatives

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, pt
real (r14), intent(out), optional :: ct_sa_sa, ct_sa_pt, ct_pt_pt

real (r14) :: ct_pt_l, ct_pt_u, ct_sa_l, ct_sa_u, dpt, dsa, pt_l, pt_u
real (r14) :: sa_l, sa_u

if (present(ct_sa_sa)) then

    dsa  = 1d-3        ! increment of absolute salinity is 0.001 g/kg.
    sa_l = max(sa - dsa, 0d0)
    sa_u = sa + dsa

    call gsw_ct_first_derivatives(sa_l,pt,ct_sa=ct_sa_l)
    call gsw_ct_first_derivatives(sa_u,pt,ct_sa=ct_sa_u)

    ct_sa_sa = (ct_sa_u - ct_sa_l)/(sa_u - sa_l)

end if

if (present(ct_sa_pt) .or. present(ct_pt_pt)) then

    dpt  = 1d-2        ! increment of potential temperature is 0.01 degrees c.
    pt_l = pt - dpt
    pt_u = pt + dpt

    if (present(ct_sa_pt) .and. present(ct_pt_pt)) then

        call gsw_ct_first_derivatives(sa,pt_l,ct_sa_l,ct_pt_l)
        call gsw_ct_first_derivatives(sa,pt_u,ct_sa_u,ct_pt_u)

    else if (present(ct_sa_pt) .and. .not. present(ct_pt_pt)) then

        call gsw_ct_first_derivatives(sa,pt_l,ct_sa=ct_sa_l)
        call gsw_ct_first_derivatives(sa,pt_u,ct_sa=ct_sa_u)

    else if (.not. present(ct_sa_pt) .and. present(ct_pt_pt)) then

        call gsw_ct_first_derivatives(sa,pt_l,ct_pt=ct_pt_l)
        call gsw_ct_first_derivatives(sa,pt_u,ct_pt=ct_pt_u)

    end if

    if (present(ct_sa_pt)) ct_sa_pt = (ct_sa_u - ct_sa_l)/(pt_u - pt_l)

    if (present(ct_pt_pt)) ct_pt_pt = (ct_pt_u - ct_pt_l)/(pt_u - pt_l)

end if

return
end subroutine

!--------------------------------------------------------------------------
