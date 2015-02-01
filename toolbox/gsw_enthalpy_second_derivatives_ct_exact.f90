!==========================================================================
elemental subroutine gsw_enthalpy_second_derivatives_ct_exact (sa, ct, p, &
                                                 h_sa_sa, h_sa_ct, h_ct_ct)
!==========================================================================
!
!  Calculates three second-order derivatives of specific enthalpy (h).
!
!  Note that this function uses the full Gibbs function.  There is an 
!  alternative to calling this function, namely 
!  gsw_enthalpy_second_derivatives which uses the computationally
!  efficient 48-term expression for density in terms of sa, ct and p 
!  (IOC et al., 2010).   
!
!  sa  =  Absolute Salinity                                        [ g/kg ]
!  ct  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  h_sa_sa  =  The second derivative of specific enthalpy with respect to 
!              Absolute Salinity at constant ct & p.    [ J/(kg (g/kg)^2) ]
!  h_sa_ct  =  The second derivative of specific enthalpy with respect to 
!              sa and ct at constant p.                  [ J/(kg K(g/kg)) ]
!  h_ct_ct  =  The second derivative of specific enthalpy with respect to 
!              ct at constant sa and p.                      [ J/(kg K^2) ]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_gibbs, gsw_pt_from_ct, gsw_pt_from_t

use gsw_mod_teos10_constants, only : gsw_cp0, gsw_t0

use gsw_mod_kinds

implicit none

real (r8), intent(in) :: sa, ct, p
real (r8), intent(out), optional :: h_sa_sa, h_sa_ct, h_ct_ct

real (r8) :: factor, gs_pt0, gst_pt0, gst_t, part, pt0, h_ct_ct_val
real (r8) :: rec_abs_pt0, rec_gtt_pt0, rec_gtt_t, t, temp_ratio

integer, parameter :: n0=0, n1=1, n2=2
real (r8), parameter :: pr0 = 0.0_r8, sa_small = 1e-100_r8

pt0 = gsw_pt_from_ct(sa,ct)
rec_abs_pt0 = 1.0_r8/(gsw_t0 + pt0)
t = gsw_pt_from_t(sa,pt0,pr0,p)
temp_ratio = (gsw_t0 + t)*rec_abs_pt0

rec_gtt_pt0 = 1.0_r8/gsw_gibbs(n0,n2,n0,sa,pt0,pr0)
rec_gtt_t = 1.0_r8/gsw_gibbs(n0,n2,n0,sa,t,p)
gst_pt0 = gsw_gibbs(n1,n1,n0,sa,pt0,pr0)
gst_t = gsw_gibbs(n1,n1,n0,sa,t,p)
gs_pt0 = gsw_gibbs(n1,n0,n0,sa,pt0,pr0)

! h_ct_ct is naturally well-behaved as sa approaches zero. 
h_ct_ct_val = gsw_cp0*gsw_cp0* &
    (temp_ratio*rec_gtt_pt0 - rec_gtt_t)*(rec_abs_pt0*rec_abs_pt0)

if (present(h_ct_ct)) h_ct_ct = h_ct_ct_val

part = (temp_ratio*gst_pt0*rec_gtt_pt0 - gst_t*rec_gtt_t)*rec_abs_pt0
factor = gs_pt0/gsw_cp0

if (present(h_sa_sa)) then

    ! h_sa_sa has a singularity at sa = 0, and blows up as sa approaches zero.  
    h_sa_sa = gsw_gibbs(n2,n0,n0,sa,t,p) &
        - temp_ratio*gsw_gibbs(n2,n0,n0,sa,pt0,pr0)  &
        + temp_ratio*gst_pt0*gst_pt0*rec_gtt_pt0  &
        - gst_t*gst_t*rec_gtt_t  &
        - 2.0_r8*gs_pt0*part + (factor*factor)*h_ct_ct_val

end if
if (.not. present(h_sa_ct)) return

! h_sa_ct should not blow up as sa approaches zero.  The following lines
! of code ensure that the h_sa_ct output of this function does not blow
! up in this limit.  That is, when sa < 1e-100 g/kg, we force the h_sa_ct 
! output to be the same as if sa = 1e-100 g/kg.  
if (sa .lt. sa_small) then
    rec_gtt_pt0 = 1.0_r8/gsw_gibbs(n0,n2,n0,sa_small,pt0,pr0)
    rec_gtt_t = 1.0_r8/gsw_gibbs(n0,n2,n0,sa_small,t,p)
    gst_pt0 = gsw_gibbs(n1,n1,n0,sa_small,pt0,pr0)
    gst_t = gsw_gibbs(n1,n1,n0,sa_small,t,p)
    gs_pt0 = gsw_gibbs(n1,n0,n0,sa_small,pt0,pr0)
    part = (temp_ratio*gst_pt0*rec_gtt_pt0 - gst_t*rec_gtt_t)*rec_abs_pt0
    factor = gs_pt0/gsw_cp0
end if

h_sa_ct  = gsw_cp0*part - factor*h_ct_ct_val

return
end subroutine

!--------------------------------------------------------------------------
