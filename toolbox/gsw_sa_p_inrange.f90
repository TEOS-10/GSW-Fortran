!==========================================================================
elemental function gsw_sa_p_inrange (sa, p)
!==========================================================================
!
!  Check for any values that are out of the TEOS-10 range ...
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar ) 
!---------------------------------------------------------------------------

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, p

logical :: gsw_sa_p_inrange

gsw_sa_p_inrange = .true.

if (p.gt.10000 .or. sa.gt.120 .or. &
    p + sa*71.428571428571402d0.gt.13571.42857142857d0) &
    gsw_sa_p_inrange = .false.

return
end function

!---------------------------------------------------------------------------
