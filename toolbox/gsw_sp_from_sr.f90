!==========================================================================
function gsw_sp_from_sr(sr)  
!==========================================================================

! Calculates Practical Salinity, sp, from Reference Salinity, sr. 
!
! sr     : Reference Salinity                              [g/kg]
!
! gsw_sp_from_sr  : Practical Salinity                     [unitless]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sr, gsw_sp_from_sr

gsw_sp_from_sr = 0.995306702338459d0*sr;

if (gsw_sp_from_sr.gt.1d10) then
    gsw_sp_from_sr = 9d15
end if

return
end function

!--------------------------------------------------------------------------

