!==========================================================================
function gsw_sp_from_sk(sk)       
!==========================================================================

! Calculates Practical Salinity, SP, from SK
!
!  SK    : Knudsen Salinity                        [parts per thousand, ppt]
!
! gsw_sp_from_sk  : Practical Salinity                              [unitless]

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sk, gsw_sp_from_sk

gsw_sp_from_sk = (sk - 0.03d0)*(1.80655d0/1.805d0) 

! This line ensures that SP is non-negative.
if (gsw_sp_from_sk.lt.0d0) then
	gsw_sp_from_sk = 9d15
end if

return
end function

!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! salinity and temperature conversions
!--------------------------------------------------------------------------

