!==========================================================================
function gsw_entropy_part_zerop(sa,pt0)
!==========================================================================

! entropy part evaluated at the sea surface
!
! sa     : Absolute Salinity                               [g/kg]
! pt0    : insitu temperature                              [deg C]
! 
! gsw_entropy_part_zerop : entropy part at the sea surface

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sa, pt0, sfac, x2, x, y, g03, g08, gsw_entropy_part_zerop

sfac = 0.0248826675584615d0

x2 = sfac*sa
x = sqrt(x2)
y = pt0*0.025d0

g03 = y*(-24715.571866078d0 + y*(2210.2236124548363d0 + &
    y*(-592.743745734632d0 + y*(290.12956292128547d0 + &
    y*(-113.90630790850321d0 + y*21.35571525415769d0)))))

g08 = x2*(x*(x*(y*(-137.1145018408982d0 + y*(148.10030845687618d0 + &
    y*(-68.5590309679152d0 + 12.4848504784754d0*y)))) + &
    y*(-86.1329351956084d0 + y*(-30.0682112585625d0 + y*3.50240264723578d0))) + &
    y*(1760.062705994408d0 + y*(-675.802947790203d0 + &
    y*(365.7041791005036d0 + y*(-108.30162043765552d0 + 12.78101825083098d0*y)))))

gsw_entropy_part_zerop = -(g03 + g08)*0.025d0

return
end function

!--------------------------------------------------------------------------

